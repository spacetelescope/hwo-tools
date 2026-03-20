import os, time, shutil, datetime
import numpy as np
import yaml 
from astropy.io import fits
import astropy.units as u
import EXOSIMS.MissionSim as ems
from copy import deepcopy


class ExosimsWrapper:
    """
    Takes in a config dict specifying desired targets and their corresponding earth equivalent
    insolation distances (eeid), Earth-equivalent planet-star flux ratios (eepsr) and exo-zodi levels.
    Main function is run_exosims which given an input json file returns the necessary integration times and various
    output count rates.

    Parameters
    ----------
    config: dict
        dictionary of configuration parameters.
    """
    def __init__(self, config):
        # pull relevant values from the config
        self.working_angles = config["working_angles"]
        hip_numbers = [str(config['targets'][star]['HIP']) for star in config['targets']]
        self.target_list = [f"HIP {n}" for n in hip_numbers]
        self.eeid = [config['targets'][star]['eeid'] for star in config['targets']]
        self.eepsr = [config['targets'][star]['eepsr'] for star in config['targets']]
        self.exo_zodi = [config['targets'][star]['exo_zodi'] for star in config['targets']]

        # intitialize output arrays
        self.C_p = np.empty((len(self.target_list), len(self.working_angles)))
        self.C_b = np.empty((len(self.target_list), len(self.working_angles)))
        self.C_sp = np.empty((len(self.target_list), len(self.working_angles)))
        self.C_sr = np.empty((len(self.target_list), len(self.working_angles)))
        self.C_z = np.empty((len(self.target_list), len(self.working_angles)))
        self.C_ez = np.empty((len(self.target_list), len(self.working_angles)))
        self.C_dc = np.empty((len(self.target_list), len(self.working_angles)))
        self.C_rn = np.empty((len(self.target_list), len(self.working_angles)))
        self.C_star = np.empty((len(self.target_list), len(self.working_angles)))
        self.int_time = np.empty((len(self.target_list), len(self.working_angles))) * u.d

    def run_exosims(self):
        """
        Run EXOSIMS to generate results, including exposure times 
        required for reaching specified SNR.  

        """

        wa_coefs = self.config["working_angles"]
        n_angles = len(wa_coefs)
        target_list = self.target_list
        eeid = self.eeid
        eepsr = self.eepsr
        exo_zodi = self.exo_zodi
        sim = ems.MissionSim(use_core_thruput_for_ez=False
                             , **deepcopy(self.exosims_pars_dict))
        
        sInds = np.array([np.where(sim.TargetList.Name == t)[0][0] for t 
                         in target_list])
        
        # assemble information needed for integration time calculation:
        
        # we have only one observing mode defined, so use that
        mode = sim.OpticalSystem.observingModes[0]
        
        # use the nominal local zodi and exozodi values
        fZ = sim.ZodiacalLight.fZ0
        
        # now we loop through the targets of interest and compute integration 
        # times for each:
        int_time = np.empty((len(target_list), n_angles))*u.d
        for j, sInd in enumerate(sInds):
            # choose angular separation for coronagraph performance
            # this doesn't matter for a flat contrast/throughput, but
            # matters a lot when you have real performane curves
            # target planet deltaMag (evaluate for a range):
            WA = np.array(wa_coefs)*eeid[j]
            dMags = 5.0*np.log10(np.array(wa_coefs)) - 2.5*np.log10(eepsr[j])
            self.int_time[j] = sim.OpticalSystem.calc_intTime(
                sim.TargetList,
                [sInd] * n_angles,
                [fZ.value] * n_angles * fZ.unit,
                [exo_zodi[j]*sim.ZodiacalLight.fEZ0.value] * n_angles 
                    * sim.ZodiacalLight.fEZ0.unit,
                dMags,
                WA * u.arcsec,
                mode
            )
            counts = sim.OpticalSystem.Cp_Cb_Csp(
                sim.TargetList,
                [sInd] * n_angles,
                [fZ.value] * n_angles * fZ.unit,
                [exo_zodi[j]*sim.ZodiacalLight.fEZ0.value] * n_angles 
                    * sim.ZodiacalLight.fEZ0.unit,
                dMags,
                WA * u.arcsec,
                mode,
                True
            )

            self.C_p[j] = counts[0].value
            self.C_b[j] = counts[1].value
            self.C_sp[j] = counts[2].value
            self.C_sr[j] = counts[3]["C_sr"].value
            self.C_z = counts[3]["C_z"].value
            self.C_ez = counts[3]["C_ez"].value
            self.C_dc= counts[3]["C_dc"].value
            self.C_rn = counts[3]["C_rn"].value
            self.C_star = counts[3]["C_star"].value

        return self.int_time, self.C_p, self.C_b, self.C_sp, self.C_sr, self.C_z, self.C_ez, self.C_dc, self.C_rn, self.C_star


class ErrorBudget(ExosimsWrapper):
    """
    Exposure time calculator incorporating dynamical wavefront errors and
    WFS&C. Can also implement the Markov-chain-Monte-Carlo exploration of coronagraphic parameters.

    Parameters
    ----------
    config_file: str
        Name of the YAML configuration file.
    """
    def __init__(self, config_file):
        self.config_file = config_file
        with open(config_file, 'r') as config:
            self.config = yaml.load(config, Loader=yaml.FullLoader)
        super().__init__(self.config)

        self.input_dir = self.config["paths"]["input"]
        self.temp_dir = self.config["paths"]["temporary"]

        self.exosims_pars_dict = None


    def initialize_for_exosims(self):
        """Intitializes the EXOSIMS parameter dict with values from the config."""

        config = self.config

        self.exosims_pars_dict = config['initial_exosims']

        self.exosims_pars_dict['cherryPickStars'] = self.target_list

        for key in self.exosims_pars_dict['scienceInstruments'][0].keys():
            if key != 'optics' in dir(self):
                setattr(self, key, self.exosims_pars_dict['scienceInstruments'][0][key])
        for key in self.exosims_pars_dict['starlightSuppressionSystems'][0].keys():
            if key in dir(self):
                setattr(self, key, self.exosims_pars_dict['starlightSuppressionSystems'][0][key])


    def update_attributes(self, subsystem=None, name=None, value=None):
        """Updates the EXOSIMS parameter dict for the subsystem and name with the value.

        For example if subsystem = scienceInstruments and name = QE then
        self.exosims_pars_dict["scienceInstruments"]["QE] will take on the value.

        Parameters
        ----------
        subsystem: str
            Name of subsystem to update.
        name: str
            Name of the variable to update.
        value: float or int or str or bool
            Value to update the variable to.
        """
        if name is not None:
            try:
                self.exosims_pars_dict[subsystem][0][name] = value
            except KeyError:
                self.exosims_pars_dict[name][0] = value

    
    def run(self, subsystem=None, name=None, value=None, clean_files=True):
        """Main method for running EBS not in MCMC mode.

        Updates the variable defined by the subsystem and name with the given value before
        generating all the necessary files to run EXOSIMS as specified by the config and input
        CSV files. If no subsystem or name is given then EXOSIMS is just run with the input CSV files and
        config.

        Parameters
        ----------
        subsystem: str
            Name of subsystem to update.
        name: str
            Name of the variable to update.
        value: float or int or str or bool
            Value to update the variable to.
        clean_files: bool
            If True will remove temporary intermediate files after they are used.
        """
        self.initialize_for_exosims()
        self.update_attributes(subsystem=subsystem, name=name, value=value)

        self.run_exosims()

