
import astropy.units as u 
import numpy as np 


def calc_count_rates_imaging(int_time, input_aperture):
    epsilon = 1.0e-10 # mean astrophysical planet-to-star flux ratio over the bandpass 
    quantum_eff = 0.9 
    raw_contrast = 1e-10 
    eta_p = 0.2 
    sigma_DeltaC = 1e-11 
    telescope_aperture = input_aperture * u.m 

    # Sun at 10 pc in V band (V_abs_Sun = 4.8) in units photons cm^-2 s^-1 micron^-1 - using https://irsa.ipac.caltech.edu/data/SPITZER/docs/dataanalysistools/tools/pet/magtojy/ 
    sun_flux = 4.51e-7 * u.erg / u.cm**2 / u.s / u.micron
    ph_energy = 6.626196e-27 * u.erg * u.s * (2.9979e10 * u.cm / u.s) / (5500 * u.Angstrom).to(u.cm) / u.ph 
    sun_phot = sun_flux / ph_energy 
 
    phi_S = sun_phot / 12.**2 * 10.**2 # correct photon rate to 12 pc for a Solar twin - should be photons / s/ cm^2 / micron 

    telescope_aperture = 8. * u.m 
    collecting_area = np.pi * (telescope_aperture.to(u.cm) / 2.)**2 # square cm, 6 meters, filled aperture, this will have to come in from the telescope object in the end. 

    dlambda = 0.55 / 5 * u.micron # bandpass taken from 0.55 micron central wavelength and spectral resolution per Table 2 
    
    Tr = 0.3 # Total throughput 
    N_s = phi_S * collecting_area * quantum_eff * Tr * dlambda # total number of stellar photons per sec in the absence of any corongraph - Eq (8) 

    # the various count rates 
    r_pl = eta_p * epsilon * N_s # planet photon rate (ph / s) detected in a circular photometric aperture Eq (3) 
    r_sp = eta_p * raw_contrast * N_s      # photon rate from residual starlight in speckles Eq (4) 
    r_sz = 1. * u.ph / u.s / u.micron * dlambda * Tr # read off the upper right panel of Figure 23. 
    r_xz = 7. * u.ph / u.s / u.micron * dlambda * Tr # read off the upper right panel of Figure 23. 

    r_n = r_pl + r_sp + r_sz + r_xz # total noise rate, sum of terms just below Eq (2) 
    post_term = eta_p * sigma_DeltaC * N_s * int_time.to(u.s) / u.ph #<---- see Eq 13

    # the factor of 2 on r_n here in the denominator is from Eq 13 and accounts for differential imaging 
    signal_to_noise = (r_pl * int_time.to(u.s) / u.ph) / (2. * r_n * int_time.to(u.s)/u.ph  + post_term**2 )**0.5 

    return signal_to_noise


print(calc_count_rates_imaging(2. * u.hr, 8. * u.m)) 
