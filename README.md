# hwo-tools

This repo contains the science simulation tools for the Habitable Worlds Observatory. 

## Setting up a conda environment to use what is on the main branch
Try the following steps to set up the basic functionality using a supplied conda environment. From this environment you can execute the code examples contained in the notebooks directory.

> [!CAUTION]
> There are two repositories that are still private which the code committed to main currently depends on.
> You will not be able to fully clone all depenencies unless you have access to them.


The $INSTALL_DIR is the name of any local directory of your choice. After you go into
that directory, in order to create a conda environment with the specifications to use
what is on the main branch:

The following directions assume that you have ssh keys setup with GitHub.
```
cd $INSTALL_DIR
git clone git@github.com:spacetelescope/hwo-tools.git
conda env create -file hwotools.yml
conda activate hwotools
pip install git+ssh://git@github.com/spacetelescope/syotools.git
```
Make sure that your environment variables point to the appropriate place to find:

NOTE: If you already have $PYSYN_CDBS set in your environment, you *do NOT need* to set the environment
variable to the subset that comes with SYOTools - a full checkout will work just as well.
```
export PYSYN_CDBS= "the full path to where syotools/reference_data/pysynphot_data lives in your env"
export SYOTOOLS_DATA_DIR="full path to syotools/reference_data in your env" 
export HWOME_DATA_PATH=$INSTALL_DIR/hwome_data
export YIP_CORO_DIR="$INSTALL_DIR/yip"
```
To test that everything is working as expected, open and run BasicRun.ipynb

If this gives you S/N values, you have it working.

## Setting yourself up for local development
- Clone the hwo-tools repo:
```git clone https://github.com/spacetelescope/hwo-tools.git```

For Imaging and Spectroscopy (camera_etc, uvspec_etc, ifs_etc):
- Clone the hwome-core repo:  **(This repo is NOT currently public)**
```
   cd $INSTALL_DIR
   git clone git@github.com:HWO-Project/hwome-core.git
   cd hwome-core
   pip install .
```
- Clone hwome_data repo: **(This repo is NOT currently public)**
```   cd $INSTALL_DIR
   git clone git@github.com:HWO-Project/hwome_data.git
   set environment variable: export HWOME_DATA_PATH=$INSTALL_DIR/hwome_data
```
- Clone the SYOTools repo:
```
   cd $INSTALL_DIR
   git clone git@github.com:spacetelescope/syotools.git
   cd syotools
   pip install .
```

And for Coronagraphy (coron_imaging, coron_spec):
- Clone the EACy repo:
```
   cd $INSTALL_DIR
   git clone git@github.com:curriem/eacy.git
   cd eacy
   pip install .
```

- For the Coronagraphic Spectroscopy and Imaging ETCs, you will also need a YIP file for the coronagraph; put it in $INSTALL_DIR/yip

- Clone the pyEDITH repo:
```
   cd $INSTALL_DIR
   git clone git@github.com:HabitableWorldsObservatory/pyEDITH.git
   cd pyEDITH
   pip install .
```

The above commands should have pulled in all necessary dependencies except bokeh and jupyter:
```  pip install bokeh jupyter```


More specific, complete environments for conda-forge users.
   - The hwotools_linuxx86-64.yml files are for an 64-bit Intel Linux computer, and come with MKL-accelerated numpy and scipy.
   - The hwotools_macarm.yml files are for Apple Silicon computers, and come with Apple Accelerate-accelerated numpy and scipy.

Install SYOTools (will pull in hwome-core):
```
pip install git+ssh://git@github.com/spacetelescope/syotools.git```
conda activate hwotools
```
Add to your .bashrc / .bash_profile:
```
export PYSYN_CDBS=/Users/tumlinson/anaconda3/envs/hwotools/lib/python3.12/site-packages/syotools/reference_data/pysynphot_data
export SYOTOOLS_DATA_DIR=/Users/tumlinson/anaconda3/envs/hwotools/lib/python3.12/site-packages/syotools/reference_data
export HWOME_DATA_PATH=$INSTALL_DIR/hwome_data
export YIP_CORO_DIR="$INSTALL_DIR/yip"
```
These exact pathnames will vary by system, please do your best.

Your root directory (here, /Users/tumlinson/anaconda3/) may vary, and your python version may as well. 
In the likely event that you already have PYTHONPATH set in your .*rc file, append these:
```export PYTHONPATH=$PYTHONPATH:$INSTALL_DIR/hwo-tools/```

Open and run BasicRun.ipynb. If this gives you S/N values, you have it working.

## Using the tools and example notebooks
Once you have this basic test working, you can run the flat python script 
BasicRun.py or create your own. 

The camera_wrapper and uvspec_wrapper notebooks show the simplest way to call these tools, using bare python wrappers around the SYOTools API. These allow you to get SNR results with one import and one line of code. Use these if you need to run in 'batch mode', which the online GUI tools will not do. 

For a deeper illustration of how the tools work, try one of the notebooks in the notebooks directory, like Camera_ETC_Tutorial and UVSpec_ETC_Tutorial.  


## Migrating from SYOTools 1.3 and earlier to SYOTools 1.4 (HWOME version)
### 1. You need to set the "unknown" last, now.

Previous versions of SYOTools had defaults sourced from LUVOIR-era telescopes and had the necessary information to do calculations immediately. With HWOME and its more detailed EACs, all the necessary properties must be set before any calculations can be done.

Previous versions of SYOTools defaulted to computing for SNR and would start to do so immediately, recalculating any time a parameter changed. The new SYOTools WILL recompute on all changes, but only after the SourceExposure.unknown is set to "snr", "exptime", or "magnitude" for the first time.

This is the new order of operations:
* Create a Telescope() 
* Load an EAC into the Telescope with Telescope.set_from_hwome()
* Select an instrument from the Telescope's list
  * Manually
  * Runtime discovery using Telescope.find_instrument_with()
* Create the correct class of SourceExposure() (or use Instrument.create_exposure())
* Create a Source()
* Set the source SED (optional, there is still a default flat source)
* Add the Source to the SourceExposure
* Add the SourceExposure to the Instrument
* Set exptime and/or snr in the SourceExposure
* Set the unknown to be solved for, in the SourceExposure. **(must be done last)**

### 2. You do not create a Camera, Spectrograph, or IFS any more. 
Instead, you select the instrument from the Telescope.instruments dict.

The function that populates the Telescope with an EAC definition from HWOME *populates the entire EAC for you*, and you simply select the already-loaded instrument from the Telescope.instruments list.

This new behavior supports SCDDs that use more than one instrument, and running all SCDDs through a single telescope configuration for DISRA uses.

### 3. There are more than two instruments now.
Instrument now means Instrument Channel, and there are now many, entirely defined by data from HWOME.

Previous versions of SYOTools had one Camera, "HRI"; one Spectograph, "UVI"; and one IFU, "IFS". Now (for instance) there are four data-defined cameras: HRI_S.HRI_S_UVIS, HRI_S.HRI_S_NIR, HRI_A.HRI_A_VIS_Imager, and UV_MOS.FUV_IMG_Imager channels, each with their own filters (which themselves no longer have simple names like "V" and "R"), detectors, and so on. All are now considered Instruments; Telescope no longer has separate attributes to hold a single camera, spectrograph, or ifs.

It is recommended to use functions like telescope.find_instrument_with() to *discover* the filter bandpass you need, rather than rely on specific names; the entire system is data-driven and in flux.


Some more minor notes:
* The outputs of calculate(), calculate_snr(), calculate_exptime(), and calculate_magnitude() are always lists now, even if they only have one element
* The outputs of imaging, spectroscopic, and IFU calculations can now be just a single filter/disperser OR all of them.
  * In previous versions of SYOTools, camera calculations ran every filter; spectroscopy and ifs calculations ran one selected spectroscopic band. In the new version, all instrument types can do either.
  * To run just one filter (or spectroscopic bandpass): set instrument.band or pass a `custom_band=<name of band>` keyword argument into calculate(), calculate_snr(), calculate_exptime(), or calculate_magnitude(). 
  * To run all filters or all spectroscopic bandpasses, either don't set instrument.band, set instrument.band to `None`, or pass `custom_band=None` to calculate(), calculate_snr(), etc.
* SNR is now set with an attribute just called "snr". Previous SYOTools used the name "snr_goal" as the input SNR name.

