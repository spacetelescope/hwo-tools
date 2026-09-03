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
git clone git+ssh://git@github.com:spacetelescope/hwo-tools.git
conda create --file hwotools.yml
conda activate hwotools
pip install git_ssh://git@github.com:spacetelescope/syotools.git
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
   git clone https://github.com/HWO-Project/hwome-core.git
   cd hwome-core
   pip install .
```
- Clone hwome_data repo: **(This repo is NOT currently public)**
```   cd $INSTALL_DIR
   git clone https://github.com/HWO-Project/hwome_data.git
   set environment variable: export HWOME_DATA_PATH=$INSTALL_DIR/hwome_data
```
- Clone the SYOTools repo:
```
   cd $INSTALL_DIR
   git clone https://github.com/spacetelescope/syotools.git
   cd syotools
   pip install .
```

And for Coronagraphy (coron_imaging, coron_spec):
- Clone the EACy repo:
```
   cd $INSTALL_DIR
   git clone https://github.com/curriem/eacy.git
   cd eacy
   pip install .
```

- For the Coronagraphic Spectroscopy and Imaging ETCs, you will also need a YIP file for the coronagraph; put it in $INSTALL_DIR/yip

- Clone the pyEDITH repo:
```
   cd $INSTALL_DIR
   git clone https://github.com/HabitableWorldsObservatory/pyEDITH.git
   cd pyEDITH
   pip install .
```

The above commands should have pulled in all necessary dependencies except bokeh and jupyter:
```  pip install bokeh jupyter```


More specific, complete environments
   - The hwotools_linuxx86-64.yml files are for an 64-bit Intel Linux computer, and come with MKL-accelerated numpy and scipy.
   - The hwotools_macarm.yml files are for Apple Silicon computers, and come with Apple Accelerate-accelerated numpy and scipy.

Install SYOTools (will pull in hwome-core):
```
pip install git+ssh://git@www.github.com/spacetelescope/syotools```
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
