# hwo-tools

This repo contains the science simulation tools for the Habitable Worlds Observatory. 

Try these steps to set up the basic functionality using a supplied conda environment. 
The $INSTALL_DIR is just any local directory of your choice. 

- cd $INSTALL_DIR
- Clone the hwo-tools repo:
   git clone https://github.com/spacetelescope/hwo-tools.git

For Imaging and Spectroscopy (camera_etc, uvspec_etc, ifs_etc):
- Clone the hwome-core repo:  **(This repo is NOT currently public)**
   cd $INSTALL_DIR
   git clone https://github.com/HWO-Project/hwome-core.git
   cd hwome-core
   pip install .
- Clone hwome_data repo: **(This repo is NOT currently public)**
   cd $INSTALL_DIR
   git clone https://github.com/HWO-Project/hwome_data.git
   set environment variable: export HWOME_DATA_PATH=$INSTALL_DIR/hwome_data
- Clone the SYOTools repo:
   cd $INSTALL_DIR
   git clone https://github.com/spacetelescope/syotools.git
   cd syotools
   pip install .

And for Coronagraphy (coron_imaging, coron_spec):
- Clone the EACy repo:
   cd $INSTALL_DIR
   git clone https://github.com/curriem/eacy.git
   cd eacy
   pip install .
- For the Coronagraphic Spectroscopy and Imaging ETCs, you will also need a YIP file for the coronagraph; put it in $INSTALL_DIR/yip
- Clone the pyEDITH repo:
   cd $INSTALL_DIR
   git clone https://github.com/HabitableWorldsObservatory/pyEDITH.git
   cd pyEDITH
   pip install .
   
The above commands should have pulled in all necessary dependencies except bokeh and jupyter:
   pip install bokeh jupyter

- Alternatively, set up a conda-forge installation and then run:
    cd $INSTALL_DIR/hwo-tools; conda env create -f hwotools.yml
- More specific, complete environments
   - The hwotools_linuxx86-64.yml files are for an 64-bit Intel Linux computer, and come with MKL-accelerated numpy and scipy.
   - The hwotools_macarm.yml files are for Apple Silicon computers, and come with Apple Accelerate-accelerated numpy and scipy. 


- conda activate hwotools
- add to your .bashrc / .bash_profile:
   - export PYSYN_CDBS=/Users/tumlinson/anaconda3/envs/hwotools/lib/python3.12/site-packages/syotools/reference_data/pysynphot_data
   - export SYOTOOLS_DATA_DIR=/Users/tumlinson/anaconda3/envs/hwotools/lib/python3.12/site-packages/syotools/reference_data
   - export HWOME_DATA_PATH=$INSTALL_DIR/hwome_data
   - export YIP_CORO_DIR="$INSTALL_DIR/yip"
   - These exact pathnames will vary by system, please do your best.
   - Your root directory (here, Users/tumlinson/anaconda3/) may vary, and your python version may as well. 
   - in the likely event that you already have PYTHONPATH set in your .*rc file, append these:
        - export PYTHONPATH=$PYTHONPATH:$INSTALL_DIR/hwo-tools/
- open and run BasicRun.ipynb
- If this gives you S/N values, you have it working.

Once you have this basic test working, you can run the flat python script 
BasicRun.py or create your own. 

The camera_wrapper and uvspec_wrapper notebooks show the simplest way to call these tools, 
using bare python wrappers around the SYOTools API. These allow you to get SNR results
with one import and one line of code. Use these if you need to run in 'batch mode', which the 
online GUI tools will not do. 

For a deeper illustration of how the tools work, try the other notebooks 
Camera_ETC_Tutorial and UVSpec_ETC_Tutorial.  
