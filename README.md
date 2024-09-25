# hwo-tools

This repo contains the science simulation tools for the Habitable Worlds Observatory. 

Try these steps to set up the basic functionality using a supplied conda environment. 
The $INSTALL_DIR is just any local directory of your choice. 

- cd $INSTALL_DIR
- Clone the hwo-tools repo:
   git clone https://github.com/spacetelescope/hwo-tools.git
- The Sci-Eng-Interface material is now bundled with SYOtools 
- Set up the conda repo with the hwotools dependencies: \
    cd $INSTALL_DIR/hwo-tools; conda env create -f hwotools.yml
- conda activate hwotools
- add to your .bashrc / .bash_profile:
   - export PYTHONPATH=$INSTALL_DIR/Sci-Eng-Interface/:$INSTALL_DIR/hwo-tools/
   - export PYSYN_CDBS=/Users/tumlinson/anaconda3/envs/hwotools/lib/python3.12/site-packages/syotools/pysynphot_data
   - in the likely event that you already have PYTHONPATH set in your .*rc file, append these:
        - export PYTHONPATH=$PYTHONPATH:$INSTALL_DIR/hwo-tools/
- open and run BasicRun.ipynb
- If this gives you S/N values, you have it working.

Once you have this basic test working, you can run the flat python script 
BasicRun.py or create your own. 

The camera_wrapper and uvspec_wrapper notebooks show the simplest way to call these tools, 
using bare python wrappers around the SWOTools API. These allow you to get SNR results
with one import and one line of code. Use these if you need to run in 'batch mode', which the 
online GUI tools will not do. 

For a deeper illustration of how the tools work, try the other notebooks 
Camera_ETC_Tutorial and UVSpec_ETC_Tutorial.  
