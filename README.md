# hwo-tools

This repo contains the science simulation tools for the Habitable Worlds Observatory. 

Try these steps to set up the basic functionality using a supplied conda environment. 
The $INSTALL_DIR is just any local directory of your choice. 

- cd $INSTALL_DIR
- Clone the hwo-tools repo:
   git clone https://github.com/spacetelescope/hwo-tools.git
- Clone the Sci-Eng-Interface repo:
   git clone https://github.com/HWO-GOMAP-Working-Groups/Sci-Eng-Interface.git
- Set up the conda repo with the hwotools dependencies: \
    cd $INSTALL_DIR/hwo-tools; conda env create -f hwotools.yml
- conda activate hwotools 
- add to your .bashrc / .bash_profile:
   - export PYTHONPATH=$INSTALL_DIR/Sci-Eng-Interface/:$INSTALL_DIR/hwo-tools/
   - export SCI_ENG_DIR=$INSTALL_DIR/Sci-Eng-Interface/ 
   - export PYSYN_CDBS=$INSTALL_DIR/hwo-tools/pysynphot_data/grp/redcat/trds/
   - in the likely event that you already have PYTHONPATH set in your .*rc file, append these:
        - export PYTHONPATH=$PYTHONPATH:$INSTALL_DIR/Sci-Eng-Interface/:$INSTALL_DIR/hwo-tools/
- open and run BasicRun.ipynb
- If this gives you S/N values, you have it working.

Once you have this basic test working, you can run the flat python script 
BasicRun.py or create your own. 
