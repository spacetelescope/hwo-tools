# hwo-tools

This repo contains the science simulation tools for the Habitable Worlds Observatory. 

Try these steps to set up the basic functionality. 

- Clone the hwo-tools repo:
   git clone https://github.com/spacetelescope/hwo-tools.git
- Clone the Sci-Eng-Interface repo:
   git clone https://github.com/HWO-GOMAP-Working-Groups/Sci-Eng-Interface.git
- Set up the conda repo with the hwotools dependencies: 
      cd hwo-tools; conda env create -f hwotools.yml
- conda activate hwotools
- add to your .bashrc / .bash_profile:
      export PYTHONPATH=~/Sci-Eng-Interface/
      export SCI_ENG_DIR=~/Sci-Eng-Interface/
- open and run BasicRun.ipynb
- If this gives you S/N values, you have it working.
