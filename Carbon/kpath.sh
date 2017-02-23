#!/bin/bash

# Dependencies: 
#  1. Key exchange between carbon and your local machine
#  2. ssh-daemon has to be running on Carbon with your key added
#  3. A folder on your local machine with kpath.py script located in it (e.g. ~/carbon_scratch)
#  4. Change the python executable location on line 16 to match that of the Python executable on your local machine
#  5. Instead of deadpool, have an entry in your ~/.ssh/config file for your own local machine

# Purpose: Goes to a folder on your local machine to run a script
# that uses pymatgen for writing KPOINTS files with high symmetry
# lines for band structure calculations. 

fold="~/carbon_scratch"
scp POSCAR deadpool:$fold/
ssh -A deadpool << EOF
    cd $fold
    /home/nicholas/CODE/git/aiida/bin/python kpath.py
EOF
scp deadpool:$fold/KPOINTS.bands ./
