#!/bin/bash

# Dependencies: 
#  1. Key exchange between carbon and your local machine
#  2. ssh-daemon has to be running on Carbon with your key added
#  3. A folder on your local machine with calc_elasts.py script located in it (e.g. ~/carbon_scratch)
#  4. Change the python executable location on line 16 to match that of the Python executable on your local machine
#  5. Instead of deadpool, have an entry in your ~/.ssh/config file for your own local machine

# Purpose: Extracts matrix of elastic constants in Voigt notation as well as Voigt average bulk modulus
#          from the output of 24 calculations setup by elast_setup

folder="~/carbon_scratch"
scp -r POSCAR POTCAR poscar* deadpool:$folder/
ssh -A deadpool << EOF
    cd $folder
    /home/nicholas/CODE/git/aiida/bin/python calc_elasts.py
EOF
scp deadpool:$folder/elasts.txt ./
