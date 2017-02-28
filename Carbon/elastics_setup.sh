#!/bin/bash

# Dependencies: 
#  1. Key exchange between carbon and your local machine
#  2. ssh-daemon has to be running on Carbon with your key added
#  3. A folder on your local machine with setup_elastic_POSCARs.py script located in it (e.g. ~/carbon_scratch)
#  4. Change the python executable location on line 16 to match that of the Python executable on your local machine
#  5. Instead of deadpool, have an entry in your ~/.ssh/config file for your own local machine
#  6. cjob.py has to be located in your path on Carbon

# Purpose: Creates strained POSCAR files for extracting elastic constant matrix in line with Materials Project workflow.
# Then uses those POSCARs with provided INCAR.static to setup and queue 24 calculations. 

# See also extract_elast_consts.py for post-processing


folder="~/carbon_scratch"
scp POSCAR deadpool:$folder/
ssh -A deadpool << EOF
    cd $folder
    /home/nicholas/CODE/git/aiida/bin/python setup_elastic_POSCARs.py
EOF
scp deadpool:$folder/POSCAR* ./
for i in {1..24}
do
    mkdir poscar$i
    mv POSCAR$i poscar$i/POSCAR
    cp POTCAR KPOINTS.static INCAR.static poscar$i/
    cd poscar$i
    cjob.py -jt elast -n 3_elast${i}_13
    qsub r*.job > JID
    cd ..
done
