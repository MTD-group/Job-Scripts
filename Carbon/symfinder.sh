#!/bin/bash
# Requires: ssh-agent running with key
# symmetry_finder.py on local machine

# Usage: uses symmetry_finder.py to get symmetrized cif from POSCAR
fold="~/carbon_scratch"
scp POSCAR POTCAR deadpool:$fold/ # POTCAR sent because a POTCAR's elements will overwrite POSCAR's in cif if present
ssh deadpool << EOF
 cd $fold
 /home/nicholas/CODE/git/aiida/bin/python ./symmetry_finder.py
EOF
name=$(head -n 1 POSCAR | awk '{print $1;}')
scp deadpool:$fold/$name.cif ./
