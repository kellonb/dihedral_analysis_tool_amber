# dihedral_analysis_tool_amber

This page describes a simple dihedral analysis tool. This tool can be useful based on the following scenarios:

Scenario 1: You have two Xray structures of the beginning state and ending state of a pathway in PDB format and 
you will like to know what are the dihedral changes that occurs when going from one state to another. 
Additionally, you will like to know whether your trajectory exhibit the dihedral changes between the two states.
Scenario 2: You have two trajectories and you want to compare the differences between the dihedral changes
Scenario 3: You want to compare the change in dihedral of one trajectory to a reference structure in PDB format.

It has three python scripts. The diih_analysis.py is a python module I wrote that has all the functions, 
the analysis.examples.py is a script with some examples of how to use the functions in dih_analysis.py and 
the compare_traj.py is a script that plots the deviation for all 7 types of dihedral. If you want to test it 
right away you can use the compare_traj.py by making the changes below. 
 
The python script uses pytraj which is a python interface/API to cpptraj. you can download it here as well [https://amber-md.github.io/pytraj/latest/installation.html]. 
Also, this python script was tested and developed using python 3 which can be downloaded from here as well [https://www.anaconda.com/download/#linux]
 
 
