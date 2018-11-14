import pytraj as pt
import dih_analysis as dap 
import os 
import sys
import pandas as pd


# parameters 
tot_res = 10 # total residues 

# first create an instance of traj
top = "1hhp.18SB.strip.parm7" # path to your topology file
crd = "9md.strip.nc"  # path to your trajectory file 
traj = pt.iterload([crd], top=top)
ref_pdb = 'PDB/1hhp.pdb' # path to reference pdb

# plottin parameters
residue_range = [2,8] # 2 is residue 2 and 8 is residue 8
traj_range = [1, 50] # range of trajectory you want to look at
# range of deviations you want to look at
vmin = 0
vmax = 40 
# dictionary to hold all values for a given dihedral type in the trajectory
dev = {}
dihtype = {}
# list of dihedrals 
dih_list = ['phi', 'psi', 'chi1', 'chi2', 'chi3', 'chi4', 'chi5']
print ("DOING compare_traj2PDB")
for dih in dih_list:
  dev[dih], dihtype[dih] = dap.compare_traj2PDB(traj, ref_pdb, tot_res, dih, residue_range)
  print (dev['phi']) # an example
  # Example 12: plot a heatmap of the deviation of a dihedral
  dap.plot_dev_heatmap(dev[dih], vmin, vmax, "%s_dev" %(dih), traj_range)

#dap.plot_dev_heatmap(dev_df, vmin, vmax, fig_title,  x_lim)
#example of plotting the deviation in heatmap



