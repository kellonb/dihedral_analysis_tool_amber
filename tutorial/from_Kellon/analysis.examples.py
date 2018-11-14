import pytraj as pt
import dih_analysis as dap 
import os 
import sys
import pandas as pd


# parameters 
tot_res = 10 # total residues 


# Example 1: download PDB from a list of pdb codes 
# http://simmerlinglab.org/wiki/index.php?title=Dihedral_analysis_tool&action=edit&section=3
# use sys to open pdb_list and save it to an array name pdb_list
# use the download_pdb function to download the pdbs and save them in PDB dirctory
#dap.download_pdb("pdb_list")


# Example 2: how to get a atom mask for a dihedral. limited options 
dih_mask = dap.dihmask(1, 'chi1', 'ARG')
print (dih_mask)

# Example 3: calculate a dihdral across a trajectory for a given mask
# first create an instance of traj
top = "1hhp.18SB.strip.parm7" # your topology file
crd = "9md.strip.nc"  # trajectory file 
traj = pt.iterload([crd], top=top)
# zero dih_mask and create new mask for R8 chi1
dih_mask = 0 
dih_mask = dap.dihmask(8, 'chi1', 'ARG')
R8_chi1 = dap.cal_dih(traj, dih_mask) 
print (R8_chi1)

# Example 4: create a library of dihedral values for Xray structures
pdb = dap.create_PDBdf("./pdb_list", tot_res)
#if I want an individual data frame, for example 1hhp
print ("1hhp")
print (pdb['1hhp'])
print ("2qnp")
print (pdb['2qnp'])

# Example 5: take two Xray structures from the library pdb and creat a dataframe with the#            difference in their dihedral value
diff_df = dap.compare_2PDB(pdb['1hhp'], pdb['2qnp'], tot_res, "1hhp_2qnp_diff")
print ("difference between 1hhp and 2qnp")
print (diff_df.to_string())

# Example 6: how to get the atom mask for all dihedrals that differ by a given threshold
threshold = 30 
diff_atom_mask, res_diff = dap.PDB_diff_mask(diff_df, threshold, tot_res, "1hhp_2qnp_t30")
print (diff_atom_mask)
print (res_diff)

# Example 7: plotting a scatter plot of dihedrals based on atom mask, using one reference 
# we will plot and let the function take care of the axis
ref_pdb = 'PDB/1hhp.pdb'
min_array = [0, 20] # min values for x and y axis
max_array = [10, 180] # max values for x and y axis
dap.plot_dih(top, crd, ref_pdb, diff_atom_mask, res_diff, "ff14SB")
dap.plot_dih(top, crd, ref_pdb, diff_atom_mask, res_diff, "ff14SB_limit", min_array, max_array)

# Example 8: plotting a scatter plot of dihedrals based on atom mask, using two references 
ref_pdb2 = 'PDB/2qnp.pdb'
dap.plot_dih2ref(top, crd, ref_pdb, ref_pdb2, diff_atom_mask, res_diff, "ff14SB_2ref")

# Example 9: get all the dihedral mask in a pdb
all_mask, ref_mask = dap.get_allmasks(pdb['1hhp'], tot_res)
print (all_mask)

# Example 10: compare two trajectories dihedral values and plot only the dihedral that differ 
thresh = 0.9
plot_thresh = 0.30
#define new traj 
top2 = "2qnp.18SB.strip.parm7" # your topology file
crd2 = "9md.2qnp.strip.nc"  # trajectory file 
traj2 = pt.iterload([crd2], top=top2)
label1 = '1hhp' ; label2 = '2qnp'
dap.plot_scat_2traj(traj, traj2, ref_pdb, all_mask, ref_mask, thresh, plot_thresh, label1, label2, "1hhpv2qnp")

# Example 11: get the deviations of a dihedral in a trajectory from a reference structure. 
# dictionary to hold all values for a given dihedral type in the trajectory
dev = {}
dihtype = {}
dih_list = ['phi', 'psi', 'chi1', 'chi2', 'chi3', 'chi4', 'chi5']
print ("DOING compare_traj2PDB")
for dih in dih_list:
  dev[dih], dihtype[dih] = dap.compare_traj2PDB(traj, ref_pdb, 10, dih, [2,8])
print (dev['phi'])

# Example 12: plot a heatmap of the deviation of a dihedral
#dap.plot_dev_heatmap(dev_df, vmin, vmax, fig_title,  x_lim)
#example of plotting the deviation in heatmap
dap.plot_dev_heatmap(dev['phi'], 0, 40, "phi_dev", [1,50])



