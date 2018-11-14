import pytraj as pt
import numpy as np
import pandas as pd
import csv
import sys
import seaborn as sns
from pytraj import matrix
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import os

# this script uses pytraj so download it here: 
'''https://amber-md.github.io/pytraj/latest/installation.html'''

# function to download pdb
def download_pdb(pdb_list):
  '''
  function takes list of pdb code, you can download from pdb website
  pdb_list is like 1HHP, 3TKG, 3TKW, 3TL9
  '''
  data=[]
  with open(pdb_list) as prlist:
    pr = csv.reader(prlist, delimiter=' ')
    for row in pr:
      data.append(row)
  num_of_pdb = len(data[0])
  # make a directory PDB
  os.system("mkdir PDB")
  #get each pdb 
  for i in np.arange(0, num_of_pdb, 1):
    t = data[0][i].strip(',') #get the fourlettercode
    lower = t.lower() #change to lower case
    pdb = pt.io.loadpdb_rcsb("%s" %lower)
    #write pdb
    pt.write_traj("PDB/%s.pdb" %lower, pdb, overwrite=True)

# function gets the mask for the dihedral
def dihmask(res, dih, name):
    '''
    res is the residue number e.g 8
    dih is the dihedral e.g 'phi' 'psi' 'chi1' 'chi2' 'chi3' etc
    name is the three letter code of amino acid
    '''
    # do phi and psi
    if dih == 'phi': 
        if (res -  1) != 0: # phi cannot be on 1st residue
            beres = res - 1
            return ':%s@C :%s@N :%s@CA :%s@C' %(beres,res,res,res)
        else: 
            return 'empty' 
    if dih == 'psi':
        afres = res + 1
        return ':%s@N :%s@CA :%s@C :%s@N' %(res,res,res,afres)
    #list for chi1
    list1 = ['ARG', 'ASN', 'ASP', 'GLN', 'GLU', 'HIS', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'TRP', 'TYR']
    #list for chi2 and chi3 
    list2 = ['ARG', 'GLN', 'GLU', 'LYS', 'PRO']
    list3 = ['TYR', 'TRP', 'LEU', 'PHE']
    
    # do the chi1
    if dih == 'chi1': 
        if name in list1:
            return ':%s@N :%s@CA :%s@CB :%s@CG' %(res,res,res,res)
        elif (name == 'ILE') or (name == 'VAL'):
            return ':%s@N :%s@CA :%s@CB :%s@CG1' %(res,res,res,res)
        elif name == 'CYS':
            return ':%s@N :%s@CA :%s@CB :%s@SG' %(res,res,res,res)
        elif name == 'SER':
            return ':%s@N :%s@CA :%s@CB :%s@OG' %(res,res,res,res)
        elif name == 'THR':
            return ':%s@N :%s@CA :%s@CB :%s@OG1' %(res,res,res,res)
        else: 
            return 'empty' 
    # do the chi2
    elif dih == 'chi2':
        if name in list2: 
            return ':%s@CA :%s@CB :%s@CG :%s@CD' %(res,res,res,res)    
        elif name in list3:
            return ':%s@CA :%s@CB :%s@CG :%s@CD1' %(res,res,res,res)
        elif (name == 'ASN') or (name == 'ASP'):
            return ':%s@CA :%s@CB :%s@CG :%s@OD1' %(res,res,res,res)
        elif name == 'HIS':
            return ':%s@CA :%s@CB :%s@CG :%s@ND1' %(res,res,res,res)
        elif name == 'ILE':
            return ':%s@CA :%s@CB :%s@CG1 :%s@CD' %(res,res,res,res)        
        elif name == 'MET':
            return ':%s@CA :%s@CB :%s@CG :%s@SD' %(res,res,res,res)
        else: 
            return 'empty' 
    # do the chi3
    elif dih == 'chi3':
        if name == 'MET': 
            return ':%s@CB :%s@CG :%s@SD :%s@CE' %(res,res,res,res)
        elif (name == 'GLN') or (name == 'GLU'):
            return ':%s@CB :%s@CG :%s@CD :%s@OE1' %(res,res,res,res)  
        elif name == 'ARG':
            return ':%s@CB :%s@CG :%s@CD :%s@NE' %(res,res,res,res)
        elif name == 'LYS':
            return ':%s@CB :%s@CG :%s@CD :%s@CE' %(res,res,res,res)
        elif name == 'PRO':
            return ':%s@CB :%s@CG :%s@CD :%s@N' %(res,res,res,res)
        elif name == 'TYR':
            return ':%s@CB :%s@CG :%s@CD1 :%s@CE1' %(res,res,res,res)
        else: 
            return 'empty' 
    # do the chi4
    elif dih == 'chi4':
        if name == 'ARG':
            return ':%s@CG :%s@CD :%s@NE :%s@CZ' %(res,res,res,res)
        elif name == 'LYS':
            return ':%s@CG :%s@CD :%s@CE :%s@NZ' %(res,res,res,res)
        elif name == 'TYR':
            return ':%s@CG :%s@CD1 :%s@CE1 :%s@CZ' %(res,res,res,res)
        else: 
            return 'empty' 
    elif dih == 'chi5':
        if name == 'ARG':
            return ':%s@CD :%s@NE :%s@CZ :%s@NH1' %(res,res,res,res)
        elif name == 'TYR':
            return ':%s@CD1 :%s@CE1 :%s@CZ :%s@OG' %(res,res,res,res)
        else: 
            return 'empty' 

# function calculate a dihedral for each frame
def cal_dih(traj, mask):
  '''mask is dihedral mask'''
  data = [] 
  if mask != 'empty': # not all residues have same chis, so return empty for those using the dihmask function above
    data = pt.dihedral(traj, mask)
  return data    

# function to create data frame of all dihedral value for a pdb in a list of pdb 
def create_PDBdf(pdb_list, tot_res):
  '''
  pdb_list is a list of pdb e.g ['2qnp', '1hhp']. They must be stored in a directory name PDB
  tot_res is the total number of residues
  function return a dictionary name pdb. it can be called as pdb[pdbname]
  e.g pdb['1hhp']
  '''
  # populate a dataframe for each PDB that contains all the dihedrals 
  pdb = {} #dictionary 
  #get each pdb
  # strip the less using comma delimiter
  list_pdb = np.genfromtxt(pdb_list, dtype=str, delimiter=',')
  for i in list_pdb:
    # convert to lower case
    i = i.lower()
    i = i.strip()
    # get the PDB path
    #print (i) 
    j = 'PDB/%s.pdb' %(i)
    # define PDB as topology using pytraj
    top = pt.load_topology(j)
    traj = pt.iterload(j)
    # set up a dataframe that have residue#, residue name, and the dihedrals up to chi5
    data_df = pd.DataFrame(index = np.arange(1, tot_res, 1), columns=["res#", "resname", "phi", "psi", "chi1", "chi2", "chi3", "chi4", "chi5"])
    # populate the column residue #
    data_df['res#'] = np.arange(1, tot_res, 1)
    # we will calculate the dihedrals for each residued and make a map, k is residue #
    for k in np.arange(1, tot_res, 1): 
      # cheap method to get the residues, a hack but works
      residues, ss, _ = pt.dssp(traj, ":%s" %(k))
      res = residues[0][0:3]
      #print (res)
      # put the residue name in the dataframe
      data_df.at[k, 'resname'] = '%s' %(res) 
      # do phi first
      if k != 1 : #no phi for first residue 
        # get the index for phi 
        indx = '%s' %(dihmask(int(k), 'phi', res)) 
        # calculate the dihedral value, and store in data
        data = cal_dih(traj, indx)
        data_df.at[k, 'phi'] = data[0]
        #print (data)
      # do psi after
      if k != tot_res: # no psi for last residue
        # get the index for pytraj dihedral function
        indx = '%s' %(dihmask(int(k), 'psi', res)) 
        # calculate the dihedral value
        data = pt.dihedral(traj, indx)
        data_df.at[k, 'psi'] = data[0]
        #print (data)
      # now do chi's
      chis = ['chi1', 'chi2', 'chi3', 'chi4', 'chi5']
      for chi in chis: 
        indx = '%s' %(dihmask(int(k), chi, res)) 
        # calculate the dihedral value
        if indx != 'empty': 
          data = pt.dihedral(traj, indx)
          data_df.at[k, chi] = data[0]
    #make a directory to store reference value
    os.system("mkdir reference_values")
    data_df.to_csv("./reference_values/%s_ref.dat" %(i), float_format='%.4f') 
    pdb[i] = data_df # set pdb dataframes
  return pdb

# function 3 compare two pdb and print the difference in dihedrals
# user can use the mask for these dihedrals and check if these differences are seen in a trajectory
# this is good if you have a pdb of structure at the begining and end of a reaction pathways 
# such as in closed and semi-open in HIV-1 Pr 
def compare_2PDB(pdb1, pdb2, tot_res, filename):
  '''
  pdb1 and pdb2 are two pdb structures in the PDB directory
  pdb is a dictionary that hold all the PDBs dihedral. It is the output of running compare_PDB 
  returns diff_df an array 
  '''
  # set the diff dataframe to be same as anyone of them, will change values
  
  diff_df = pdb1
  dih = ['phi', 'psi', 'chi1', 'chi2', 'chi3', 'chi4', 'chi5']
  for name in dih: 
    # name is the dihedral 
    diff=[]
    for val, val1 in zip(pdb1[name], pdb2[name]):
      if np.sign(val) == np.sign(val1): # so they have the same sign
        diff1 = abs(val) - abs(val1)
        diff1 = abs(diff1)
        diff.append(diff1)
      else: # sign is different ,to deal with values on other end of the -180 to 180 spectrum
  
        # sum the distance from 180 for each number 
        diff1 = (180 - abs(val)) + (180 - abs(val1))
        # sum the distance from 0 for each number
        diff2 = (abs(val) - 0) + (abs(val1 - 0))
        # take the smalleer difference
        if diff1 <= diff2:
          diff.append(diff1)
        else:
          diff.append(diff2)
    diff_df[name] = diff
    diff = 0.0
  os.system("mkdir difference_values")
  diff_df.to_csv("./difference_values/%s.dat" %(filename), float_format='%.4f') 
  return diff_df

#function to get all dihedral masks in a reference pdb
def get_allmasks(pdb, tot_res):
  # function will break for psi with last reside got to find a way to delete last psi 
  df = pdb
  row_indx = np.arange(0, tot_res, 1) # this is the row in diff_df dataframe
  dih = ['phi', 'psi', 'chi1', 'chi2', 'chi3', 'chi4', 'chi5'] #columns in diff_df
  indx = []; res= []
  # write a file
  # for each dihedral name
  for name in dih:
    # for each cell
    for ind, val in zip(row_indx, df[name]):
      resn =  (df.iloc[ind]['res#']) # find the cell with row ind and column res#
      resnam = (df.iloc[ind]['resname'])
      indx_traj = '%s' %(dihmask(int(resn), name, resnam))
      if indx_traj != 'empty':
        indx.append(indx_traj)
        dihedral = '%s%s_%s' %(resnam,resn,name)
        res.append(dihedral)
  return indx, res

# this function takes the difference and return the mask for the dihedrals that differ by some thereshold
def PDB_diff_mask(diff_df, threshold, tot_res, filename):
  ''' diff_df is a dataframe with the difference between two PDB dihedrals value.
      it is the output of compare_2PDB'''

  row_indx = np.arange(0, tot_res, 1) # this is the row in diff_df dataframe
  dih = ['phi', 'psi', 'chi1', 'chi2', 'chi3', 'chi4', 'chi5'] #columns in diff_df
  indx_diff = []; res_diff = []
  # write a file
  f = open("./difference_values/%s.dat" %(filename), "w")
  f.write("# The following dihedrals differ between two structures by %s degrees\n" %(threshold))
  f.write("#Dihedral      Dihedral value      Dihedral mask\n")
  # for each dihedral name 
  for name in dih: 
    # for each cell
    for ind, val in zip(row_indx, diff_df[name]):
      if val > threshold: 
        resn =  (diff_df.iloc[ind]['res#']) # find the cell with row ind and column res#
        resnam = (diff_df.iloc[ind]['resname'])
        indx_traj = '%s' %(dihmask(int(resn), name, resnam)) 
        indx_diff.append(indx_traj)
        dihedral = '%s%s_%s' %(resnam,resn,name)
        res_diff.append(dihedral)
        f.write("%s       %.4f         %s\n" %(dihedral, val, indx_traj))
  return indx_diff, res_diff
      
# plot a bunch of dihedral time series of a trajectory based on mask 
def plot_dih(top, crd, ref_pdb, mask_array, res_array, fig_title, min_array = [0, -180], max_array = [1, 180]):
  '''
  top is topology path
  crd is trajectory path
  ref_pdb is the reference PDB path 
  mask_array is an array with dihedral masks you want to plot
  fig_title is the title for plot 
  min_array is the array that has the minimum value on the x and y axis, allows the ability to zoom in on a plot 
    default to x = 0, y = -180 
  max_array is the array that has the minimum value on the x and y axis, allows the ability to zoom in on a plot 
    default to x = len of trajectory , y = 180 
  '''
  os.system("mkdir images")
  import seaborn as sns; sns.set(style="ticks")
 # create traj and reference traj
  traj = pt.iterload([crd], top=top)    
  ref_traj = pt.load(ref_pdb)   

  # if the maximum limit on the x axis is not specified change it to length of traj array  
  if max_array[0] == 1: # so if it is default
    max_array[0] = len(traj)
  # histogram parameters
  min_val=min_array[1]; max_val = max_array[1]; bin_width=2.0
  # loop through all mask in mask_array and plot dihedral
  for indx, restitle in zip(mask_array, res_array): 
    # calculate dihedral of mask indx for trajectory and reference 
    data1 = cal_dih(traj, indx)
    ref_val = cal_dih(ref_traj, indx)

    # plot data
    fig = plt.figure(figsize=(12,6))
    #plot my data
    Yarray = data1
    x=np.arange(0, len(Yarray), 1)

    #gs = GridSpec(1, 2) #no of rows and no of columns in grid, it is equal at default
    # if you want to adjust it you can use the ratios, ratio of 2n, 1m means n will be
    # 2 times longer than m. n is the first subplot and m the second
    gs = GridSpec(1, 2, width_ratios=[4,1], height_ratios=[1])
    # subplot with scatter information      
    ax1 = plt.subplot(gs[0])
    ax1.scatter(x, Yarray, s=0.5, color="blue")
    ax1.axhline(y=ref_val, color='black', linewidth=2)
    ax1.set_ylabel("%s" %(restitle), fontsize=28)
    ax1.set_xlabel("frame #", fontsize=28)
    ax1.tick_params(axis='both', which='both', labelsize=18, labelcolor='purple')
    ax1.set_xlim(min_array[0], max_array[0])
    ax1.set_ylim(min_array[1], max_array[1])
    plt.grid(b=True, which='both', color='gray', linewidth=0.13, linestyle='--')
    
    # subplot with histogram
    ax2 = plt.subplot(gs[1]) #, sharey=ax1) # neeto share axis with 1
    ax2.hist(Yarray[min_array[0]:max_array[0]], bins=np.arange(min_val, max_val, bin_width), orientation='horizontal')
    ax2.axhline(y=ref_val, color='black', linewidth=2)
    ax2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) # turn off labels and ticks
    ax2.yaxis.set_ticklabels([]) #hide the labels 
    ax2.set_ylim(min_array[1], max_array[1])
  
    plt.grid(b=True, which='both', color='gray', linewidth=0.13, linestyle='--')
    plt.savefig('./images/%s_%s.png' %(restitle,fig_title))

# use this plot if you want to plot two ref pdb value that u know differe 
def plot_dih2ref(top, crd, ref_pdb1, ref_pdb2, mask_array, res_array, fig_title, min_array = [0, -180], max_array = [1, 180]):
  '''
  top is topology path
  crd is trajectory path
  ref_pdb is the reference PDB path 
  mask_array is an array with dihedral masks you want to plot
  fig_title is the title for plot 
  min_array is the array that has the minimum value on the x and y axis, allows the ability to zoom in on a plot 
    default to x = 0, y = -180 
  max_array is the array that has the minimum value on the x and y axis, allows the ability to zoom in on a plot 
    default to x = len of trajectory , y = 180 
  '''
  os.system("mkdir images")
  import seaborn as sns; sns.set(style="ticks")
  # histogram parameters
  min_val=min_array[1]; max_val = max_array[1]; bin_width=2.0
  # create traj and reference traj
  traj = pt.iterload([crd], top=top)    
  ref_traj1 = pt.load(ref_pdb1)    
  ref_traj2 = pt.load(ref_pdb2)    
  # if the maximum limit on the x axis is not specified change it to length of traj array  
  if max_array[0] == 1: # so if it is default
    max_array[0] = len(traj)
   # histogram parameters
  min_val=min_array[1]; max_val = max_array[1]; bin_width=2.0
 
  # loop through all mask in mask_array and plot dihedral
  for indx, restitle in zip(mask_array, res_array): 
    # calculate dihedral of mask indx for trajectory and reference 
    data1 = cal_dih(traj, indx)
    ref_val1 = cal_dih(ref_traj1, indx)
    ref_val2 = cal_dih(ref_traj2, indx)
    # plot data
    fig = plt.figure(figsize=(12,6))
    #plot my data
    Yarray = data1
    x=np.arange(0, len(Yarray), 1)

    #gs = GridSpec(1, 2) #no of rows and no of columns in grid, it is equal at default
    # if you want to adjust it you can use the ratios, ratio of 2n, 1m means n will be
    # 2 times longer than m. n is the first subplot and m the second
    gs = GridSpec(1, 2, width_ratios=[4,1], height_ratios=[1])

    ax1 = plt.subplot(gs[0])
    ax1.scatter(x, Yarray, s=0.5, color="blue") #, label="%s" % (label1))
    ax1.axhline(y=ref_val1, color='black', linewidth=2)
    ax1.axhline(y=ref_val2, color='red', linewidth=2)
    ax1.set_ylabel("%s" %(restitle), fontsize=28)
    ax1.set_xlabel("frame #", fontsize=28)
    ax1.tick_params(axis='both', which='both', labelsize=18, labelcolor='purple')
    ax1.set_xlim(min_array[0], max_array[0])
    ax1.set_ylim(min_array[1], max_array[1])
    plt.grid(b=True, which='both', color='gray', linewidth=0.13, linestyle='--')

    ax2 = plt.subplot(gs[1]) #, sharey=ax1) # neeto share axis with 1
    ax2.hist(Yarray[min_array[0]:max_array[0]], bins=np.arange(min_val, max_val, bin_width), orientation='horizontal')
    ax2.axhline(y=ref_val1, color='black', linewidth=2)
    ax2.axhline(y=ref_val2, color='red', linewidth=2)
    ax2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) # turn off labels and ticks
    ax2.yaxis.set_ticklabels([]) #hide the labels 
    ax2.set_ylim(min_array[1], max_array[1])
  
    plt.grid(b=True, which='both', color='gray', linewidth=0.13, linestyle='--')
    plt.savefig('./images/%s_%s.png' %(restitle,fig_title))
 

# plot dihedrals that differ between two trajectories
def plot_scat_2traj(traj1, traj2, ref_pdb, mask_array, res_array, thresh, plot_thresh, label_traj1, label_traj2, title):
  '''
   thresh allow the count, so if a given bin have normalize values that differ by 0.8 then add 1 to the counter
   recommend using a number close to 1 (0.9), if you looking for a trajectory that might sample a given range 
   of dihedral (-60 to -90) but the  the other trajectory do not sample that dihedral range at all. 
   use 0.5 if you want to know if both trajectory sample the same dihedral range but one trajectory sampled
   more of that dihedral range than the other.
   ref_pdb is pdb path e.g ./PDB/1hhp.pdb
   #only plot data that have an overall_frac of 0.30, if plot_thresh is 0.000001 then all dihedrals will be plot
   plot_thresh is the percentage of bins that are different.
   plot_thresh = 0.30  #(57/144 bins are completely different)
  '''
 
  ref_traj = pt.load(ref_pdb)    
  
  os.system("mkdir images_traj_diff")
  for indx, restitle in zip(mask_array, res_array):
    print (indx)
    # get the dihedral data based on atom mask indx
    data1 = cal_dih(traj1, indx)
    data2 = cal_dih(traj2, indx)
    ref_val = cal_dih(ref_traj, indx)
  
    # create histogram of both dataset 
    min_val=-180; max_val = 180; bin_width=2
    h1 = np.histogram(data1, bins=np.arange(min_val, max_val, bin_width), density=True)
    h2 = np.histogram(data2, bins=np.arange(min_val, max_val, bin_width), density=True)
   
    # get the number of bins in any one the dataset. they are both equal
    num_bins = len(h1[1])
    # the count for the difference between values in each bin
    diff_count = 0  
    # loop through values in each bins, bins have to be similar
    # h1[0] and h2[0] are arrays that have the values in each bin 
    for val1, val2 in zip(h1[0], h2[0]): 
      # take the max of the two counts in each bin
      max_val = max(val1, val2) # e.g val1 = 50 and val2 = 100
      # normalize by dividing both values for a given bin by max val
      val1 = val1 / max_val # e.g 50/100
      val2 = val2 / max_val # e.g 100/100
      # find the difference between normalize values
      val_diff = abs(val1-val2) # e.g |0.5-1|
      # if the value is greater than threshold then count it as a bin that differs 
      if val_diff > thresh: 
        diff_count += 1
      else: 
        diff_count = diff_count
   
    print ("difference_count: " + str(diff_count))
    # normalize the diff count
    overall_frac = diff_count/num_bins
    print ("overall_fraction: " + str(overall_frac))
    # An overall_frac of 1 means all bins are different, 0 means all bins are the same (same based on thresh)

    if overall_frac > plot_thresh:
      #plot my data
      import seaborn as sns; sns.set(style="ticks")
      fig = plt.figure(figsize=(16, 10))
      x1=np.arange(0, len(data1), 1)
      x2=np.arange(0, len(data2), 1)
      gs = GridSpec(2, 4) #no of rows and no of columns in grid
  
      ax1 = plt.subplot(gs[0, :-1])
      ax1.scatter(x1, data1, s=5, color="blue") #, label="%s" % (label1))
      ax1.axhline(y=ref_val, color='black', linewidth=2)
      ax1.set_title("%s" %label_traj1, fontsize=28)
      ax1.set_ylabel("%s" %(restitle), fontsize=28)
      ax1.set_yticks(np.arange(-180, 180, 40))
      ax1.tick_params(axis='both', which='both', labelsize=20, labelcolor='purple')
      ax1.set_ylim(-180, 180)
      ax1.set_xlim(0, max(x1))

   
      ax2 = plt.subplot(gs[1, :-1])
      ax2.scatter(x2, data2, s=5, color="red") #, label="%s" % (label2))
      ax2.axhline(y=ref_val, color='black', linewidth=2)
      ax2.set_title("%s" %label_traj2, fontsize=28)
      ax2.set_ylabel("%s" %(restitle), fontsize=28)
      ax2.set_xlabel("Frame #", fontsize=28)
      ax2.set_yticks(np.arange(-180, 180, 30))
      ax2.tick_params(axis='both', which='both', labelsize=20, labelcolor='purple')
      ax2.set_ylim(-180, 180)
      ax2.set_xlim(0, max(x2))

      ax3 = plt.subplot(gs[1, -1]) #sharey=ax2) # need to share axis with 2
      ax3.hist(data2, bins=np.arange(-180, 180, 1), orientation='horizontal')
      ax3.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
      ax3.yaxis.set_ticklabels([]) #hide the labels 
      ax3.set_yticks(np.arange(-180, 180, 30))
      ax3.set_ylim(-180, 180)

      ax4 = plt.subplot(gs[:-1, -1]) # sharey=ax1) # neeto share axis with 1
      ax4.hist(data1, bins=np.arange(-180, 180, 1), orientation='horizontal')
      ax4.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) # turn off labels and ticks
      ax4.yaxis.set_ticklabels([]) #hide the labels 
      ax4.set_yticks(np.arange(-180, 180, 40))
      ax4.set_ylim(-180, 180)
  
      plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
      plt.savefig('./images_traj_diff/%s_%s.png' %(restitle, title))


#plot the deviations of a trajectory from a reference PDB for a given dihedral dih. return a file with frame number at first deviation 
def compare_traj2PDB(traj, ref_pdb, tot_res, dih, residue_array=[1,0]):
  '''
  traj is your pytraj trajectory to compare to the pdb structure 
   generated in the python script as traj = pt.iterload([crd], top=top)    
  tot_res is total residue 
  ref_pdb is pdb path e.g ./PDB/1hhp.pdb
  residue_array is an array has the minimum and maximum residue number
  if residue_array is not given then function use from residue 1 to total residue
  '''
  #traj = traj[0:99] # 1 ns
  r0 = residue_array[0] # the lowest residue #
  # the last or highest residue is total residue if the residue_array is not specify
  if residue_array[1] == 0:
    r1 = tot_res
  else:
    r1 = residue_array[1] # the highest residue #

  # reference traj which is the pdb, if you are using a trajectory frame, use cpptraj and make that 
  # frame into a pdb  
  ref_traj = pt.load(ref_pdb)    
  # length of trajectory
  traj_len = len(traj)
  # set up a dataframe that frames on the y axis (index) and dihedral, dih, of a residue# on x axis (columns)
  data_df = pd.DataFrame(index = np.arange(1, traj_len+1, 1)) 
  # dataframe that holds the reference dihedral value 
  data_df_ref = pd.DataFrame()
  # for each residue in ....
  for k in np.arange(r0, r1+1, 1): 
    # for phi
    if (dih == 'phi' and k != 1): # no phi for first residue 
      # cheap method to get the residues, a hack but works
      # res is the residue name and k is the residue number 
      residues, ss, _ = pt.dssp(ref_traj, ":%s" %(k))
      res = residues[0][0:3]
      # get the index for phi 
      indx = '%s' %(dihmask(int(k), dih, res)) 
      # calculate the phi dihedral value for residue k in the trajectory and store in data1
      data1 = cal_dih(traj, indx)
      # create a column name for the residue. e.g gly17_phi and store values in column
      data_df['%s%s_%s'%(res,k,dih)] = data1
      # repeat for reference, should remove out of the loop?
      data_ref1 = cal_dih(ref_traj, indx)
      data_df_ref['%s%s_%s'%(res,k,dih)] = data_ref1
    # set to nan for the 1st residue if it is a phi dihedral
    elif (dih == 'phi' and k ==1): 
      # res is the residue name and k is the residue number 
      residues, ss, _ = pt.dssp(ref_traj, ":%s" %(k))
      res = residues[0][0:3]
      data_df['%s%s_%s'%(res,k,dih)] = "Nan"
      data_df_ref['%s%s_%s'%(res,k,dih)] = "Nan"
    # for psi 
    elif (dih == 'psi' and k != tot_res): # no psi for last residue
      # cheap method to get the residues, a hack but works
      # res is the residue name and k is the residue number 
      residues, ss, _ = pt.dssp(ref_traj, ":%s" %(k))
      res = residues[0][0:3]
      # get the index for psi 
      indx = '%s' %(dihmask(int(k), dih, res)) 
      # calculate the psi dihedral value for residue k in the trajectory and store in data1
      data1 = cal_dih(traj, indx)
      # create a column name for the residue. e.g gly17_psi and store values in column
      data_df['%s%s_%s'%(res,k,dih)] = data1
      # repeat for reference, should remove out of the loop?
      data_ref1 = cal_dih(ref_traj, indx)
      data_df_ref['%s%s_%s'%(res,k,dih)] = data_ref1
    # set to nan for the last residue if it is a psi dihedral
    elif (dih == 'phi' and k == tot_res): 
      # res is the residue name and k is the residue number 
      residues, ss, _ = pt.dssp(ref_traj, ":%s" %(k))
      res = residues[0][0:3]
      data_df['%s%s_%s'%(res,k,dih)] = "Nan"
      data_df_ref['%s%s_%s'%(res,k,dih)] = "Nan"
    # dihedral is a chi dihedral 
    else: 
      # cheap method to get the residues, a hack but works
      # res is the residue name and k is the residue number 
      residues, ss, _ = pt.dssp(ref_traj, ":%s" %(k))
      res = residues[0][0:3]
      # get the index for psi 
      indx = '%s' %(dihmask(int(k), dih, res)) 
      # calculate the chi dihedral value for residue k in the trajectory and store in data1
      data1 = cal_dih(traj, indx)
      if ( (len(data1)) != 0): # if there is values in data1 then there is a chi
        # create a column name for the residue. e.g pro79_chi1 and store values in column
        data_df['%s%s_%s'%(res,k,dih)] = data1
        # repeat for reference, should remove out of the loop?
        data_ref1 = cal_dih(ref_traj, indx)
        data_df_ref['%s%s_%s'%(res,k,dih)] = data_ref1
      else: # enter Nan for the missing chi 
        data_df['%s%s_%s'%(res,k,dih)] = "Nan"
        data_df_ref['%s%s_%s'%(res,k,dih)] = "Nan"
  # create dih data frame by appending data_df with data_df_ref
  dih_val = data_df.append(data_df_ref)
  # name the index frame
  dih_val.index.name = "Frame"
  # rename the reference row ref 
  dih_val.rename(index={0:'ref'},inplace=True)
  # write out dihedral value to 4 decimal place
  os.system("touch %s_values.dat" %(dih))
  dih_val.to_csv("%s_values.dat" %(dih), float_format='%.4f') 
  # now let us get the deviation 
  dev_df = pd.DataFrame(index = np.arange(1, traj_len+1, 1)) # deviation dataframe
  # for values in given residue calculate the deviation to a ref value 
  for k in np.arange(r0, r1+1, 1): 
    diff=[]
    residues, ss, _ = pt.dssp(ref_traj, ":%s" %(k))
    res = residues[0][0:3]
    val1 = data_df_ref['%s%s_%s'%(res,k,dih)].values # reference for dihedral of residue k 
    # for value at each frame, calculate the deviation to the reference value val1
    for val in data_df['%s%s_%s'%(res,k,dih)].values:
      if (val != 'Nan'): #check one, either reference or traj val should be always Nan

      else: # sign is different ,to deal with values on other end of the -180 to 180 spectrum
        if np.sign(val) == np.sign(val1): # so they have the same sign
          diff1 = abs(val) - abs(val1)
          diff1 = abs(diff1)
          diff.append(diff1[0])
        else: # sign is different ,to deal with values on other end of the -180 to 180 spectrum
          # sum the distance from 180 for each number 
          diff1 = (180 - abs(val)) + (180 - abs(val1))
          # sum the distance from 0 for each number
          diff2 = (abs(val) - 0) + (abs(val1 - 0))
          if diff1 <= diff2:
            diff.append(diff1[0])
          else:
            diff.append(diff2[0])
      else: 
        diff.append('nan')
    dev_df['%s%s_%s'%(res,k,dih)] = diff
    diff = 0.0 # zero diff since we are appending 
  os.system("touch %s_values_dev.dat" %(dih))
  dev_df.to_csv("%s_values_dev.dat" %(dih), float_format='%.4f') 

  return dev_df, dih_val 
   
# function to plot deviations  as heat map
def plot_dev_heatmap(dev_df, vmin, vmax, fig_title,  x_lim):
  '''
  dev_df is the deviation dataframe which is one of the output from compare_traj2PDB
  vmin is the lowest limit in the deviation you want  
  vmax is the highest limit in the deviation you want
  x_lim is an array that set the limits on the frames you want to focus on it take ymin, ymax
  # if you want to focus on the residues (y axis) set it in compare_traj2PDB
  '''
  #http://seaborn.pydata.org/generated/seaborn.heatmap.html
  import seaborn as sns
  os.system("mkdir images")
  xmin = x_lim[0]; xmax = x_lim[1]
  # deal with Nan values convert them to 0
  for col in dev_df:
    dev_df[col] = pd.to_numeric(dev_df[col], errors='coerce')

  # transpose x and y 
  dev_df_transpose = dev_df.transpose()
  # plot figures
  fig = plt.figure(figsize=(8, 8))
  heatmap = sns.heatmap(dev_df_transpose, vmin=vmin, vmax=vmax, cmap='Reds', cbar=True)
  heatmap.set_xlim(xmin, xmax)
  plt.savefig('./images/%s.png' %(fig_title))

