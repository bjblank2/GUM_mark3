__author__ = 'brian'
import calc_params as Cp
import mc_functions as mc
import mc_functions_2 as mc2
import mc_supercell as ms
import time
import sys

# This is the main file that calculates the fit and runs the
# MonteCarlo code
#--------------------------------------------------------------#
# Code for finding the fitting parameters.
# First read in the DFT data, generate the appropriate data
# structures and finally write the final parameters to the
# the "Js" list

root_dir = '/Volumes/TOURO/Ni-Fe-Ga/Data_Pts'
data_file = './NiMnIn_Data'
beg_file = './BEG_rules'
cluster_file = './Cluster_Rules'
j_file = './J_Rules'
num_comps = 3
num_species = 3
#alphas = [.000001,.00001,.001,.01,.05,.075,.1,.25,.5,1,5,10,100,1000]
alphas = [.0001,.0005,.0075,.001,.00125,.0025,.00375,.005,.01,.075,.1,.5,1,5,10,100]
# Determine if rules and structure data is defined
Data_file_exists = True
BEG_file_exists = True
Cluster_rules_exist = True
J_ruels_exist = True
Js_exist = False
# Set all rules and structure data
if Data_file_exists is False:
    Cp.import_data(num_species, root_dir, data_file)
if BEG_file_exists is False:
    Cp.write_beg_rules(beg_file)
if Cluster_rules_exist is False:
    Cp.write_cluster_rules(cluster_file)
if J_ruels_exist is False:
    Cp.write_j_rules(j_file)

# Read the data and calculate the fit
BEG_rules = Cp.read_beg_rules(beg_file)
Cluster_rules = Cp.read_cluster_rules(cluster_file)
J_rules = Cp.read_j_rules(j_file)
if Js_exist == True:
    Js = Cp.read_Js(len(J_rules)+len(Cluster_rules)+len(BEG_rules))
else:
    M_structures = Cp.read_m_structure_data(data_file, num_species, len(BEG_rules), len(Cluster_rules), len(J_rules))
    Cp.calculate_sums(M_structures, BEG_rules, Cluster_rules, J_rules)
    Js = Cp.ridge_simple(M_structures,1)
    Cp.write_data(M_structures, 200, Js)
    Cp.write_output(M_structures, BEG_rules, Cluster_rules, J_rules, Js, 200)
    Cp.plot_data3(M_structures,BEG_rules,Cluster_rules,J_rules,Js,200)
print('#######################')
#print(Cp.CV_score(Js,M_structures))
#Cp.CV_score(M_structures,BEG_rules,Cluster_rules,J_rules,alphas)
#print(Cp.CV_score2(M_structures))
print('#######################\n')
#--------------------------------------------------------------#

#--------------------------------------------------------------#
# Code for running the actual MonteCarlo simulation.
# First create and initialize an array of atom objects for the simulation
# then run MonteCarlo.
temp_data = open('Temp_data','a')
temp_data.write('size,Temp,passes,H_avg,mimj_avg,mag_avg,absmag_avg,phase_avg,absphase_avg\n')
temp_data.close()


x_pts = 2 #|
y_pts = 2 #|- Dimensions of the mc_supercellObj (simulation supercell)
z_pts = 4 #|
comp = [8,6,2]
phase_init = 'rand' # initial phase configuration (aust, mart, rand)
spin_init = 'rand' # initial spin configuration
species_init = 'rand'
species_flips = True
num_passes = 30 # number of cluster/wolf moves done
num_sub_passes = 30 # number of spin/species flips done per cluster/wolf move
Temp0 = 280 # initial temperature in K
TempF = 450 # final temperature in K
Temp_inc = 5 # temperature increase per pass in K
Mag_field = 0

# Initialize an array of atoms with ms.mc_supercellObj(size,species,composition)
# size is (x,y,z)dimensions, species is types of atoms allowed (0=Ni,1=Mn,2=In)
# composition is number of each atom (#Ni,#Mn,#In)

#sys.setrecursionlimit(2000)
#lattice = ms.mc_supercellObj((x_pts,y_pts,z_pts),(0,1,2),comp,phase_init,spin_init,species_init)#(64,48,16))
#sys.setrecursionlimit(lattice.num_sites+2)
# To actually run the simulation use
# mc.run_montecarlo(reference_to_atom_array,number_of_passes,starting_temp, BEG_rules,Cluster_rules,J_rules,plot_figs=TRUE)
# BEG_rules,Cluster_rules,J_rules are objects that determine when and how the fitted parameters are applied

#print("Beginning MonteCarlo\n")
##mc2.run_WA_MCA(lattice,num_passes,num_sub_passes,Temp0,Temp_inc,TempF,Mag_field,Cluster_rules,J_rules,Js,species_flips)
#t1 = time.time()
#site_eval = mc2.eval_site_new([0,0,0],lattice,Cluster_rules,J_rules,Js,10,0)
#print(site_eval)
#t2 = time.time()
#print(t2-t1)