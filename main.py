__author__ = 'brian'
import calc_params as Cp
import mc_functions as mc
import mc_supercell as ms
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
# Determine if rules and structure data is defined
Data_file_exists = True
BEG_file_exists = True
Cluster_rules_exist = True
J_ruels_exist = True
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
M_structures = Cp.read_m_structure_data(data_file, num_species, len(BEG_rules), len(Cluster_rules), len(J_rules))
Cp.calculate_sums(M_structures, BEG_rules, Cluster_rules, J_rules)
Js = Cp.ridge_simple(M_structures,1)
Cp.write_data(M_structures, 200, Js)
Cp.write_output(M_structures, BEG_rules, Cluster_rules, J_rules, Js, 200)
Cp.plot_data3(M_structures,BEG_rules,Cluster_rules,J_rules,Js,200)
print('#######################\n')
#print(Cp.CV_score(Js,M_structures))
#print(Cp.CV_score2(M_structures))
print('#######################\n')
#--------------------------------------------------------------#

#--------------------------------------------------------------#
# Code for running the actual MonteCarlo simulation.
# First create and initialize an array of atom objects for the simulation
# then run MonteCarlo.
temp_data = open('Temp_data','w')
temp_data.write('Temp  H_avg  mag_avg  mag2_avg  phase_avg  phase2_avg\n')
temp_data.close()
x_pts = 10
y_pts = 10
z_pts = 20
# Initialize an array of atoms with ms.mc_supercellObj(size,species,composition)
# size is (x,y,z)dimensions, species is types of atoms allowed (0=Ni,1=Mn,2=In)
# composition is number of each atom (#Ni,#Mn,#In)
lattice = ms.mc_supercellObj((x_pts,y_pts,z_pts),(0,1,2),[64,64,0])#(64,48,16))
sys.setrecursionlimit(lattice.num_sites+1)
# To actually run the simulation use
# mc.run_montecarlo(reference_to_atom_array,number_of_passes,starting_temp, BEG_rules,Cluster_rules,J_rules,plot_figs=TRUE)
# BEG_rules,Cluster_rules,J_rules are objects that determine when and how the fitted parameters are applied
print("Beginning MonteCarlo\n")
mc.run_WA_MCA(lattice,100,50 ,100,Cluster_rules,J_rules,Js,do_figs=True)