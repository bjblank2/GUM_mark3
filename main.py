__author__ = 'brian'
import compile_vasp_structures as cvs
import postprocess_vasp as ppv
import clustermag_rules as cmr  
import calc_fitting_params as cfp
#import mc_functions as mc
import mc_functions_2 as mc2
import mc_supercell as ms
import sys
#
############################
## Some notes from Elif's Modifications to Parameter Fitting
##
## *) As implemented right now,
##      several parts of this code assume that the VASP POSCAR/CONTCAR always enter Ni,Mn,In in the
##      same order.  You will be in trouble on many levels if this is not actually the case.
##      For instance, atom.py assumes the species order is NiMnIn always.
##
## *) Generally, would be good to slightly overhaul m_structure so that it ONLY contains the
##      information necessary to do the fitting and to run MC simulations.  It is carrying around
##      a lot of information we don't need.  Cleaner to leave it behind when its not needed.
##
## *) Look at m_structure.py line 86 for minor question.
##
## *) When doing the fitting, is selecting an intercept a problem? cluster expansion model inherently
##      expects non-centered data?  Otherwise it spends it energy pulling out the linear dependence on composition.
##      Note: fitting is currently modified so intercept is not fit, need to change cfp.ridge_simple if we want
##      to include it.  But I think we do not.
##################################################

aust_tol = 0.02
spin_style = ['threshold','threshold','threshold']  # options for spin_tol. Assuming [Ni Mn In]. choose 'threshold' or 'factor'
spin_tol = [0.05,2,0]                                # insert spin parameters here, this assumes [Ni Mn In ]
species = ['Ni','Mn','In']                          # this is the order that the post-processed data is reported, NEEDS TO BE Heusler format Ni2MnIn, Ni2FeGa.

root_dir = './NiMnIn_Vasp_Data'                     # where the VASP directories are  #root_dir = '/Users/brian/Desktop/folder'
vasp_data_file = './NiMnIn_Data'               # generated in compile_vasp_structures>import_vasp, summarizes output of all VASP calculations
vasp_data_file_pp = './NiMnIn_Data_pp'              # post-processed version of VASP results with spins, positions selected
cluster_file = './Cluster_Rules'                    # cluster expansion rules
j_file = './J_Rules'                                # heisenberg rules
#fitting_structures_file = './'

# Determine what needs to be generated from scratch
vasp_data_exists = True                            # should be False if I want to regenrate files
vasp_pp_exists = False                              # postprocessing of VASP results or no?
Cluster_rules_exist = True                          # define cluster rules
J_rules_exist = True                                # define heisenberg rules                                  # results of fitting model
Fitting_params_exist = False

# write Cluster and J rules file if doesn't exist, then read the rules
if Cluster_rules_exist is False:                    # writes cluster rules if doesn't already exist
    cmr.write_cluster_rules(cluster_file)
if J_rules_exist is False:                          # writes j_rules if doesn't already exist
    cmr.write_j_rules(j_file)
Cluster_rules = cmr.read_cluster_rules(cluster_file)
J_rules = cmr.read_j_rules(j_file)
if Fitting_params_exist == False:
    # summarize VASP data
    if vasp_data_exists is False:                    # will make summary of VASP results if it doesn't already exist
        cvs.import_vasp(root_dir, vasp_data_file,species)

    # Read the VASP datafile and initialize a structure object for each
    # postprocess according to user selected parameters above and the cluster and j rules
    # calculation of sums and checking for duplicates occurs in here now
    # if a given structure is considered a duplicate then it is not added to the structure_list
    M_structures = ppv.generate_m_structure(vasp_data_file, len(Cluster_rules), len(J_rules), aust_tol, spin_style, spin_tol, Cluster_rules, J_rules)

    ppv.write_structures_processedvasp(M_structures,vasp_data_file_pp)
    ppv.summarize_classification(M_structures)
    warning_threshold = 0.5
    ppv.summarize_fitting_structures(M_structures,warning_threshold)

    ## Ridge Regression Fitting with Regularization
    ppv.scale_enrg(M_structures)
    Js,intercept = cfp.ridge_simple(M_structures,1,Cluster_rules,J_rules)
    cfp.write_fitting_parameters(M_structures, Cluster_rules, J_rules, Js, intercept, 200)
    cfp.plot_data3(M_structures,Cluster_rules,J_rules,Js,intercept,200)

else:
    paramiters_file = open('FittingParameters','r')
    lines = paramiters_file.readlines()
    reading_params = True
    line_index = 1
    Js = []
    while reading_params == True:
        line = lines[line_index]
        if 'Actual Energy' in line or line.split() == []:
            reading_params = False
        else:
            line = line.split()
            Js.append(float(line[2]))
            line_index += 1


##--------------------------------------------------------------#
## Code for running the actual MonteCarlo simulation.
## First create and initialize an array of atom objects for the simulation
## then run MonteCarlo.
## temp_data = open('Temp_data','a')
## temp_data.write('size,Temp,passes,H_avg,mimj_avg,mag_avg,absmag_avg,phase_avg,absphase_avg\n')
## temp_data.close()

x_pts = 4 #|
y_pts = 4 #|- Dimensions of the mc_supercellObj (simulation supercell)
z_pts = 8 #|
phase_init = 'mart' # initial phase configuration
spin_init = 'rand' # initial spin configuration
species_init = 'rand'
num_passes = 20 # number of cluster/wolf moves done
num_sub_passes = 30 # number of spin/species flips done per cluster/wolf move
Temp0 = 150 # initial temperature in K
TempF = 500 # final temperature in K
Temp_inc = 5 # temperature increase per pass in K

rule1 = J_rules[1]
rule2 = Cluster_rules[1]
print(rule1.tag)
print(rule2.tag)

## Initialize an array of atoms with ms.mc_supercellObj(size,species,composition)
## size is (x,y,z)dimensions, species is types of atoms allowed (0=Ni,1=Mn,2=In)
## composition is number of each atom (#Ni,#Mn,#In)
#lattice = ms.mc_supercellObj((x_pts,y_pts,z_pts),(0,1,2),[64,48,16],phase_init,spin_init,species_init)
#sys.setrecursionlimit(lattice.num_sites+2)
## To actually run the simulation use
## mc.run_montecarlo(reference_to_atom_array,number_of_passes,starting_temp, BEG_rules,Cluster_rules,J_rules,plot_figs=TRUE)
## BEG_rules,Cluster_rules,J_rules are objects that determine when and how the fitted parameters are applied
print("Beginning MonteCarlo\n")
#mc2.run_WA_MCA_nested(lattice,num_passes,num_sub_passes,Temp0,Temp_inc,TempF,Cluster_rules,J_rules,Js,do_figs=True)
