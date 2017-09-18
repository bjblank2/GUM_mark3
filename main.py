__author__ = 'brian'
import compile_fitting_structures as cfs
import clustermag_rules as cmr  
import calc_fitting_params as cfp
#import mc_functions as mc
#import mc_functions_2 as mc2
#import mc_supercell as ms
import sys

############################
## Some notes from Elif's Modifications to Parameter Fitting
##
## *) We need to modify CFS.IMPORT_DATA which writes NiMnIn_Data to avoid
##      assigning phase and magnetism (avoid human intervention at this stage).
##      This function should only summarize the results of DFT simulations compactly.
##      For now I modified the NiMnIn_Data file the way I want it to appear
##      but have not yet modified the function because I do not have Brian's directory.
##
## *) I eliminated BEG rules entirely.
##
## *) We can now set the tolerance whether a given structure is austenite or not
##      here using variable 'aust_tol'. I think it was hard-coded in, and set too loose before (0.1).
##      Right now I'm going with 0.025 but this also may not be optimal.  I wonder if there is a way
##      to assess this "intelligently" by comparing how well our fitting goes for different tolerances.
##
## *) Adjustment of spins.  We do some assessment when selecting how to define the spins. It is now
##      also no longer hard coded in. This also would be nice to evaluate in an intelligent manner. I
##      did this in a sloppy way though, and generally this could be cleaned up.
##
## *) M_structures now stores the array of species like [Ni,Mn,In]. As implemented right now,
##      several parts of this code assume that the VASP POSCAR/CONTCAR always enter Ni,Mn,In in the
##      same order.  You will be in trouble on many levels if this is not actually the case.
##      For instance, atom.py assumes the species order is NiMnIn always.
##
## *) Generally, would be good to slightly overhaul m_structure so that it ONLY contains the
##      information necessary to do the fitting and to run MC simulations.  It is carrying around
##      a lot of information we don't need.  Cleaner to leave it behind when its not needed.
##
## *) When M_structures are set, right now it just calls everything 'sd' = spin-disordered.
##      Would be good to have the code guess AFM, FM, etc.
##
## *) I added a function in cfs to check for duplicate data sets in NiMnIn_Data based on cluster summing.
##      Note that Brian's data set did have some duplicates, which I have removed from NiMnIn_Data.
##      For now you get warned and you have to manually remove them. In the future, might be nice to
##      just have it done internally in the code. Would be better to leave NiMnIn_Data represent EXACTLY
##      what came from VASP directories, and just have duplicates not included in the structures.
##      Some manual work is required though because need to check that the structures are
##      actually duplicates -- i.e. energies, latt consts are similar.
##
## *) Now, after creating M_structures, a summary of all fitting structures and how they have been assigned
##      is output to file 'summary_fitting_structures'.  Here is where we have post-processed the VASP results
##      to assign phase and spin, and summed the cluster/spin rules for each structure, and removed duplicates.
##
## *) we should make an option to just read cluster and j sums from summary rather than regenerating each time
##
## *) Look at m_structure.py line 86 for minor question.
##
## *) Now we create plots of fitted parameters vs. regularization, and the score vs. the regularization
############################

## *) manually play with spin description. spin selection rules need eyeballing - see Ni, seems to step up in 0.15 units , in atom.py
##      Seems that the way I did it made it worse?
##
## *) SHOULD OUTPUT HERE THE POST-PROCESSED VASP WITHOUT PRIOR TO DOING ANY SUM CALCULATIONS WHICH DEPEND ON THE RULES
##      want to break up cfs into post-process vasp data and then generate fitting structures
##
## *) I have not yet made sure the cluster and j rules are being calculated properly.  I also think we need
##      to make it easy to change them to get the best possible description.
##
## *) need to assess degree of overfitting in model using cross-validation in sklearn. Can we improve the model?
##
## *) When doing the fitting, is selecting an intercept a problem? cluster expansion model inherently
##      expects non-centered data?  Otherwise it spends it energy pulling out the linear dependence on composition.
##      Note: fitting is currently modified so intercept is not fit, need to change cfp.ridge_simple if we want
##      to include it.  But I think we do not.

aust_tol = 0.025
spin_style = ['threshold','threshold','threshold']  # options for spin_tol. Assuming [Ni Mn In]. choose 'threshold' or 'factor'
spin_tol = [0.3,2.5,0]                        # insert spin parameters here, this assumes [Ni Mn In ]

root_dir = '/Volumes/TOURO/Ni-Fe-Ga/Data_Pts'  # where the VASP directories are
data_file = './NiMnIn_Data_NoDups'             # generated in calc_params>import_data, summarizes output of all VASP calculations
data_file_pp = './NiMnIn_Data_NoDups_pp'       # post-processed version of VASP results with spins, positions selected
cluster_file = './Cluster_Rules'               # cluster expansion rules
j_file = './J_Rules'                           # heisenberg rules
fitting_structures_file = './'

# Determine what needs to be generated from scratch
vasp_summary_exists = True                     # summarize VASP results from VASP directories or no?
Cluster_rules_exist = True                     # define cluster rules
J_rules_exist = True                           # define heisenberg rules
Js_exist = True                                # results of fitting model

# set all rules and summarize VASP data
if vasp_summary_exists is False:               # will make summary of VASP results if it doesn't already exist
    cfs.import_vasp(root_dir, data_file)
if Cluster_rules_exist is False:               # writes cluster rules if doesn't already exist
    cmr.write_cluster_rules(cluster_file)
if J_rules_exist is False:                     # writes j_rules if doesn't already exist
    cmr.write_j_rules(j_file)

# Read the summarized VASP data and rules files, and initialize a structure object for each
# data set without doing sum rules yet.
Cluster_rules = cmr.read_cluster_rules(cluster_file)
J_rules = cmr.read_j_rules(j_file)
M_structures = cfs.generate_m_structure(data_file, len(Cluster_rules), len(J_rules), aust_tol, spin_style, spin_tol)
cfs.write_structures_processedvasp(M_structures,data_file_pp)

# Evaluate cluster and spin sums here, and check for duplicates
# Seems like there should be the option to read the sums from the
# summary_fitting_structures file to avoid doing this summing each time.
cfs.calculate_sums(M_structures, Cluster_rules, J_rules, spin_style, spin_tol)
if (cfs.check_duplicate_structures(M_structures)=='True'):
    print ('Based on summed cluster and spin rules, fitting structures appear to contain duplicates. \n')
cfs.summarize_fitting_structures(M_structures)

## Ridge Regression Fitting with Regularization
Js,intercept = cfp.ridge_simple(M_structures,1)
cfp.write_fitting_parameters(M_structures, Cluster_rules, J_rules, Js, intercept, 200)
cfp.plot_data3(M_structures,Cluster_rules,J_rules,Js,intercept,200)

#print('#######################\n')
##print(cp.CV_score(Js,M_structures))
##print(cp.CV_score2(M_structures))
#print('#######################\n')
##--------------------------------------------------------------#








##--------------------------------------------------------------#
## Code for running the actual MonteCarlo simulation.
## First create and initialize an array of atom objects for the simulation
## then run MonteCarlo.
#temp_data = open('Temp_data','a')
#temp_data.write('size,Temp,passes,H_avg,mimj_avg,mag_avg,absmag_avg,phase_avg,absphase_avg\n')
#temp_data.close()

#x_pts = 4 #|
#y_pts = 4 #|- Dimensions of the mc_supercellObj (simulation supercell)
#z_pts = 4 #|
#phase_init = 'mart' # initial phase configuration
#spin_init = 'FM' # initial spin configuration
#species_init = 'ordered'
#num_passes = 5 # number of cluster/wolf moves done
#num_sub_passes = 5 # number of spin/species flips done per cluster/wolf move
#Temp0 = 100 # initial temperature in K
#TempF = 100 # final temperature in K
#Temp_inc = 100 # temperature increase per pass in K

## Initialize an array of atoms with ms.mc_supercellObj(size,species,composition)
## size is (x,y,z)dimensions, species is types of atoms allowed (0=Ni,1=Mn,2=In)
## composition is number of each atom (#Ni,#Mn,#In)
#lattice = ms.mc_supercellObj((x_pts,y_pts,z_pts),(0,1,2),[32,32,0],phase_init,spin_init,species_init)
##sys.setrecursionlimit(lattice.num_sites+2)
## To actually run the simulation use
## mc.run_montecarlo(reference_to_atom_array,number_of_passes,starting_temp, BEG_rules,Cluster_rules,J_rules,plot_figs=TRUE)
## BEG_rules,Cluster_rules,J_rules are objects that determine when and how the fitted parameters are applied
#print("Beginning MonteCarlo\n")
##mc.run_WA_MCA(lattice,num_passes,num_sub_passes,Temp0,Temp_inc,TempF,Cluster_rules,J_rules,Js,do_figs=True)
#mc2.run_WA_MCA(lattice,num_passes,num_sub_passes,Temp0,Temp_inc,TempF,Cluster_rules,J_rules,Js,do_figs=True)
