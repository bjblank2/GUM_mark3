__author__ = 'brian'
import calc_params as Cp
import mc_functions as mc
import mc_functions_2 as mc2
import mc_supercell as ms
import sys


def read_inputs(file_path):
    x_pts = 4 #|
    y_pts = 4 #|- Dimensions of the mc_supercellObj (simulation supercell)
    z_pts = 8 #|
    phase_init = 'mart' # initial phase configuration
    spin_init = 'rand' # initial spin configuration
    species_init = 'rand' # inital species configuration 
    num_passes = 100 # number of cluster/wolf moves done
    num_sub_passes = 100 # number of spin/species flips done per cluster/wolf move
    Temp0 = 200 # initial temperature in K
    Temp_inc = 0 # temperature increase per pass in K
    species_flips = False
    comp = [0]*3
    input_file = open(file_path+'/starting_params','r')
    lines = input_file.readlines()
    for i in range(len(lines)):
        line = lines[i]
        line = line.split()
        if "x_" in line[0]:
            x_pts = int(line[1])
        if "y_" in line[0]:
            y_pts = int(line[1])
        if "z_" in line[0]:
            z_pts = int(line[1])
        if "phase" in line[0]:
            phase_init = line[1]
        if "spin" in line[0]:
            spin_init = line[1]
        if "num_passes" in line[0]:
            num_passes = int(line[1])
        if "num_sub_passes" in line[0]:
            num_sub_passes = int(line[1])
        if "Temp0" in line[0]:
            Temp0 = float(line[1])
        if "Temp_inc" in line[0]:
            Temp_inc = float(line[1])
        if "Comp" in line[0]:
            comp[0] = int(line[1])
            comp[1] = int(line[2])
            comp[2] = int(line[3])
        if "species" in line[0]:
            species_init = line[1]
        if "allow_element_flips" in line[0]:
            if 'T' in line[1]:
                species_flips = True
            if 'F' in line[1]:
                species_flips = False
    return x_pts,y_pts,z_pts,phase_init,spin_init,species_init,num_passes,num_sub_passes,Temp0,Temp_inc,comp,species_flips


# This is the main file that calculates the fit and runs the
# MonteCarlo code
#--------------------------------------------------------------#
# Code for finding the fitting parameters.
# First read in the DFT data, generate the appropriate data
# structures and finally write the final parameters to the
# the "Js" list
if __name__ in '__main__':
    root_dir = '/Volumes/TOURO/Ni-Fe-Ga/Data_Pts'
    data_file = '../../../../../../../../../../../home/bjblank2/GUM_mark3/NiMnIn_Data'
    beg_file = '../../../../../../../../../../../home/bjblank2/GUM_mark3/BEG_rules'
    cluster_file = '../../../../../../../../../../../home/bjblank2/GUM_mark3/Cluster_Rules'
    j_file = '../../../../../../../../../../../home/bjblank2/GUM_mark3/J_Rules'
    num_comps = 3
    num_species = 3
    # Determine if rules and structure data is defined
    Data_file_exists = True
    BEG_file_exists = True
    Cluster_rules_exist = True
    J_ruels_exist = True
    Js_exist = True
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

#    x_pts = 4 #|
#    y_pts = 4 #|- Dimensions of the mc_supercellObj (simulation supercell)
#    z_pts = 8 #|
#    phase_init = 'mart' # initial phase configuration
#    spin_init = 'rand' # initial spin configuration
#    num_passes = 2 # number of cluster/wolf moves done
#    num_sub_passes = 100 # number of spin/species flips done per cluster/wolf move
#    Temp0 = 200 # initial temperature in K
#    Temp_inc = 10 # temperature increase per pass in K

    x_pts,y_pts,z_pts,phase_init,spin_init,species_init,num_passes,num_sub_passes,Temp0,Temp_inc,comp,species_flips = read_inputs(sys.argv[1])

    # Initialize an array of atoms with ms.mc_supercellObj(size,species,composition)
    # size is (x,y,z)dimensions, species is types of atoms allowed (0=Ni,1=Mn,2=In)
    # composition is number of each atom (#Ni,#Mn,#In)
    lattice = ms.mc_supercellObj((x_pts,y_pts,z_pts),(0,1,2),comp,phase_init,spin_init,species_init)#(64,48,16))
    #sys.setrecursionlimit(lattice.num_sites+2)
    # To actually run the simulation use
    # mc.run_montecarlo(reference_to_atom_array,number_of_passes,starting_temp, BEG_rules,Cluster_rules,J_rules,plot_figs=TRUE)
    # BEG_rules,Cluster_rules,J_rules are objects that determine when and how the fitted parameters are applied
    print("Beginning MonteCarlo\n")
    #mc.run_WA_MCA(lattice,num_passes,num_sub_passes,Temp0,Temp_inc,Cluster_rules,J_rules,Js,do_figs=True)
    mc2.run_WA_MCA(lattice,num_passes,num_sub_passes,Temp0,Temp_inc,Cluster_rules,J_rules,Js,species_flips,do_figs=True)
