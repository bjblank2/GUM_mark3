__author__ = 'brian'
import calc_params as Cp
import mc_functions as mc
import mc_supercell as ms
#--------------------------------------------------------------#
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
BEG_rules = Cp.read_beg_rules(beg_file)
Cluster_rules = Cp.read_cluster_rules(cluster_file)
J_rules = Cp.read_j_rules(j_file)
M_structures = Cp.read_m_structure_data(data_file, num_species, len(BEG_rules), len(Cluster_rules), len(J_rules))
Cp.linearize(M_structures)
# Calculate all sums
Cp.calculate_sums(M_structures, BEG_rules, Cluster_rules, J_rules)

# Do weighted least squares
Cp.find_weights(M_structures, [8, 6, 4], 1)
#Cp.find_weights_2(M_structures, [8,6,4], 5)
Js_w = Cp.do_weighted_ls(M_structures, 5.00)
#Js_r = Cp.do_robust_ls(M_structures)
#Js_rancom = Cp.ransacom(M_structures,1,30,3000)
#Js_ran = Cp.ransac(M_structures,1,50,1000)
#Js = Js_r
Js = Js_w
#Js = Js_ran
#Js = Js_rancom
# Display data
Cp.write_data(M_structures, 200, Js)
Cp.write_output(M_structures, BEG_rules, Cluster_rules, J_rules, Js, 200)
#Cp.plot_data()
Cp.plot_data2()
#--------------------------------------------------------------#

#--------------------------------------------------------------#
temp_data = open('Temp_data','w')
temp_data.write('Temp  H_avg  mag_avg  mag2_avg  phase_avg  phase2_avg\n')
temp_data.close()
x_pts = 4
y_pts = 4
z_pts = 8
lattice = ms.mc_supercellObj((x_pts,y_pts,z_pts),(0,1,2),[64,64,0])#(64,48,16))
#lattice = ms.mc_supercellObj((x_pts,y_pts,z_pts),(0,1,2),(8,6,2))#(64,48,16))
lattice.find_neighbors()
mc.run_montecarlo(lattice,500,1,BEG_rules,Cluster_rules,J_rules,Js,do_figs=True)
#for i in range(150,1700,50):
#    mc.run_montecarlo(lattice,900,i,BEG_rules,Cluster_rules,J_rules,Js,do_figs=False)
#    print(i)
