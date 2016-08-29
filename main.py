__author__ = 'brian'
import calc_params as Cp
import mc_functions as mc
import numpy as np
import matplotlib.pyplot as plt
#--------------------------------------------------------------#
root_dir = '/Volumes/TOURO/Ni-Fe-Ga/Data_Pts'
data_file = './NiMnIn_Data'
beg_file = './BEG_rules'
cluster_file = './Cluster_rules'
j_file = './J_rules'
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
# Calculate all sums
Cp.calculate_sums(M_structures, BEG_rules, Cluster_rules, J_rules)
# Do weighted least squares
Cp.find_weights(M_structures, [8, 6, 4], 1)
Js = Cp.do_weighted_ls(M_structures, 200)
# Display data
Cp.write_data(M_structures, 200, Js)
Cp.write_output(M_structures, BEG_rules, Cluster_rules, J_rules, Js, 200)
#Cp.plot_data()
#Cp.plot_data2()
#--------------------------------------------------------------#

#--------------------------------------------------------------#
pts = 4 # this refers to sub-lattice
T = 1001
Kb = .000086173324 #8.6173324(78)×10−5 eV*K^-1
supercell_list,supercell = mc.init_supercell(pts)
neighbor_list,plane_list = mc.calc_neighbors(supercell)
H_total,H_mag = mc.eval_supercell(supercell_list,neighbor_list,plane_list,BEG_rules,Cluster_rules,J_rules,Js)
print(H_total,H_mag)

# for i in range(1000):
#     supercell_list,supercell = mc.init_supercell(pts)
#     neighbors,neighbor_plain = mc.calc_neighbors(supercell)
#     ham,mag_ham = mc.eval_supercell(supercell_list,neighbors,neighbor_plain,BEG_rules,Cluster_rules,J_rules,Js)
#     plt.plot(i,ham/pts**3,lw=3,marker='o',color='b')
# plt.show()

plt.figure(1)
plt.plot(0,H_total/np.size(supercell),lw=3,marker='o',color='b')
inc = 0
for passes in range(1,1500):
    m = 0
    m2 = 0
    for i in range(len(supercell_list)):
        home_site = supercell_list[i]
        H_new = 0
        H_old = 0
        H_old,H_mag_old = mc.eval_site(supercell_list[i],supercell_list,neighbor_list[i],plane_list[i],BEG_rules,Cluster_rules,J_rules,Js)
        old_home_site = mc.flip_phase(i,supercell_list)
        H_new,H_mag_old = mc.eval_site(supercell_list[i],supercell_list,neighbor_list[i],plane_list[i],BEG_rules,Cluster_rules,J_rules,Js)
        if H_new > H_old:
            rand = np.random.random()
            prob = np.exp(-1/(Kb*T)*(H_new-H_old))
            if rand > prob:
                supercell_list[i] = old_home_site
            else:
                #H_total += H_new-H_old
                x = 0
        else:
            #H_total += H_new-H_old
            x = 0
        H_new = 0
        H_old = 0
        H_old,H_mag_old = mc.eval_site(supercell_list[i],supercell_list,neighbor_list[i],plane_list[i],BEG_rules,Cluster_rules,J_rules,Js)
        old_home_site = mc.flip_spin(i,supercell_list)
        H_new,H_mag_old = mc.eval_site(supercell_list[i],supercell_list,neighbor_list[i],plane_list[i],BEG_rules,Cluster_rules,J_rules,Js)
        if H_new > H_old:
            rand = np.random.random()
            prob = np.exp(-1/(Kb*T)*(H_new-H_old))
            if rand > prob:
                supercell_list[i] = old_home_site
            else:
                #H_total += H_new-H_old
                x = 0
        else:
            #H_total += H_new-H_old
            x = 0
        m2 += home_site[5]**2/np.size(supercell)
        m += home_site[5]/np.size(supercell)
    H_total,H_mag = mc.eval_supercell(supercell_list,neighbor_list,plane_list,BEG_rules,Cluster_rules,J_rules,Js)
    # inc +=1
    # if inc >= 100:
    #     T -= 100
    #     inc = 0
    #     if T <= 0:
    #         T = 1
    #     print(passes)
    #     print(T)
    T -= 1
    if T <= 0:
        T = 1
    inc += 1
    if inc >= 100:
        inc = 0
        print(passes)
        print(T)
        
    plt.figure(1)
    plt.plot(passes,H_total/np.size(supercell),lw=3,marker='o',color='b')
    plt.figure(2)
    plt.plot(passes,m/np.size(supercell),lw=3,marker='o',color='b')
    plt.figure(3)
    plt.plot(passes,m2/np.size(supercell),lw=3,marker='o',color='b')
plt.show()

