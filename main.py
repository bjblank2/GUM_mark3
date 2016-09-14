__author__ = 'brian'
import calc_params as Cp
import mc_functions as mc
import mc_supercell as ms
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
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

Js_r = Cp.do_robust_ls(M_structures)
# Do weighted least squares
#Cp.find_weights(M_structures, [8, 6, 4], 1)
Js = Cp.do_weighted_ls(M_structures, 200)
print(Js)
Js = Js_r
# Display data
Cp.write_data(M_structures, 200, Js)
Cp.write_output(M_structures, BEG_rules, Cluster_rules, J_rules, Js, 200)
#Cp.plot_data()
Cp.plot_data2()
#--------------------------------------------------------------#

#--------------------------------------------------------------#
T = 1001
Kb = .000086173324 #8.6173324(78)×10−5 eV*K^-1
x_pts = 2
y_pts = 2
z_pts = 4
lattice = ms.mc_supercellObj((x_pts,y_pts,z_pts),(0,1),(.5,.5))
lattice.find_neighbors()
H_total = mc.eval_supercell(lattice,BEG_rules,Cluster_rules,J_rules,Js)
print(H_total/np.size(lattice.supercell))
print(np.size(lattice.supercell))

# for i in range(4000):
#     supercell_list,supercell = mc.init_supercell(pts)
#     neighbors,neighbor_plain = mc.calc_neighbors(supercell)
#     ham,mag_ham = mc.eval_supercell(supercell_list,neighbors,neighbor_plain,BEG_rules,Cluster_rules,J_rules,Js)
#     plt.plot(i,ham/np.size(supercell),lw=3,marker='o',color='b')
# plt.show()

plt.figure(1)
plt.plot(0,H_total/np.size(lattice.supercell),lw=3,marker='o',color='b')
inc = 0
for passes in range(1,4500):
    for i in range(x_pts):
        for j in range(y_pts):
            for k in range(z_pts):
                home_site = lattice.supercell[i,j,k]
                H_new = 0
                H_old = 0
                H_old = mc.eval_site(lattice.supercell,(i,j,k),BEG_rules,Cluster_rules,J_rules,Js)
                old_home_site = mc.flip_phase(lattice.supercell,(i,j,k))
                H_new = mc.eval_site(lattice.supercell,(i,j,k),BEG_rules,Cluster_rules,J_rules,Js)
                if H_new > H_old:
                    rand = np.random.random()
                    prob = np.exp(-1/(Kb*T)*(H_new-H_old))
                    if rand > prob:
                        lattice.supercell[i,j,k] = old_home_site
                    else:
                        #H_total += H_new-H_old
                        x = 0
                else:
                    #H_total += H_new-H_old
                    x = 0
                H_new = 0
                H_old = 0
                H_old = mc.eval_site(lattice.supercell,(i,j,k),BEG_rules,Cluster_rules,J_rules,Js)
                old_home_site = mc.flip_spin(lattice.supercell,(i,j,k))
                H_new = mc.eval_site(lattice.supercell,(i,j,k),BEG_rules,Cluster_rules,J_rules,Js)
                if H_new > H_old:
                    rand = np.random.random()
                    prob = np.exp(-1/(Kb*T)*(H_new-H_old))
                    if rand > prob:
                        lattice.supercell[i,j,k] = old_home_site
                    else:
                        #H_total += H_new-H_old
                        x = 0
                else:
                    #H_total += H_new-H_old
                    x = 0
    H_total = mc.eval_supercell(lattice,BEG_rules,Cluster_rules,J_rules,Js)
    # inc +=1
    # if inc >= 100:
    #     T -= 100
    #     inc = 0
    #     if T <= 0:
    #         T = 1
    #     print(passes)
    #     print(T)
    T -= .25
    if T <= 0:
        T = 1
    inc += 1
    if inc >= 100:
        inc = 0
        print(passes)
        print(T)

    plt.figure(1)
    plt.plot(passes,H_total/np.size(lattice.supercell),lw=3,marker='o',color='b')
plt.show()

h = 0
supercell = lattice.supercell
for i in range(lattice.i_length):
    for j in range(lattice.j_length):
        for k in range(lattice.k_length):
            h_site = mc.eval_site(supercell,(i,j,k),BEG_rules,Cluster_rules,J_rules,Js)
            print(h_site)
            h += float(h_site)
print('last energy')
print(h)
