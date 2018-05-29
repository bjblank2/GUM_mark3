__author__ = 'brian'
# This file is the meat of the MonteCarlo code. The functions in
# this file opperate on MStructureObj objects. These objects contain
# information about the supercell as well as an array of mc_siteObj objects
# that represent the individual atomic sites and all the properties
# they can have (position,spin,species,phase,neighbors).
# These classes are defined in mc_structure.py

# The functions that I am currently using for the MonteCarlo have a #-#
# in front of their comment.

import numpy as np
cimport numpy as np
from mc_supercellCY cimport mc_supercellObj
import mpmath as math
import matplotlib.pyplot as plt
#cimport matplotlib.pyplot as plt
import mc_supercell as mcs
from copy import deepcopy
from mpl_toolkits.mplot3d import Axes3D



cpdef void do_nothing():
    print('test')

cdef list calc_BEG_params(list site,mc_supercellObj supercell,Cluster_rules,J_rules,list Js,float T):
    cdef float BEG_J = 0
    cdef float BEG_K = 0
    cdef float Kb = .000086173324
    cdef int neighbor
    cdef int neighbor_spin
    cdef int rule
    cdef mc_supercellObj supercell_obj = <mc_supercellObj>supercell
    cdef int home_spin = supercell_obj.get_site_spin(site)

    for neighbor in range(supercell_obj.get_number_of_neighbors(site)):
        neighbor_spin = supercell_obj.get_neighbor_spin(site,neighbor)
        for rule in range(len(Cluster_rules)):
            if Cluster_rules[rule].neighbor_order != 0:
                if supercell_obj.get_site_species(site) in Cluster_rules[rule].home_atom_list:
                    if supercell_obj.get_neighbor_order(site,neighbor) == Cluster_rules[rule].neighbor_order:
                        if supercell_obj.get_neighbor_plain(site,neighbor) == Cluster_rules[rule].plane or Cluster_rules[rule].plane == 'ALL':
                            if supercell_obj.get_neighbor_species(site,neighbor) in Cluster_rules[rule].neighbor_atom_list:
                                if Cluster_rules[rule].neighbor_arrangement == 'PERM':
                                    if supercell_obj.get_site_species(site) != supercell_obj.get_neighbor_species(site,neighbor):
                                        if Cluster_rules[rule].phase == 'mart':
                                            BEG_J += Js[rule]
                                        if Cluster_rules[rule].phase == 'aust':
                                            BEG_K += Js[rule]
                                if Cluster_rules[rule].neighbor_arrangement == 'COMB':
                                    if Cluster_rules[rule].phase == 'mart':
                                        BEG_J += Js[rule]
                                    if Cluster_rules[rule].phase == 'aust':
                                        BEG_K += Js[rule]
        for rule in range(len(J_rules)):
            if J_rules[rule].neighbor_order != 0:
                if supercell_obj.get_site_species(site) in J_rules[rule].home_atom_list:
                    if supercell_obj.get_neighbor_order(site,neighbor) == J_rules[rule].neighbor_order:
                        if supercell_obj.get_neighbor_plain(site,neighbor) == J_rules[rule].plane or J_rules[rule].plane == 'ALL':
                            if supercell_obj.get_neighbor_species(site,neighbor) in J_rules[rule].neighbor_atom_list:
                                if J_rules[rule].neighbor_arrangement == 'PERM':
                                    if supercell_obj.get_site_species(site) != supercell_obj.get_neighbor_species(site,neighbor):
                                        if J_rules[rule].phase == 'mart':
                                            BEG_J += Js[rule+len(Cluster_rules)]*home_spin*neighbor_spin
                                        if J_rules[rule].phase == 'aust':
                                            BEG_K += Js[rule+len(Cluster_rules)]*home_spin*neighbor_spin
                                if J_rules[rule].neighbor_arrangement == 'COMB':
                                    if J_rules[rule].phase == 'mart':
                                        BEG_J += Js[rule+len(Cluster_rules)]*home_spin*neighbor_spin
                                    if J_rules[rule].phase == 'aust':
                                        BEG_K += Js[rule+len(Cluster_rules)]*home_spin*neighbor_spin
    for rule in range(len(Cluster_rules)):
        if Cluster_rules[rule].neighbor_order == 0:
            if supercell_obj.get_site_species(site) in Cluster_rules[rule].home_atom_list:
                if Cluster_rules[rule].phase == 'mart':
                    BEG_J += Js[rule]
                if Cluster_rules[rule].phase == 'aust':
                    BEG_K += Js[rule]
    for rule in range(len(J_rules)):
        if J_rules[rule].neighbor_order == 0:
            if supercell_obj.get_site_species(site) in Cluster_rules[rule].home_atom_list:
                if J_rules[rule].phase == 'mart':
                    BEG_J += Js[rule+len(Cluster_rules)]
                if J_rules[rule].phase == 'aust':
                    BEG_K += Js[rule+len(Cluster_rules)]
    return [BEG_J,BEG_K]


#-# Determine the total energy of the entire lattice and return that energy
### COMMENT FROM ELIF: IS THIS THE ENTIRE LATTICE OR IS THIS A GIVEN SITE SPECIFIC CONTRIBUTION TO THE ENERGY??????

cdef float eval_site_new(list site, mc_supercellObj supercell, Cluster_rules,J_ruels, list Js, float T):
    cdef mc_supercellObj supercell_obj
    cdef float Kb = .000086173324
    cdef float total_Ham = 0
    cdef int site_phase
    cdef int neighbor_phase
    cdef list BEG_params
    cdef float J
    cdef float K
    cdef int neighbor

    supercell_obj = <mc_supercellObj>supercell
    site_phase = supercell_obj.get_site_phase(site)
    BEG_params = calc_BEG_params(site,supercell_obj,Cluster_rules,J_ruels,Js,T)
    J = BEG_params[0]
    K = BEG_params[1]
    for neighbor in range(supercell_obj.get_number_of_neighbors(site)):
        if supercell_obj.get_neighbor_order(site,neighbor) == 1:
            neighbor_phase = supercell_obj.get_neighbor_phase(site,neighbor)
            total_Ham += (J*(site_phase*neighbor_phase)+K*(1-site_phase**2)*(1-neighbor_phase**2))/8
    total_Ham += Kb*T*np.log(8)*(site_phase**2)
    return total_Ham


#-# Determine the total energy of the entire lattice and return that energy
cdef list eval_lattice_new(mc_supercellObj supercell, Cluster_rules, J_rules, list Js, float T):
    cdef float total_Ham = 0
    cdef float total_phase = 0
    cdef float total_phase2 = 0
    cdef float total_spin = 0
    cdef float total_spin2 = 0
    cdef mc_supercellObj supercell_obj
    cdef int i
    cdef int j
    cdef int k
    cdef list site

    ######### START ELIF COMMENT #############
    # I am assuming that this is now summing over all lattice points and each lattice point neighbors, so that each interaction does get summed twice???
    ######### END ELIF COMMENT #############

    supercell_obj = <mc_supercellObj>supercell
    for i in range(supercell_obj.i_length):
        for j in range(supercell_obj.j_length):
            for k in range(supercell_obj.k_length):
                site = [i,j,k]
                total_Ham += eval_site_new(site,supercell_obj,Cluster_rules,J_rules,Js,T)
                total_phase += supercell_obj.get_site_phase(site)
                total_phase2 += np.absolute(supercell_obj.get_site_phase(site))
                total_spin += supercell_obj.get_site_spin(site)
                total_spin2 += np.absolute(supercell_obj.get_site_spin(site))
    return [total_Ham,total_phase/supercell_obj.num_sites,total_phase2/supercell_obj.num_sites,total_spin/supercell_obj.num_sites,total_spin2/supercell_obj.num_sites]

#-# Randomly change the phase of a specific element in the lattice and return the value of the original phase
cdef list flip_phase(list site, int neighbor,mc_supercellObj supercell):
    cdef mc_supercellObj supercell_obj
    cdef int old_neighbor_phase
    cdef int old_phase
    cdef bint phase_changed
    cdef float rand
    cdef int phase

    supercell_obj = <mc_supercellObj>supercell
    old_neighbor_phase = supercell_obj.get_neighbor_phase(site,neighbor)
    old_phase = supercell_obj.get_site_phase(site)
    phase_changed = False
    while phase_changed == False:
        rand = np.random.random()
        if rand <= 1/3.0:
            phase = 0
        elif rand > 1/3.0 and rand <= 2/3.0:
            phase = -1
        elif rand >2/3.0:
            phase = 1
        if phase != old_phase:
            phase_changed = True
    supercell_obj.set_site_phase(site,phase)
    supercell_obj.set_neighbor_phase(site,neighbor,phase)
    return [old_phase, old_neighbor_phase]

#-# Randomly change the species of a specific element in the lattice and return the value of the original species
cdef list flip_species(list site_1, list site_2, mc_supercellObj supercell):
    cdef mc_supercellObj supercell_obj
    cdef int old_species_1
    cdef int old_species_2

    supercell_obj = <mc_supercellObj>supercell
    old_species_1 = supercell_obj.get_site_species(site_1)
    old_species_2 = supercell_obj.get_site_species(site_2)
    supercell_obj.set_site_species(site_1,old_species_2)
    supercell_obj.set_site_species(site_2,old_species_1)
    return [old_species_1,old_species_2]

#-# Randomly change the spin of a specific element in the lattice and return the value of the original spin
cdef int flip_spin(list site,mc_supercellObj supercell):
    cdef mc_supercellObj supercell_obj
    cdef int old_spin
    cdef bint spin_changed
    cdef float rand
    cdef int spin

    supercell_obj = <mc_supercellObj>supercell
    old_spin = supercell_obj.get_site_spin(site)
    spin_changed = False
    while spin_changed == False:
        rand = np.random.random()
        if rand <= 1/3.0:
            spin = 0
        elif rand > 1/3.0 and rand <= 2/3.0:
            spin = -1
        elif rand >2/3.0:
            spin = 1
        if spin != old_spin:
            spin_changed = True
    supercell_obj.set_site_spin(site,spin)
    return old_spin

#-# ELIF: THIS ALSO APPEARS TO BE USED
cdef int get_new_phase(list site, mc_supercellObj supercell):
    cdef mc_supercellObj supercell_obj
    cdef int old_phase
    cdef bint phase_changed
    cdef float rand
    cdef int phase

    supercell_obj = <mc_supercellObj>supercell
    old_phase = supercell_obj.get_site_phase(site)
    phase_changed = False
    while phase_changed == False:
        rand = np.random.random()
        if rand <= 1/3.0:
            phase = -1
        elif rand > 1/3.0 and rand <= 2/3.0:
            phase = 1
        elif rand >2/3.0:
            phase = 0
        if phase != old_phase:
            phase_changed = True
    return phase

#-# ELIF: THIS ALSO APPEARS TO BE USED
cdef float calc_avg_spin(list site, mc_supercellObj supercell):
    cdef mc_supercellObj supercell_obj
    cdef float M
    cdef int count
    cdef int i
    cdef int neighbor

    supercell_obj = <mc_supercellObj>supercell
    M = 0
    count = 0
    for i in range(supercell_obj.get_number_of_neighbors(site)):
        neighbor = i
        if supercell_obj.get_neighbor_order(site,neighbor) == 2:
            M += supercell_obj.get_site_spin(site)*supercell_obj.get_neighbor_spin(site,neighbor)
            count += 1
    return M/count

#-# Runs the Wolf/Mixed Cluster Algorithm
######### START ELIF COMMENT #############
# So this seems to potentially increment T after only a single attempt at growing a cluster?
# sub_passses seems to be some conditioning/single site flip tests?
# then we loop through numb_passes, but increment T each time.
# weird.
# I changed x-axis from "T" to "passes"
######### END ELIF COMMENT #############

cpdef void run_WA_MCA(mc_supercellObj supercell, int numb_passes, int num_sub_passes,float temp,float temp_inc,float tempf,Cluster_rules,J_rules,list Js, bint do_figs=True):
    cdef mc_supercellObj supercell_obj
    cdef float T = temp
    cdef float Kb = .000086173324
    cdef float M = 0
    cdef float H_total,total_phase,total_phase2,total_spin,total_spin2,old_ham,new_ham,rand,prob,H_cluster_old
    cdef float H_cluster_new, X_axis
    cdef int inc_down,inc_up,inc_not = 0
    cdef int passes,i,j,k,new_spin,old_site_species,old_randsite_species,seed_phase,new_phase
    cdef list ghost_Js,site,random_site,cluster,seed
    cdef str c
    cdef bint random_site_not_0
    cdef int neighbor
    cdef list neighbor_site
    cdef int old_neighbor_species
    cdef int site_species

    supercell_obj = <mc_supercellObj>supercell
    ghost_Js = apply_diffusion_ghost_field(2,Cluster_rules,J_rules,Js)
    H_total,total_phase,total_phase2,total_spin,total_spin2 = eval_lattice_new(supercell_obj,Cluster_rules,J_rules,Js,T)

    while T<=tempf:
        print('CURRENT TEMP = ',T)
        print('starting details of phase: total phase = ',total_phase,' ; total |phase| = ',total_phase2)
        print('starting details of magnetization: total spin = ',total_spin,' ; total |spin| = ',total_spin2,'. energy = ',H_total,'\n' )
        for passes in range(numb_passes):
            #Flip spins and Species
            print('initiating pass no. ',passes,'\n')
            print('...start subpasses')
            for sub_passes in range(num_sub_passes):
                M = 0
                for i in range(supercell_obj.i_length):
                    for j in range(supercell_obj.j_length):
                        for k in range(supercell_obj.k_length):
                            site = [i,j,k]
                            old_Ham = eval_site_new(site,supercell_obj,Cluster_rules,J_rules,Js,T)
                            old_spin = flip_spin(site,supercell_obj)
                            new_Ham = eval_site_new(site,supercell_obj,Cluster_rules,J_rules,Js,T)
                            if new_Ham < old_Ham:
                                inc_down += 1
                            else:
                                rand = np.random.random()
                                prob = math.exp(-1/(Kb*T)*(new_Ham-old_Ham))
                                if rand < prob:
                                    inc_up += 1
                                else:
                                    supercell_obj.set_site_spin(site,old_spin)
                                    inc_not += 1
                            M += calc_avg_spin(site,supercell_obj)
                            ##############
                            # FLIP SPECIES
                            ##############
                            site_species = supercell_obj.get_site_species(site)
                            if site_species != 0:
                                for neighbor in range(supercell_obj.get_number_of_neighbors(site)):
                                    if supercell_obj.get_neighbor_order(site,neighbor) != 0:
                                        if supercell_obj.get_neighbor_order(site,neighbor) != 1:
                                            if supercell_obj.get_neighbor_species(site,neighbor) != site_species:
                                                neighbor_site = supercell_obj.get_neighbor_pos(site,neighbor)
                                                old_Ham = eval_site_new(site,supercell_obj,Cluster_rules,J_rules,ghost_Js,T)
                                                old_Ham += eval_site_new(neighbor_site,supercell_obj,Cluster_rules,J_rules,ghost_Js,T)
                                                old_site_species,old_neighbor_species = flip_species(site,neighbor_site,supercell_obj)
                                                new_Ham = eval_site_new(site,supercell_obj,Cluster_rules,J_rules,ghost_Js,T)
                                                new_Ham += eval_site_new(neighbor_site,supercell_obj,Cluster_rules,J_rules,ghost_Js,T)
                                                if new_Ham < old_Ham:
                                                    inc_down += 1
                                                else:
                                                    rand = np.random.random()
                                                    prob = math.exp(-1/(Kb*T)*(new_Ham-old_Ham))
                                                    if rand < prob:
                                                        inc_up += 1
                                                    else:
                                                        supercell_obj.set_site_species(site,old_site_species)
                                                        supercell_obj.set_site_species(neighbor_site,old_neighbor_species)
                                                        inc_not += 1
#                            if supercell_obj.get_site_species(site) != 0:
#                                random_site_not_0 = False
#                                while random_site_not_0 != True:
#                                    random_site_not_0 = False
#                                    random_site = [np.random.randint(0,supercell_obj.i_length),np.random.randint(0,supercell_obj.j_length),np.random.randint(0,supercell_obj.k_length)]
#                                    if supercell_obj.get_site_species(random_site) != 0:
#                                        random_site_not_0 = True
#                                old_Ham = eval_site_new(site,supercell_obj,Cluster_rules,J_rules,ghost_Js,T)
#                                old_Ham += eval_site_new(random_site,supercell_obj,Cluster_rules,J_rules,ghost_Js,T)
#                                old_site_species,old_randsite_species = flip_species(site,random_site,supercell_obj)
#                                new_Ham = eval_site_new(site,supercell_obj,Cluster_rules,J_rules,ghost_Js,T)
#                                new_Ham += eval_site_new(random_site,supercell_obj,Cluster_rules,J_rules,ghost_Js,T)
#                                if new_Ham < old_Ham:
#                                    inc_down += 1
#                                else:
#                                    rand = np.random.random()
#                                    prob = math.exp(-1/(Kb*T)*(new_Ham-old_Ham))
#                                    if rand < prob:
#                                        inc_up += 1
#                                    else:
#                                        supercell_obj.set_site_species(site,old_site_species)
#                                        supercell_obj.set_site_species(random_site,old_randsite_species)
#                                        inc_not += 1
            H_total,total_phase,total_phase2,total_spin,total_spin2 = eval_lattice_new(supercell_obj,Cluster_rules,J_rules,Js,T)
            print('details of phase: total phase = ',total_phase,' ; total |phase| = ',total_phase2)
            print('details of magnetization: total spin = ',total_spin,' ; total |spin| = ',total_spin2,'. energy = ',H_total,'\n' )

            #Randdom Seed
            print('...sub-passes done, start cluster growth!')
            cluster = []
            seed =[np.random.randint(0,supercell_obj.i_length),np.random.randint(0,supercell_obj.j_length),np.random.randint(0,supercell_obj.k_length)]
            seed_phase = supercell_obj.get_site_phase(seed)
            new_phase = get_new_phase(seed,supercell_obj)
            grow_cluster(seed,supercell_obj,seed_phase,new_phase,cluster,Cluster_rules,J_rules,Js,T)
            ### Track size here (print(len(cluster))
            print('\tcluster length = ',len(cluster))
            if seed_phase*new_phase == -1:
                print('\tenter Wolff')
                flip_cluster(supercell_obj,seed_phase,new_phase,cluster)
                print('\taccepting Wolff cluster flip')
            else:
                print('\tenter Mixed Cluster')
                H_cluster_old = eval_cluster(supercell_obj,seed_phase,new_phase,cluster,Cluster_rules,J_rules,Js,T)
                flip_cluster(supercell_obj,seed_phase,new_phase,cluster)
                H_cluster_new = eval_cluster(supercell_obj,seed_phase,new_phase,cluster,Cluster_rules,J_rules,Js,T)
                print('new Ham = ',H_cluster_new,' ; old Ham = ',H_cluster_old)
                if H_cluster_new <= H_cluster_old:
                    inc_down += 1
                    print('accepting MC cluster flip: new energy < old energy')
                else:
                    rand = np.random.random()
                    prob = math.exp(-1/(Kb*T)*(H_cluster_new-H_cluster_old))
                    if rand < prob:
                        print('accepting MC cluster flip: prob is ',prob,' ... rand is ',rand)
                        inc_up += 1
                    else:
                        print('rejecting MC cluster flip: prob is ',prob,' ... rand is ',rand)
                        flip_cluster(supercell_obj,new_phase,seed_phase,cluster)
                        inc_not += 1

            H_total,total_phase,total_phase2,total_spin,total_spin2 = eval_lattice_new(supercell_obj,Cluster_rules,J_rules,Js,T)
            print('details of phase: total phase = ',total_phase,' ; total |phase| = ',total_phase2)
            print('details of magnetization: total spin = ',total_spin,' ; total |spin| = ',total_spin2,'. energy = ',H_total,'\n' )

            print('...finish cluster moves, run subpasses')
            for sub_passes in range(num_sub_passes):
                M = 0
                for i in range(supercell_obj.i_length):
                    for j in range(supercell_obj.j_length):
                        for k in range(supercell_obj.k_length):
                            site = [i,j,k]
                            old_Ham = eval_site_new(site,supercell_obj,Cluster_rules,J_rules,Js,T)
                            old_spin = flip_spin(site,supercell_obj)
                            new_Ham = eval_site_new(site,supercell_obj,Cluster_rules,J_rules,Js,T)
                            if new_Ham < old_Ham:
                                inc_down += 1
                            else:
                                rand = np.random.random()
                                prob = math.exp(-1/(Kb*T)*(new_Ham-old_Ham))
                                if rand < prob:
                                    inc_up += 1
                                else:
                                    supercell_obj.set_site_spin(site,old_spin)
                                    inc_not += 1
                            M += calc_avg_spin(site,supercell_obj)

            H_total,total_phase,total_phase2,total_spin,total_spin2 = eval_lattice_new(supercell_obj,Cluster_rules,J_rules,Js,T)
            print('details of phase: total phase = ',total_phase,' ; total |phase| = ',total_phase2)
            print('details of magnetization: total spin = ',total_spin,' ; total |spin| = ',total_spin2,'\n' )

            temp_output = open('Temp_data','a')
            temp_output.write(str(supercell_obj.i_length)+','+str(T)+','+str(passes)+','+str(H_total)+','+str(M/supercell_obj.num_sites)+','+str(total_spin)+','+str(total_spin2)+','+str(total_phase)+','+str(total_phase2)+'\n')
            temp_output.close()

        T += temp_inc

        if temp_inc == 0:
            X_axis = passes
        else: X_axis = T
        if supercell_obj.get_site_phase([0,0,0]) == 0:
            c = 'r'
        else: c = 'b'
        H_total,total_phase,total_phase2,total_spin,total_spin2 = eval_lattice_new(supercell_obj,Cluster_rules,J_rules,Js,T)
        plt.figure(2)
        plt.errorbar(X_axis,H_total,lw=3,marker='o',color=c)
        plt.figure(3)
        plt.subplot(311)
        plt.xlabel("Temp", fontsize=10)
        plt.ylabel("Magnetic Order Parameter", fontsize=10)
        plt.errorbar(X_axis,M/supercell_obj.num_sites,lw=3,marker='o',color=c)
        plt.subplot(312)
        plt.xlabel("Temp", fontsize=10)
        plt.ylabel("(Average Mag)^2", fontsize=10)
        plt.errorbar(X_axis,total_spin2,lw=3,marker='o',color=c)
        plt.figure(4)
        plt.subplot(411)
        plt.xlabel("Temp", fontsize=10)
        plt.ylabel("Average Phase", fontsize=10)
        plt.errorbar(X_axis,total_phase,lw=3,marker='o',color=c)
        plt.subplot(412)
        plt.xlabel("Temp", fontsize=10)
        plt.ylabel("Average Phase^2", fontsize=10)
        plt.errorbar(X_axis,total_phase2,lw=3,marker='o',color=c)

    if do_figs is True:
        plt.figure(2)
        plt.xlabel("Temp", fontsize=20)
        plt.ylabel("Energy of lattice (eV)", fontsize=20)
        plt.savefig('Enrg.pdf')
        plt.figure(3)
        plt.savefig('Mag.pdf')
        plt.figure(4)
        plt.savefig('Phase.pdf')

    fig = plt.figure(5)
    ax = fig.add_subplot(111, projection='3d')
    xs = []
    ys = []
    zs = []
    cs = []
    us = []
    vs = []
    ws = []
    for i in range(supercell_obj.i_length):
        for j in range(supercell_obj.j_length):
            for k in range(supercell_obj.k_length):
                if np.mod(k,2) == 0:
                    offset = 0
                else:
                    offset = .5
                site = [i,j,k]
                pos = supercell_obj.get_site_pos(site)
                xs.append(pos[0]+offset)
                ys.append(pos[1]+offset)
                zs.append(pos[2]*.5)
                us.append(0)
                vs.append(0)
                ws.append(supercell_obj.get_site_spin(site))
                if supercell_obj.get_site_species(site) == 0:
                    cs.append('g')
                if supercell_obj.get_site_species(site) == 1:
                    cs.append('r')
                if supercell_obj.get_site_species(site) == 2:
                    cs.append('b')
    ax.quiver(xs,ys,zs,us,vs,ws,pivot='middle',length=.5)
    ax.scatter(xs,ys,zs,c=cs,marker='o',s=50)
    plt.savefig('3D_plt.png')
    plt.show()

#-# Grows the clusters
cdef void grow_cluster(list site, mc_supercellObj supercell, int seed_phase,int new_phase,list links, Cluster_rules, J_rules, list Js, float T): # Recursive function
    cdef mc_supercellObj supercell_obj
    supercell_obj = <mc_supercellObj>supercell

    cdef float Kb = .000086173324
    cdef float B = 1/(Kb*T)
    cdef int site_phase = supercell_obj.get_site_phase(site)
    cdef list BEG_params
    cdef float J
    cdef float K
    cdef float BEG_K
    cdef float BEG_M
    cdef int neighbor
    cdef float rand
    cdef float prob
    cdef list new_site

    BEG_params = calc_BEG_params(site,supercell_obj,Cluster_rules,J_rules,Js,T)
    BEG_params[0] = J
    BEG_params[1] = K
######### START ELIF MODIFICATIONS #############
# I think here is where we need to translate the Entel terminology into the cluster algorithms terminology.
# I think Brian did not switch the J,K to K,M.  I will do it including this swap so I don't go crazy.
    #BEG_J = 2*B*J
    #BEG_K = 2*B*K
    BEG_K = -2*B*J
    BEG_M = -2*B*K
######### END ELIF MODIFICATIONS #############
    links.append(site)
    #print('site = ',site)
    # Wolff Algorithm
    if new_phase*seed_phase == -1:
        for neighbor in range(supercell_obj.get_number_of_neighbors(site)):
            #print('neighbors are ',range(supercell_obj.get_number_of_neighbors(site)))
            if supercell_obj.get_neighbor_order(site,neighbor) == 1:
                #print('here is a first neighbor')
                if supercell_obj.get_neighbor_phase(site,neighbor) == seed_phase:           ## NOT SURE WHY WE HAVE THIS HERE!!!!!
                    if supercell_obj.get_neighbor_pos(site,neighbor) not in links:
                        rand = np.random.random()
######### START ELIF MODIFICATIONS #############
                        #prob = 1-np.exp(2*BEG_J)
                        prob = 1-np.exp(-2*BEG_K)
######### START ELIF MODIFICATIONS #############
                        if rand <= prob:
                            new_site = supercell_obj.get_neighbor_pos(site,neighbor)
                            grow_cluster(new_site,supercell_obj,seed_phase,new_phase,links,Cluster_rules,J_rules,Js,T)
    # Mixed Cluster Algorithm
    if [seed_phase,new_phase] == [1,0] or [seed_phase,new_phase] == [0,-1]:
        for neighbor in range(supercell_obj.get_number_of_neighbors(site)):
            if supercell_obj.get_neighbor_order(site,neighbor) == 1:
                if supercell_obj.get_neighbor_pos(site,neighbor) not in links:
                    if supercell_obj.get_neighbor_phase(site,neighbor) == 1 or supercell_obj.get_neighbor_phase(site,neighbor) == 0:
                        if supercell_obj.get_neighbor_phase(site,neighbor) == site_phase:
                            rand = np.random.random()
######### START ELIF MODIFICATIONS #############
                            #prob = 1-np.exp(BEG_J+BEG_K/3)
                            prob = 1-np.exp(-BEG_K-BEG_M/3)
######### START ELIF MODIFICATIONS #############
                            if rand < prob:
                                new_site = supercell_obj.get_neighbor_pos(site,neighbor)
                                grow_cluster(new_site,supercell_obj,seed_phase,new_phase,links,Cluster_rules,J_rules,Js,T)
                        else:
                            rand = np.random.random()
######### START ELIF MODIFICATIONS #############
                            #prob = 1-np.exp(BEG_J-BEG_K/3)
                            prob = 1-np.exp(-BEG_K+BEG_M/3)
######### START ELIF MODIFICATIONS #############
                            if rand < prob:
                                new_site = supercell_obj.get_neighbor_pos(site,neighbor)
                                grow_cluster(new_site,supercell_obj,seed_phase,new_phase,links,Cluster_rules,J_rules,Js,T)

    if [seed_phase,new_phase] == [-1,0] or [seed_phase,new_phase] == [0,1]:
        for neighbor in range(supercell_obj.get_number_of_neighbors(site)):
            if supercell_obj.get_neighbor_order(site,neighbor) == 1:
                if supercell_obj.get_neighbor_pos(site,neighbor) not in links:
                    if supercell_obj.get_neighbor_phase(site,neighbor) == -1 or supercell_obj.get_neighbor_phase(site,neighbor) == 0:
                        if supercell_obj.get_neighbor_phase(site,neighbor) == site_phase:
                            rand = np.random.random()
######### START ELIF MODIFICATIONS #############
                            #prob = 1-np.exp(BEG_J+BEG_K/3)
                            prob = 1-np.exp(-BEG_K-BEG_M/3)
######### START ELIF MODIFICATIONS #############
                            if rand < prob:
                                new_site = supercell_obj.get_neighbor_pos(site,neighbor)
                                grow_cluster(new_site,supercell_obj,seed_phase,new_phase,links,Cluster_rules,J_rules,Js,T)
                        else:
                            rand = np.random.random()
######### START ELIF MODIFICATIONS #############
                            #prob = 1-np.exp(BEG_J-BEG_K/3)
                            prob = 1-np.exp(-BEG_K+BEG_M/3)
######### START ELIF MODIFICATIONS #############
                            if rand < prob:
                                new_site = supercell_obj.get_neighbor_pos(site,neighbor)
                                grow_cluster(new_site,supercell_obj,seed_phase,new_phase,links,Cluster_rules,J_rules,Js,T)

#-# Evaluates the total energy of the cluster
cdef float eval_cluster(mc_supercellObj supercell, int seed_phase,int new_phase,list links,Cluster_rules,J_ruels,list Js,float T):
    cdef mc_supercellObj supercell_obj
    supercell_obj = <mc_supercellObj>supercell

    cdef float Kb = .000086173324
    cdef float total_H = 0
    cdef list site
    cdef int i
    cdef int site_phase
    cdef list BEG_params
    cdef float BEG_J
    cdef float BEG_K
    cdef float total_H_inc
    cdef int inc_count
    cdef int neighbor
    cdef int neighbor_phase

    if len(links) == 1:
        site = links[0]
        total_H = eval_site_new(site,supercell_obj,Cluster_rules,J_ruels,Js,T)
    else:
        for i in range(len(links)):
            site = links[i]
            site_phase = supercell_obj.get_site_phase(site)
            BEG_params = calc_BEG_params(site,supercell_obj,Cluster_rules,J_ruels,Js,T)
            BEG_J = BEG_params[0]
            BEG_K = BEG_params[1]
            total_H_inc = 0
            inc_count = 0
            for neighbor in range(supercell_obj.get_number_of_neighbors(site)):
                if supercell_obj.get_neighbor_order(site,neighbor) == 1:
                    if supercell_obj.get_neighbor_pos(site,neighbor) in links:
                        neighbor_phase = supercell_obj.get_neighbor_phase(site,neighbor)
                        total_H_inc += BEG_J*site_phase*neighbor_phase+BEG_K*(1-site_phase**2)*(1-neighbor_phase**2)
                        inc_count += 1
            total_H += (total_H_inc/inc_count+Kb*T*np.log(8)*(site_phase**2))
    return total_H


#-# flips the cluster
cdef void flip_cluster(mc_supercellObj supercell, int seed_phase,int new_phase,list links):
    cdef mc_supercellObj supercell_obj
    supercell_obj = <mc_supercellObj>supercell

    cdef int i,old_phase

    if seed_phase*new_phase == -1:
        for i in range(len(links)):
            supercell_obj.set_site_phase(links[i],new_phase)
    else:
        for i in range(len(links)):
            if [seed_phase,new_phase] == [1,0] or [seed_phase,new_phase] == [0,-1]:
                old_phase = supercell_obj.get_site_phase(links[i])
                if old_phase == 1:
                    supercell_obj.set_site_phase(links[i],0)
                elif old_phase == 0:
                    supercell_obj.set_site_phase(links[i],-1)
            if [seed_phase,new_phase] == [-1,0] or [seed_phase,new_phase] == [0,1]:
                old_phase = supercell_obj.get_site_phase(links[i])
                if old_phase == -1:
                    supercell_obj.set_site_phase(links[i],0)
                elif old_phase == 0:
                    supercell_obj.set_site_phase(links[i],1)


cdef list apply_diffusion_ghost_field(float strength,Cluster_rules,J_ruels,list Js):
    cdef list ghost_Js = Js[:]
    cdef int i
    for i in range(len(Cluster_rules)):
        if Cluster_rules[i].neighbor_arrangement == 'COMB':
            if 0 not in Cluster_rules[i].home_atom_list:
                if 0 not in Cluster_rules[i].neighbor_atom_list:
                    if 1 not in Cluster_rules[i].neighbor_atom_list:
                        if 1 not in Cluster_rules[i].neighbor_atom_list:
                            ghost_Js[i] = ghost_Js[i]+strength
    return ghost_Js
