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
import mpmath as math
import matplotlib as mpl
import matplotlib.pyplot as plt
import mc_supercell as mcs
from copy import deepcopy
from mpl_toolkits.mplot3d import Axes3D

#-# ELIF: I THINK WE NEED THIS ONE TOO.
def calc_BEG_params(site,supercell_obj,Cluster_rules,J_rules,Js,T):
    H_BEG_J = 0
    H_BEG_K = 0
    Kb = .000086173324
    for neighbor in range(supercell_obj.get_number_of_neighbors(site)):
        for Cluster_rule in range(len(Cluster_rules)):
            if supercell_obj.get_neighbor_order(site,neighbor) == Cluster_rules[Cluster_rule].neighbor_order:
                    if supercell_obj.get_neighbor_plain(site,neighbor) == Cluster_rules[Cluster_rule].plane or Cluster_rules[Cluster_rule].plane == 'ALL':
                        if supercell_obj.get_site_species(site) in Cluster_rules[Cluster_rule].home_atom_list:
                            if supercell_obj.get_neighbor_species(site,neighbor) in Cluster_rules[Cluster_rule].neighbor_atom_list:
                                if Cluster_rules[Cluster_rule].neighbor_arrangement == 'PERM':
                                    if supercell_obj.get_site_species(site) != supercell_obj.get_neighbor_species(site,neighbor):
                                        if Cluster_rules[Cluster_rule].phase == "mart":
                                            H_BEG_J += float(Js[Cluster_rule])
                                        if Cluster_rules[Cluster_rule].phase == "aust":
                                            H_BEG_K += float(Js[Cluster_rule])
                                if Cluster_rules[Cluster_rule].neighbor_arrangement == 'COMB':
                                    if Cluster_rules[Cluster_rule].phase == "mart":
                                        H_BEG_J += float(Js[Cluster_rule])
                                    if Cluster_rules[Cluster_rule].phase == "aust":
                                        H_BEG_K += float(Js[Cluster_rule])
        for J_rule in range(len(J_rules)):
            if supercell_obj.get_neighbor_order(site,neighbor) == J_rules[J_rule].neighbor_order:
                if supercell_obj.get_neighbor_plain(site,neighbor) == J_rules[J_rule].plane or J_rules[J_rule].plane == 'ALL':
                    if supercell_obj.get_site_species(site) in J_rules[J_rule].home_atom_list:
                        if supercell_obj.get_neighbor_species(site,neighbor) in J_rules[J_rule].neighbor_atom_list:
                            if J_rules[J_rule].neighbor_arrangement == 'PERM':
                                if supercell_obj.get_site_species(site) != supercell_obj.get_neighbor_species(site,neighbor):
                                    home_spin = supercell_obj.get_site_spin(site)
                                    neighbor_spin = supercell_obj.get_neighbor_spin(site,neighbor)
                                    if J_rules[J_rule].phase == "mart":
                                        H_BEG_J += float(Js[J_rule+len(Cluster_rules)])*home_spin*neighbor_spin
                                    if J_rules[J_rule].phase == "aust":
                                        H_BEG_K += float(Js[J_rule+len(Cluster_rules)])*home_spin*neighbor_spin
                            if J_rules[J_rule].neighbor_arrangement == 'COMB':
                                home_spin = supercell_obj.get_site_spin(site)
                                neighbor_spin = supercell_obj.get_neighbor_spin(site,neighbor)
                                if J_rules[J_rule].phase == "mart":
                                    H_BEG_J += float(Js[J_rule+len(Cluster_rules)])*home_spin*neighbor_spin
                                if J_rules[J_rule].phase == "aust":
                                    H_BEG_K += float(Js[J_rule+len(Cluster_rules)])*home_spin*neighbor_spin
######### START ELIF MODIFICATIONS #############
    #K = -1*(H_BEG_K/(-1*8)-Kb*T*np.log(2)/8)
    #delta = 8*K+Kb*T*np.log(2)
    #J = -1*((H_BEG_J-delta)/(-1*8)-K)
    K = H_BEG_K/8
    J = H_BEG_J/8
# As I understand it, we are just evaluating a site-specific J, K here.
# Not an energy.
# Consequently we do not need to add these delta and kB*T*log(2) type terms in here
######## END ELIF MODIFICATIONS  ##############################
    return J,K

#-# Determine the total energy of the entire lattice and return that energy
### COMMENT FROM ELIF: IS THIS THE ENTIRE LATTICE OR IS THIS A GIVEN SITE SPECIFIC CONTRIBUTION TO THE ENERGY??????
def eval_site_new(site,supercell_obj,Cluster_rules,J_ruels,Js,T,mag_field):
    Kb = .000086173324
    g = 2
    ub = 5.7883818012*10**(-5) # Bohr Mag
    total_Ham = 0
    site_phase = supercell_obj.get_site_phase(site)
    J,K = calc_BEG_params(site,supercell_obj,Cluster_rules,J_ruels,Js,T)
    for neighbor in range(supercell_obj.get_number_of_neighbors(site)):
        if supercell_obj.get_neighbor_order(site,neighbor) == 1:
            neighbor_phase = supercell_obj.get_neighbor_phase(site,neighbor)
            total_Ham += J*(site_phase*neighbor_phase)+K*(1-site_phase**2)*(1-neighbor_phase**2)
    total_Ham += Kb*T*np.log(8)*(site_phase**2)
    total_Ham -= g*ub*mag_field*np.sign(mag_field)*supercell_obj.get_site_spin(site)
    return total_Ham

#-# Determine the total energy of the entire lattice and return that energy
def eval_lattice_new(supercell_obj,Cluster_rules,J_rules,Js,T,mag_field):
    total_Ham = 0
    total_phase = 0
    total_phase2 = 0
    total_spin = 0
    total_spin2 = 0
######### START ELIF COMMENT #############
# I am assuming that this is now summing over all lattice points and each lattice point neighbors, so that each interaction does get summed twice???
######### END ELIF COMMENT #############
    for i in range(supercell_obj.i_length):
        for j in range(supercell_obj.j_length):
            for k in range(supercell_obj.k_length):
                site = [i,j,k]
                total_Ham += eval_site_new(site,supercell_obj,Cluster_rules,J_rules,Js,T,mag_field)
                total_phase += supercell_obj.get_site_phase(site)/supercell_obj.num_sites
                total_phase2 += supercell_obj.get_site_phase(site)**2/supercell_obj.num_sites
                total_spin += supercell_obj.get_site_spin(site)/supercell_obj.num_sites
                total_spin2 += supercell_obj.get_site_spin(site)**2/supercell_obj.num_sites
    return total_Ham,total_phase,total_phase2,total_spin,total_spin2

#-# Randomly change the phase of a specific element in the lattice and return the value of the original phase
def flip_phase(site,neighbor,supercell_obj):
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
    return old_phase, old_neighbor_phase

#-# Randomly change the species of a specific element in the lattice and return the value of the original species
def flip_species(site_1,site_2,supercell_obj):
    old_species_1 = supercell_obj.get_site_species(site_1)
    old_species_2 = supercell_obj.get_site_species(site_2)
    supercell_obj.set_site_species(site_1,old_species_2)
    supercell_obj.set_site_species(site_2,old_species_1)
    return old_species_1,old_species_2

#-# Randomly change the spin of a specific element in the lattice and return the value of the original spin
def flip_spin(site,supercell_obj):
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
def get_new_phase(site,supercell_obj):
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
def calc_avg_spin(site,supercell_obj):
    M = 0
    count = 0
    for i in range(supercell_obj.get_number_of_neighbors(site)):
        neighbor = i
        if supercell_obj.get_neighbor_order(site,neighbor) == 2:
            M += supercell_obj.get_site_spin(site)*supercell_obj.get_neighbor_spin(site,neighbor)
            count += 1
    return M/count

#-# Runs the Wolf/Mixed Cluster Algorithm
def run_WA_MCA(supercell_obj,numb_passes,num_sub_passes,temp,temp_inc,tempf,mag_field,Cluster_rules,J_rules,Js,species_flips):
    T = temp
    Kb = .000086173324
    inc_down = 0
    inc_up = 0
    inc_not = 0
    M = 0
    ghost_Js = apply_diffusion_ghost_field(2,Cluster_rules,J_rules,Js)
    H_total,total_phase,total_phase2,total_spin,total_spin2 = eval_lattice_new(supercell_obj,Cluster_rules,J_rules,Js,T,mag_field)

    while T<=tempf:
        print('\nCURRENT TEMP = ',T,'\n')
        for passes in range(numb_passes):
            #Flip spins and Species
            print('initiating pass no. ',passes)
            print('...start subpasses')
            for sub_passes in range(num_sub_passes):
                M = 0
                for i in range(supercell_obj.i_length):
                    for j in range(supercell_obj.j_length):
                        for k in range(supercell_obj.k_length):
                            site = [i,j,k]
                            old_Ham = eval_site_new(site,supercell_obj,Cluster_rules,J_rules,Js,T,mag_field)
                            old_spin = flip_spin(site,supercell_obj)
                            new_Ham = eval_site_new(site,supercell_obj,Cluster_rules,J_rules,Js,T,mag_field)
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
                            if species_flips == True:
                                if supercell_obj.get_site_species(site) != 0:
                                    random_site_not_0 = False
                                    species_not_same = False
                                    while [species_not_same,random_site_not_0] != [True,True]:
                                        random_site_not_0 = False
                                        species_not_same = False
                                        random_site = [np.random.randint(0,supercell_obj.i_length),np.random.randint(0,supercell_obj.j_length),np.random.randint(0,supercell_obj.k_length)]
                                        if supercell_obj.get_site_species(random_site) != 0:
                                            random_site_not_0 = True
                                        if supercell_obj.get_site_species(random_site) != supercell_obj.get_site_species(site):
                                            species_not_same = True
                                    old_Ham = eval_site_new(site,supercell_obj,Cluster_rules,J_rules,ghost_Js,T,mag_field)
                                    old_Ham += eval_site_new(random_site,supercell_obj,Cluster_rules,J_rules,ghost_Js,T,mag_field)
                                    old_site_species,old_randsite_species = flip_species(site,random_site,supercell_obj)
                                    new_Ham = eval_site_new(site,supercell_obj,Cluster_rules,J_rules,ghost_Js,T,mag_field)
                                    new_Ham += eval_site_new(random_site,supercell_obj,Cluster_rules,J_rules,ghost_Js,T,mag_field)
                                    if new_Ham < old_Ham:
                                        inc_down += 1
                                    else:
                                        rand = np.random.random()
                                        prob = math.exp(-1/(Kb*T)*(new_Ham-old_Ham))
                                        if rand < prob:
                                            inc_up += 1
                                        else:
                                            supercell_obj.set_site_species(site,old_site_species)
                                            supercell_obj.set_site_species(random_site,old_randsite_species)
                                            inc_not += 1


            #Randdom Seed
            print('...sub-passes done, start cluster growth!')
            cluster = []
            seed =(np.random.randint(0,supercell_obj.i_length),np.random.randint(0,supercell_obj.j_length),np.random.randint(0,supercell_obj.k_length))
            seed_phase = supercell_obj.get_site_phase(seed)
            new_phase = get_new_phase(seed,supercell_obj)
            grow_cluster(seed,supercell_obj,seed_phase,new_phase,cluster,Cluster_rules,J_rules,Js,T,mag_field)
            ### Track size here (print(len(cluster))
            #print('[seed_phase, new_phase] = ',[seed_phase,new_phase])
            print('\tcluster length = ',len(cluster))
            if seed_phase*new_phase == -1:
                print('\tenter Wolff')
                flip_cluster(supercell_obj,seed_phase,new_phase,cluster)
                print('\taccepting Wolff cluster flip')
            else:
                print('\tenter Mixed Cluster')
                H_cluster_old = eval_cluster(supercell_obj,seed_phase,new_phase,cluster,Cluster_rules,J_rules,Js,T,mag_field)
                flip_cluster(supercell_obj,seed_phase,new_phase,cluster)
                H_cluster_new = eval_cluster(supercell_obj,seed_phase,new_phase,cluster,Cluster_rules,J_rules,Js,T,mag_field)
                print('new Ham = ',H_cluster_new,' ; old Ham = ',H_cluster_old)
                if H_cluster_new <= H_cluster_old:
                    inc_down += 1
                    print('\taccepting MC cluster flip: new energy < old energy')
                else:
                    rand = np.random.random()
                    prob = math.exp(-1/(Kb*T)*(H_cluster_new-H_cluster_old))
                    #print([H_cluster_new,H_cluster_old])
                    if rand < prob:
                        print('\taccepting MC cluster flip: prob is ',prob,' ... rand is ',rand)
                        inc_up += 1
                    else:
                        print('\trejecting MC cluster flip: prob is ',prob,' ... rand is ',rand)
                        flip_cluster(supercell_obj,new_phase,seed_phase,cluster)
                        inc_not += 1

            print('...finish cluster moves, run subpasses \n')
            for sub_passes in range(num_sub_passes):
                M = 0
                for i in range(supercell_obj.i_length):
                    for j in range(supercell_obj.j_length):
                        for k in range(supercell_obj.k_length):
                            site = [i,j,k]
                            old_Ham = eval_site_new(site,supercell_obj,Cluster_rules,J_rules,Js,T,mag_field)
                            old_spin = flip_spin(site,supercell_obj)
                            new_Ham = eval_site_new(site,supercell_obj,Cluster_rules,J_rules,Js,T,mag_field)
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

            H_total,total_phase,total_phase2,total_spin,total_spin2 = eval_lattice_new(supercell_obj,Cluster_rules,J_rules,Js,T,mag_field)

            temp_output = open('Temp_data','a')
            temp_output.write(str(supercell_obj.i_length)+','+str(T)+','+str(passes)+','+str(H_total/supercell_obj.num_sites)+','+str(M/supercell_obj.num_sites)+','+str(total_spin)+','+str(total_spin2)+','+str(total_phase)+','+str(total_phase2)+'\n')
            temp_output.close()

        T += temp_inc


def grow_cluster(site,supercell_obj,seed_phase,new_phase,links,Cluster_rules,J_rules,Js,T,mag_field): # Recursive function
    Kb = .000086173324
    B = 1/(Kb*T)
    site_phase = supercell_obj.get_site_phase(site)
    J,K = calc_BEG_params(site,supercell_obj,Cluster_rules,J_rules,Js,T)
    BEG_K = -2*B*J
    BEG_M = -2*B*K
    links.append(site)
    # Wolf Algorithm
    if new_phase*seed_phase == -1:
        for neighbor in range(supercell_obj.get_number_of_neighbors(site)):
            if supercell_obj.get_neighbor_order(site,neighbor) == 1:
                if supercell_obj.get_neighbor_phase(site,neighbor) == seed_phase:
                    if supercell_obj.get_neighbor_pos(site,neighbor) not in links:
                        rand = np.random.random()
                        prob = 1-np.exp(-2*BEG_K)
                        if rand <= prob:
                            new_site = supercell_obj.get_neighbor_pos(site,neighbor)
                            grow_cluster(new_site,supercell_obj,seed_phase,new_phase,links,Cluster_rules,J_rules,Js,T,mag_field)
    # Mixed Cluster Algorithm
    if [seed_phase,new_phase] == [1,0] or [seed_phase,new_phase] == [0,-1]:
        for neighbor in range(supercell_obj.get_number_of_neighbors(site)):
            if supercell_obj.get_neighbor_order(site,neighbor) == 1:
                if supercell_obj.get_neighbor_pos(site,neighbor) not in links:
                    if supercell_obj.get_neighbor_phase(site,neighbor) == 1 or supercell_obj.get_neighbor_phase(site,neighbor) == 0:
                        if supercell_obj.get_neighbor_phase(site,neighbor) == site_phase:
                            rand = np.random.random()
                            prob = 1-np.exp(-BEG_K-BEG_M/3)
                            if rand < prob:
                                new_site = supercell_obj.get_neighbor_pos(site,neighbor)
                                grow_cluster(new_site,supercell_obj,seed_phase,new_phase,links,Cluster_rules,J_rules,Js,T,mag_field)
                        else:
                            rand = np.random.random()
                            prob = 1-np.exp(-BEG_K+BEG_M/3)
                            if rand < prob:
                                new_site = supercell_obj.get_neighbor_pos(site,neighbor)
                                grow_cluster(new_site,supercell_obj,seed_phase,new_phase,links,Cluster_rules,J_rules,Js,T,mag_field)

    if [seed_phase,new_phase] == [-1,0] or [seed_phase,new_phase] == [0,1]:
        for neighbor in range(supercell_obj.get_number_of_neighbors(site)):
            if supercell_obj.get_neighbor_order(site,neighbor) == 1:
                if supercell_obj.get_neighbor_pos(site,neighbor) not in links:
                    if supercell_obj.get_neighbor_phase(site,neighbor) == -1 or supercell_obj.get_neighbor_phase(site,neighbor) == 0:
                        if supercell_obj.get_neighbor_phase(site,neighbor) == site_phase:
                            rand = np.random.random()
                            prob = 1-np.exp(-BEG_K-BEG_M/3)
                            if rand < prob:
                                new_site = supercell_obj.get_neighbor_pos(site,neighbor)
                                grow_cluster(new_site,supercell_obj,seed_phase,new_phase,links,Cluster_rules,J_rules,Js,T,mag_field)
                        else:
                            rand = np.random.random()
                            prob = 1-np.exp(-BEG_K+BEG_M/3)
                            if rand < prob:
                                new_site = supercell_obj.get_neighbor_pos(site,neighbor)
                                grow_cluster(new_site,supercell_obj,seed_phase,new_phase,links,Cluster_rules,J_rules,Js,T,mag_field)

#-# Evaluates the total energy of the cluster
def eval_cluster(supercell_obj,seed_phase,new_phase,links,Cluster_rules,J_ruels,Js,T,mag_field):
    Kb = .000086173324
    total_H = 0
    for i in range(len(links)):
        site = links[i]
        total_H += eval_site_new(site,supercell_obj,Cluster_rules,J_ruels,Js,T,mag_field)
    return total_H


#-# flips the cluster
def flip_cluster(supercell_obj,seed_phase,new_phase,links):
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


def apply_diffusion_ghost_field(strength,Cluster_rules,J_ruels,Js):
    ghost_Js = Js[:]
    for i in range(len(Cluster_rules)):
        if Cluster_rules[i].neighbor_arrangement == 'COMB':
            if 0 not in Cluster_rules[i].home_atom_list:
                if 0 not in Cluster_rules[i].neighbor_atom_list:
                    if 1 not in Cluster_rules[i].neighbor_atom_list:
                        if 1 not in Cluster_rules[i].neighbor_atom_list:
                            ghost_Js[i] = ghost_Js[i]+strength
    return ghost_Js

