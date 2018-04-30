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
import matplotlib.pyplot as plt
import mc_supercell as mcs
from copy import deepcopy
from mpl_toolkits.mplot3d import Axes3D

#-# Determine the energy contribution from a singel atomic site and return that contribution
def eval_site(site,supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T):
        Kb = .000086173324
        Ham = 0
        Ham_BEG = 0
        Ham_Clust = 0
        Ham_Spin = 0
        for i in range(supercell_obj.get_number_of_neighbors(site)):
            BEG_J = 0
            BEG_K = 0
            for j in range(len(BEG_rules)):
                supercell_comp = [x / supercell_obj.num_sites for x in supercell_obj.get_composition()]
                BEG_comp = [x / 16 for x in BEG_rules[j].composition]
                if supercell_comp == BEG_comp:
                    if supercell_obj.get_neighbor_order(site,i) == BEG_rules[j].neighbor_order:
                       if supercell_obj.get_neighbor_plain(site,i) == BEG_rules[j].plane or BEG_rules[j].plane == 'ALL':
                            if supercell_obj.get_site_species(site) in BEG_rules[j].home_atom_list:
                                if supercell_obj.get_neighbor_species(site,i) in BEG_rules[j].neighbor_atom_list:
                                    if BEG_rules[j].neighbor_arrangement == 'PERM':
                                        if supercell_obj.get_site_species(site) != supercell_obj.get_neighbor_species(site,i):
                                            if BEG_rules[j].phase == 'mart':
                                                #BEG_J = float(Js[j])
                                                BEG_J = -2
                                            if BEG_rules[j].phase == 'aust':
                                                #BEG_K = float(Js[j])
                                                BEG_K = -1
                                    if BEG_rules[j].neighbor_arrangement == 'COMB':
                                        if BEG_rules[j].phase == 'mart':
                                            #BEG_J = float(Js[j])
                                            BEG_J = -2
                                        if BEG_rules[j].phase == 'aust':
                                            #BEG_K = float(Js[j])
                                            BEG_K = -1
            site_phase = supercell_obj.get_site_phase(site)
            neighbor_phase = supercell_obj.get_neighbor_phase(site,i)
            #Ham += -BEG_J*home_phase*neighbor_phase-BEG_K*(home_phase**2)*(neighbor_phase**2)+6*BEG_K*home_phase**2-6*supercell_obj.num_sites*BEG_K/2
            Ham += (BEG_J*site_phase*neighbor_phase+BEG_K*(1-site_phase**2)*(1-neighbor_phase**2))
            #Ham += -.0001*T*(1-site_phase**2)


            for j in range(len(Cluster_rules)):
                if supercell_obj.get_neighbor_order(site,i) == Cluster_rules[j].neighbor_order:
                    if supercell_obj.get_neighbor_plain(site,i) == Cluster_rules[j].plane or Cluster_rules[j].plane == 'ALL':
                        if supercell_obj.check_site_phase(site) == Cluster_rules[j].phase:
                            if supercell_obj.get_site_species(site) in Cluster_rules[j].home_atom_list:
                                if supercell_obj.get_neighbor_species(site,i) in Cluster_rules[j].neighbor_atom_list:
                                    if Cluster_rules[j].neighbor_arrangement == 'PERM':
                                        if supercell_obj.get_site_species(site) != supercell_obj.get_neighbor_species(site,i):
                                            Ham += float(Js[j+len(BEG_rules)])
                                            Ham_Clust += float(Js[j+len(BEG_rules)])
                                    if Cluster_rules[j].neighbor_arrangement == 'COMB':
                                        Ham += float(Js[j+len(BEG_rules)])
                                        Ham_Clust += float(Js[j+len(BEG_rules)])


            for j in range(len(J_rules)):
                if supercell_obj.get_neighbor_order(site,i) == J_rules[j].neighbor_order:
                    if supercell_obj.get_neighbor_plain(site,i) == J_rules[j].plane or J_rules[j].plane == 'ALL':
                        if supercell_obj.check_site_phase(site) == J_rules[j].phase:
                            if supercell_obj.get_site_species(site) in J_rules[j].home_atom_list:
                                if supercell_obj.get_neighbor_species(site,i) in J_rules[j].neighbor_atom_list:
                                    if J_rules[j].neighbor_arrangement == 'PERM':
                                        if supercell_obj.get_site_species(site) != supercell_obj.get_neighbor_species(site,i):
                                            home_spin = supercell_obj.get_site_spin(site)
                                            neighbor_spin = supercell_obj.get_neighbor_spin(site,i)
                                            Ham += float(Js[j+len(BEG_rules)+len(Cluster_rules)])*home_spin*neighbor_spin
                                            Ham_Spin += float(Js[j+len(BEG_rules)+len(Cluster_rules)])*home_spin*neighbor_spin
                                    if J_rules[j].neighbor_arrangement == 'COMB':
                                        home_spin = supercell_obj.get_site_spin(site)
                                        neighbor_spin = supercell_obj.get_neighbor_spin(site,i)
                                        Ham += float(Js[j+len(BEG_rules)+len(Cluster_rules)])*home_spin*neighbor_spin
                                        Ham_Spin += float(Js[j+len(BEG_rules)+len(Cluster_rules)])*home_spin*neighbor_spin
        #print([Ham_BEG,Ham_Clust,Ham_Spin])
        #print(Ham)
        site_phase = supercell_obj.get_site_phase(site)
        #Ham += -.00051*T*(1-site_phase**2)
        Ham += -Kb*T*np.log(2)*(1-site_phase**2)
        return Ham

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
    J = H_BEG_J/8  ## No KbT !!!!!
    K = H_BEG_K/8
    return J,K

#-# Determine the total energy of the entire lattice and return that energy
def eval_site_new(site,supercell_obj,Cluster_rules,J_ruels,Js,T):
    Kb = .000086173324
    total_Ham = 0
    site_phase = supercell_obj.get_site_phase(site)
    J,K = calc_BEG_params(site,supercell_obj,Cluster_rules,J_ruels,Js,T)
    for neighbor in range(supercell_obj.get_number_of_neighbors(site)):
        if supercell_obj.get_neighbor_order(site,neighbor) == 1:
            neighbor_phase = supercell_obj.get_neighbor_phase(site,neighbor)
            total_Ham += J*(site_phase*neighbor_phase)+K*(1-site_phase**2)*(1-neighbor_phase**2)
    total_Ham += Kb*T*np.log(2)*(site_phase**2)
    return total_Ham

#-# Determine the total energy of the entire lattice and return that energy
def eval_lattice_new(supercell_obj,Cluster_rules,J_rules,Js,T):
    total_Ham = 0
    total_phase = 0
    total_phase2 = 0
    total_spin = 0
    total_spin2 = 0
    for i in range(supercell_obj.i_length):
        for j in range(supercell_obj.j_length):
            for k in range(supercell_obj.k_length):
                site = [i,j,k]
                total_Ham += eval_site_new(site,supercell_obj,Cluster_rules,J_rules,Js,T)
                total_phase += supercell_obj.get_site_phase(site)/supercell_obj.num_sites
                total_phase2 += supercell_obj.get_site_phase(site)**2/supercell_obj.num_sites
                total_spin += supercell_obj.get_site_spin(site)/supercell_obj.num_sites
                total_spin2 += supercell_obj.get_site_spin(site)**2/supercell_obj.num_sites
    return total_Ham,total_phase,total_phase2,total_spin,total_spin2

def run_montecarlo_new(supercell_obj,numb_passes,numb_sub_passes,temp,Cluster_rules,J_rules,Js):
    T = temp
    Kb = .000086173324 #8.6173324(78)×10−5 eV*K^-1
    H_avg = 0
    mag_avg = 0
    mag2_avg = 0
    p_avg = 0
    p2_avg = 0
    inc_down = 0
    inc_up = 0
    inc_not = 0
    inc_T = 0
    H_down = 0
    H_up = 0
    H = 0
    H_total,p,p2,mag,mag2 = eval_lattice_new(supercell_obj,Cluster_rules,J_rules,Js,T)
    plt.figure(2)
    plt.plot(-1,H_total,lw=3,marker='o',color='b')
    plt.plot(-1,H_total,lw=3,marker='o',color='r')
    print(H_total)
    for passes in range(numb_passes):
        M_plot = 0
        mag_plot = 0
        phase_plot = 0
        phase2_plot = 0
        enrg_plot = 0
        for sub_passes in range(numb_sub_passes):
            H_total = 0
            p = 0
            p2 = 0
            mag = 0
            mag2 = 0
            M = 0
            for i in range(supercell_obj.i_length):
                for j in range(supercell_obj.j_length):
                    for k in range(supercell_obj.k_length):
                        site = [i,j,k]
                        neighbor_found = False
                        neighbor = 0
                        old_Ham = 0
                        new_Ham = 0
                        old_Ham = eval_site_new(site,supercell_obj,Cluster_rules,J_rules,Js,T)
                        neighbor_site = supercell_obj.get_neighbor_pos(site,neighbor)
                        old_Ham += eval_site_new(neighbor_site,supercell_obj,Cluster_rules,J_rules,Js,T)
                        old_phase,old_neighbor_phase = flip_phase(site,neighbor,supercell_obj)
                        new_Ham = eval_site_new(site,supercell_obj,Cluster_rules,J_rules,Js,T)
                        new_Ham += eval_site_new(neighbor_site,supercell_obj,Cluster_rules,J_rules,Js,T)
                        if new_Ham < old_Ham:
                            inc_down += 1
                        else:
                            rand = np.random.random()
                            prob = np.exp(-1/(Kb*T)*(new_Ham-old_Ham))
                            if rand > prob:
                                supercell_obj.set_site_phase(site,old_phase)
                                supercell_obj.set_neighbor_phase(site,neighbor,old_neighbor_phase)
                                inc_not += 1
                            else:
                                #H_total += (new_Ham-old_Ham)
                                inc_up += 1

                        # if supercell_obj.get_site_species(site) != 0:
                        #     composition = supercell_obj.get_composition()
                        #     if composition[2] != 0:
                        #         site_found = False
                        #         while site_found == False:
                        #             site_2 = [np.random.randint(0,supercell_obj.i_length),np.random.randint(0,supercell_obj.j_length),np.random.randint(0,supercell_obj.k_length)]
                        #             if supercell_obj.get_site_species(site_2) != supercell_obj.get_site_species(site):
                        #                 if supercell_obj.get_site_species(site_2) != 0:
                        #                     site_found = True
                        #         old_Ham = eval_site(site,supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
                        #         old_Ham += eval_site(site_2,supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
                        #         old_site_state,old_site_2_state = flip_species(site,site_2,supercell_obj)
                        #         new_Ham = eval_site(site,supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
                        #         new_Ham += eval_site(site_2,supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
                        #         if new_Ham > old_Ham:
                        #             rand = np.random.random()
                        #             prob = np.exp(-1/(Kb*T)*(new_Ham-old_Ham))
                        #             if rand > prob:
                        #                 _old_species = old_site_state[0]
                        #                 _old_phase = old_site_state[1]
                        #                 _old_spin = old_site_state[2]
                        #                 supercell_obj.set_site_species(site,_old_species)
                        #                 supercell_obj.set_site_phase(site,_old_phase)
                        #                 supercell_obj.set_site_spin(site,_old_spin)
                        #                 _old_species2 = old_site_2_state[0]
                        #                 _old_phase2 = old_site_2_state[1]
                        #                 _old_spin2 = old_site_2_state[2]
                        #                 supercell_obj.set_site_species(site_2,_old_species2)
                        #                 supercell_obj.set_site_phase(site_2,_old_phase2)
                        #                 supercell_obj.set_site_spin(site_2,_old_spin2)
                        #                 inc_not += 1
                        #             else:
                        #                 #H_total += new_Ham-old_Ham
                        #                 inc_up += 1
                        #         else:
                        #             #H_total += new_Ham-old_Ham
                        #             inc_down += 1

                        old_Ham = eval_site_new(site,supercell_obj,Cluster_rules,J_rules,Js,T)
                        old_spin = flip_spin(site,supercell_obj)
                        new_Ham = eval_site_new(site,supercell_obj,Cluster_rules,J_rules,Js,T)
                        if new_Ham < old_Ham:
                            H_total -= abs(new_Ham-old_Ham)
                            inc_down += 1
                            H_down += (new_Ham-old_Ham)
                        else:
                            rand = np.random.random()
                            prob = math.exp(-1/(Kb*T)*(new_Ham-old_Ham))
                            if rand < prob:
                                H_total += abs(new_Ham-old_Ham)
                                inc_up += 1
                                H_up += (new_Ham-old_Ham)
                            else:
                                supercell_obj.set_site_spin(site,old_spin)
                                inc_not += 1
                        H_total += eval_site_new(site,supercell_obj,Cluster_rules,J_rules,Js,T)
                        M += calc_avg_spin(site,supercell_obj)
                        #plt.figure(2)
                        #plt.plot(-1+(inc_T)/(supercell_obj.i_length*supercell_obj.j_length*supercell_obj.k_length),H_total,lw=2,marker='o',color='b')
                        p += supercell_obj.get_site_phase(site)/supercell_obj.num_sites
                        p2 += supercell_obj.get_site_phase(site)**2/supercell_obj.num_sites
                        mag += supercell_obj.get_site_spin(site)/supercell_obj.num_sites
                        mag2 += supercell_obj.get_site_spin(site)**2/supercell_obj.num_sites
                        inc_T += 1
        H_total,total_phase,total_phase2,total_spin,total_spin2 = eval_lattice_new(supercell_obj,Cluster_rules,J_rules,Js,T)
        #H_total2,p,p2,mag,mag2 = eval_lattice(supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
        plt.figure(2)
        #plt.plot(passes,H_total/supercell_obj.num_sites,lw=3,marker='o',color='b')
        plt.plot(T,H_total,lw=3,marker='o',color='b')
        #plt.plot(passes,H_total,lw=3,marker='o',color='r')
        plt.figure(3)
        plt.subplot(311)
        plt.xlabel("Sweeps", fontsize=10)
        plt.ylabel("Magnetic Order Parameter", fontsize=10)
        plt.plot(T,M/supercell_obj.num_sites,lw=3,marker='o',color='g')
        plt.subplot(312)
        plt.xlabel("Sweeps", fontsize=10)
        plt.ylabel("(Average Mag)^2", fontsize=10)
        plt.plot(T,mag**2,lw=3,marker='o',color='g')###########
        plt.figure(4)
        plt.subplot(411)
        plt.xlabel("Sweeps", fontsize=10)
        plt.ylabel("Average Phase", fontsize=10)
        plt.plot(T,p,lw=3,marker='o',color='r')
        plt.subplot(412)
        plt.xlabel("Sweeps", fontsize=10)
        plt.ylabel("Average Phase^2", fontsize=10)
        plt.plot(T,p2,lw=3,marker='o',color='r')
        print([H_total,H_total])
        # inc +=1
        # if inc >= 100:
        #     T -= 100
        #     inc = 0
        #     if T <= 0:
        #         T = 1
        T += 10
        if T <= 0:
           T = 1
        print(T)
        #inc += 1
        # if inc >= 100:
        #     inc = 0

        if passes >= numb_passes*.9:
            H_avg += H_total/(numb_passes*.1)
            mag_avg += mag/(numb_passes*.1)
            mag2_avg += mag2/(numb_passes*.1)
            p_avg += p/(numb_passes*.1)
            p2_avg += p2/(numb_passes*.1)


    plt.figure(2)
    plt.xlabel("Sweeps", fontsize=20)
    plt.ylabel("Energy of lattice (eV)", fontsize=20)
    plt.savefig('Enrg.png')
    plt.figure(3)
    plt.savefig('Mag.png')
    plt.figure(4)
    plt.savefig('Phase.png')
    temp_output = open('Temp_data','a')
    temp_output.write(str(T)+'  '+str(H_avg/supercell_obj.num_sites)+'  '+str(mag_avg)+'  '+str(mag2_avg)+'  '+str(p_avg)+'  '+str(p2_avg)+'\n')
    temp_output.close()
    print("\n")
    print(inc_down)
    print(inc_up)
    print(inc_not)
    print(inc_T)
    print("\n")
    print(H_down)
    print(H_up)

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

# Determine the total energy of the entire lattice and return that energy
def eval_lattice(supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T,do_figs=True):
    total_Ham = 0
    total_phase = 0
    total_phase2 = 0
    total_spin = 0
    total_spin2 = 0
    for i in range(supercell_obj.i_length):
        for j in range(supercell_obj.j_length):
            for k in range(supercell_obj.k_length):
                site = [i,j,k]
                total_Ham += eval_site(site,supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
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

def simple_cluster_flip_phase(site,neighbors,supercell_obj):
    old_phase = supercell_obj.get_site_phase(site)
    old_neighbor_phase  = []
    for i in range(len(neighbors)):
        old_neighbor_phase.append(supercell_obj.get_site_phase(neighbors[i]))
    phase_changed = False
    while phase_changed == False:
        rand = np.random.random()
        if rand <= 1/2.0:
            phase = 0
        elif rand > 1/2.0 and rand <= 3/4.0:
            phase = -1
        elif rand >3/4.0:
            phase = 1
        if phase != old_phase:
            phase_changed = True
    supercell_obj.set_site_phase(site,phase)
    for i in range(len(neighbors)):
        supercell_obj.set_site_phase(neighbors[i],phase)
    return old_phase, old_neighbor_phase

#-# Randomly change the species of a specific element in the lattice and return the value of the original species
def flip_species(site_1,site_2,supercell_obj):
    old_species_1 = supercell_obj.get_site_species(site_1)
    old_phase_1 = supercell_obj.get_site_phase(site_1)
    old_spin_1 = supercell_obj.get_site_spin(site_1)
    old_state_1 = [old_species_1,old_phase_1,old_spin_1]

    old_species_2 = supercell_obj.get_site_species(site_2)
    old_phase_2 = supercell_obj.get_site_phase(site_2)
    old_spin_2 = supercell_obj.get_site_spin(site_2)
    old_state_2 = [old_species_2,old_phase_2,old_spin_2]

    supercell_obj.set_site_species(site_1,old_species_2)
    supercell_obj.set_site_phase(site_1,old_phase_2)
    supercell_obj.set_site_spin(site_1,old_spin_2)
    supercell_obj.set_site_species(site_2,old_species_1)
    supercell_obj.set_site_phase(site_2,old_phase_1)
    supercell_obj.set_site_spin(site_2,old_spin_2)
    return old_state_1,old_state_2

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

def get_new_phase(site,supercell_obj):
    old_phase = supercell_obj.get_site_phase(site)
    phase_changed = False
    while phase_changed == False:
        rand = np.random.random()
        if rand <= 1/3.0:
            phase = -1
        elif rand > 1/3.0 and rand <= 2/3.0:
            phase = -1
        elif rand >2/3.0:
            phase = 0
        if phase != old_phase:
            phase_changed = True
    return phase

def calc_avg_spin(site,supercell_obj):
    M = 0
    count = 0
    for i in range(supercell_obj.get_number_of_neighbors(site)):
        neighbor = i
        if supercell_obj.get_neighbor_order(site,neighbor) == 2:
            M += supercell_obj.get_site_spin(site)*supercell_obj.get_neighbor_spin(site,neighbor)
            count += 1
    return M/count

def run_montecarlo(supercell_obj,numb_passes,temp,BEG_rules,Cluster_rules,J_rules,Js,do_figs=True):
    T = temp
    Kb = .000086173324 #8.6173324(78)×10−5 eV*K^-1
    H_avg = 0
    mag_avg = 0
    mag2_avg = 0
    p_avg = 0
    p2_avg = 0
    inc_down = 0
    inc_up = 0
    inc_not = 0
    inc_T = 0
    H_down = 0
    H_up = 0
    H = 0
    numb_sub_passes = 5
    H_total,p,p2,mag,mag2 = eval_lattice(supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
    plt.figure(2)
    plt.plot(-1,H_total,lw=3,marker='o',color='b')
    plt.plot(-1,H_total,lw=3,marker='o',color='r')
    print(H_total)
    for passes in range(numb_passes):
        M_plot = 0
        mag_plot = 0
        phase_plot = 0
        phase2_plot = 0
        enrg_plot = 0
        for sub_passes in range(numb_sub_passes):
            H_total = 0
            p = 0
            p2 = 0
            mag = 0
            mag2 = 0
            M = 0
            for i in range(supercell_obj.i_length):
                for j in range(supercell_obj.j_length):
                    for k in range(supercell_obj.k_length):
                        site = [i,j,k]
                        neighbor_found = False
                        neighbor = 0
                        old_Ham = 0
                        new_Ham = 0
                        old_Ham = eval_site(site,supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
                        neighbor_site = supercell_obj.get_neighbor_pos(site,neighbor)
                        old_Ham += eval_site(neighbor_site,supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
                        old_phase,old_neighbor_phase = flip_phase(site,neighbor,supercell_obj)
                        new_Ham = eval_site(site,supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
                        new_Ham += eval_site(neighbor_site,supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
                        if new_Ham < old_Ham:
                            inc_down += 1
                        else:
                            rand = np.random.random()
                            prob = np.exp(-1/(Kb*T)*(new_Ham-old_Ham))
                            if rand > prob:
                                supercell_obj.set_site_phase(site,old_phase)
                                supercell_obj.set_neighbor_phase(site,neighbor,old_neighbor_phase)
                                inc_not += 1
                            else:
                                #H_total += (new_Ham-old_Ham)
                                inc_up += 1

                        # if supercell_obj.get_site_species(site) != 0:
                        #     composition = supercell_obj.get_composition()
                        #     if composition[2] != 0:
                        #         site_found = False
                        #         while site_found == False:
                        #             site_2 = [np.random.randint(0,supercell_obj.i_length),np.random.randint(0,supercell_obj.j_length),np.random.randint(0,supercell_obj.k_length)]
                        #             if supercell_obj.get_site_species(site_2) != supercell_obj.get_site_species(site):
                        #                 if supercell_obj.get_site_species(site_2) != 0:
                        #                     site_found = True
                        #         old_Ham = eval_site(site,supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
                        #         old_Ham += eval_site(site_2,supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
                        #         old_site_state,old_site_2_state = flip_species(site,site_2,supercell_obj)
                        #         new_Ham = eval_site(site,supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
                        #         new_Ham += eval_site(site_2,supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
                        #         if new_Ham > old_Ham:
                        #             rand = np.random.random()
                        #             prob = np.exp(-1/(Kb*T)*(new_Ham-old_Ham))
                        #             if rand > prob:
                        #                 _old_species = old_site_state[0]
                        #                 _old_phase = old_site_state[1]
                        #                 _old_spin = old_site_state[2]
                        #                 supercell_obj.set_site_species(site,_old_species)
                        #                 supercell_obj.set_site_phase(site,_old_phase)
                        #                 supercell_obj.set_site_spin(site,_old_spin)
                        #                 _old_species2 = old_site_2_state[0]
                        #                 _old_phase2 = old_site_2_state[1]
                        #                 _old_spin2 = old_site_2_state[2]
                        #                 supercell_obj.set_site_species(site_2,_old_species2)
                        #                 supercell_obj.set_site_phase(site_2,_old_phase2)
                        #                 supercell_obj.set_site_spin(site_2,_old_spin2)
                        #                 inc_not += 1
                        #             else:
                        #                 #H_total += new_Ham-old_Ham
                        #                 inc_up += 1
                        #         else:
                        #             #H_total += new_Ham-old_Ham
                        #             inc_down += 1

                        old_Ham = eval_site(site,supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
                        old_spin = flip_spin(site,supercell_obj)
                        new_Ham = eval_site(site,supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
                        if new_Ham < old_Ham:
                            H_total -= abs(new_Ham-old_Ham)
                            inc_down += 1
                            H_down += (new_Ham-old_Ham)
                        else:
                            rand = np.random.random()
                            prob = math.exp(-1/(Kb*T)*(new_Ham-old_Ham))
                            if rand < prob:
                                H_total += abs(new_Ham-old_Ham)
                                inc_up += 1
                                H_up += (new_Ham-old_Ham)
                            else:
                                supercell_obj.set_site_spin(site,old_spin)
                                inc_not += 1
                        H_total += eval_site(site,supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
                        M += calc_avg_spin(site,supercell_obj)
                        #plt.figure(2)
                        #plt.plot(-1+(inc_T)/(supercell_obj.i_length*supercell_obj.j_length*supercell_obj.k_length),H_total,lw=2,marker='o',color='b')
                        p += supercell_obj.get_site_phase(site)/supercell_obj.num_sites
                        p2 += supercell_obj.get_site_phase(site)**2/supercell_obj.num_sites
                        mag += supercell_obj.get_site_spin(site)/supercell_obj.num_sites
                        mag2 += supercell_obj.get_site_spin(site)**2/supercell_obj.num_sites
                        inc_T += 1
        H_total,total_phase,total_phase2,total_spin,total_spin2 = eval_lattice(supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
        #H_total2,p,p2,mag,mag2 = eval_lattice(supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
        plt.figure(2)
        #plt.plot(passes,H_total/supercell_obj.num_sites,lw=3,marker='o',color='b')
        plt.plot(T,H_total,lw=3,marker='o',color='b')
        #plt.plot(passes,H_total,lw=3,marker='o',color='r')
        plt.figure(3)
        plt.subplot(311)
        plt.xlabel("Sweeps", fontsize=10)
        plt.ylabel("Magnetic Order Parameter", fontsize=10)
        plt.plot(T,M/supercell_obj.num_sites,lw=3,marker='o',color='g')
        plt.subplot(312)
        plt.xlabel("Sweeps", fontsize=10)
        plt.ylabel("(Average Mag)^2", fontsize=10)
        plt.plot(T,mag**2,lw=3,marker='o',color='g')###########
        plt.figure(4)
        plt.subplot(411)
        plt.xlabel("Sweeps", fontsize=10)
        plt.ylabel("Average Phase", fontsize=10)
        plt.plot(T,p,lw=3,marker='o',color='r')
        plt.subplot(412)
        plt.xlabel("Sweeps", fontsize=10)
        plt.ylabel("Average Phase^2", fontsize=10)
        plt.plot(T,p2,lw=3,marker='o',color='r')
        print([H_total,H_total])
        # inc +=1
        # if inc >= 100:
        #     T -= 100
        #     inc = 0
        #     if T <= 0:
        #         T = 1
        T += 10
        if T <= 0:
           T = 1
        print(T)
        #inc += 1
        # if inc >= 100:
        #     inc = 0

        if passes >= numb_passes*.9:
            H_avg += H_total/(numb_passes*.1)
            mag_avg += mag/(numb_passes*.1)
            mag2_avg += mag2/(numb_passes*.1)
            p_avg += p/(numb_passes*.1)
            p2_avg += p2/(numb_passes*.1)

    if do_figs is True:
        plt.figure(2)
        plt.xlabel("Sweeps", fontsize=20)
        plt.ylabel("Energy of lattice (eV)", fontsize=20)
        plt.savefig('Enrg.png')
        plt.figure(3)
        plt.savefig('Mag.png')
        plt.figure(4)
        plt.savefig('Phase.png')
    temp_output = open('Temp_data','a')
    temp_output.write(str(T)+'  '+str(H_avg/supercell_obj.num_sites)+'  '+str(mag_avg)+'  '+str(mag2_avg)+'  '+str(p_avg)+'  '+str(p2_avg)+'\n')
    temp_output.close()
    print("\n")
    print(inc_down)
    print(inc_up)
    print(inc_not)
    print(inc_T)
    print("\n")
    print(H_down)
    print(H_up)

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
                ws.append(supercell_obj.get_site_phase(site))
                if supercell_obj.get_site_species(site) == 0:
                    cs.append('g')
                if supercell_obj.get_site_species(site) == 1:
                    cs.append('r')
                if supercell_obj.get_site_species(site) == 2:
                    cs.append('b')
    ax.quiver(xs,ys,zs,us,vs,ws,pivot='middle',length=.5)
    ax.scatter(xs,ys,zs,c=cs,marker='o',s=50)
    plt.savefig('3D_plt.png')

#-# Runs the Wolf/Mixed Cluster Algorithm
def run_WA_MCA(supercell_obj,numb_passes,num_sub_passes,temp,temp_inc,Cluster_rules,J_rules,Js,do_figs=True):
    T = temp
    Kb = .000086173324
    B = 1/(Kb*T)
    inc_down = 0
    inc_up = 0
    inc_not = 0
    M = 0
    H_total,total_phase,total_phase2,total_spin,total_spin2 = eval_lattice_new(supercell_obj,Cluster_rules,J_rules,Js,T)
    plt.figure(2)
    plt.plot(-1,H_total,lw=3,marker='o',color='b')
    plt.plot(-1,H_total,lw=3,marker='o',color='r')
    for passes in range(numb_passes):
        #Flip spins and Species
        print(T)
        for sub_passes in range(num_sub_passes):
            M = 0
            for i in range(supercell_obj.i_length):
                for j in range(supercell_obj.j_length):
                    for k in range(supercell_obj.k_length):
                        site = [i,j,k]
                        old_Ham = eval_site_new(site,supercell_obj,Cluster_rules,J_rules,Js,T)
                        old_spin = flip_spin(site,supercell_obj)
                        new_Ham = eval_site_new(site,supercell_obj,Cluster_rules,J_rules,Js,T)
                        #print([new_Ham,old_Ham])
                        if new_Ham < old_Ham:
                            inc_down += 1
                            #print('accept down')
                        else:
                            rand = np.random.random()
                            prob = math.exp(-1*B*(new_Ham-old_Ham))
                            #print(prob)
                            if rand < prob:
                                inc_up += 1
                                #print('accept up')
                            else:
                                supercell_obj.set_site_spin(site,old_spin)
                                inc_not += 1
                                #print('not accepted')
                        M += calc_avg_spin(site,supercell_obj)
                        ##############
                        # FLIP SPECIES
                        ##############
        # #Randdom Seed
        # cluster = []
        # seed =(np.random.randint(0,supercell_obj.i_length),np.random.randint(0,supercell_obj.j_length),np.random.randint(0,supercell_obj.k_length))
        # seed_phase = supercell_obj.get_site_phase(seed)
        # new_phase = get_new_phase(seed,supercell_obj)
        # grow_cluster(seed,supercell_obj,seed_phase,new_phase,cluster,Cluster_rules,J_rules,Js,T)
        # ### Track size here (print(len(cluster))
        # print([seed_phase,new_phase])
        # print(len(cluster))
        # if seed_phase*new_phase == -1:
        #     print('W')
        #     flip_cluster(supercell_obj,seed_phase,new_phase,cluster)
        # else:
        #     print('MC')
        #     H_cluster_old = eval_cluster(supercell_obj,seed_phase,new_phase,cluster,Cluster_rules,J_rules,Js,T)
        #     flip_cluster(supercell_obj,seed_phase,new_phase,cluster)
        #     H_cluster_new = eval_cluster(supercell_obj,seed_phase,new_phase,cluster,Cluster_rules,J_rules,Js,T)
        #     if H_cluster_new <= H_cluster_old:
        #         inc_down += 1
        #         print('yup')
        #     else:
        #         print(H_cluster_new-H_cluster_old)
        #         rand = np.random.random()
        #         prob = math.exp(-B*(H_cluster_new-H_cluster_old))
        #         print(prob)
        #         print([H_cluster_new,H_cluster_old])
        #         if rand < prob:
        #             print('yup')
        #             inc_up += 1
        #         else:
        #             flip_cluster(supercell_obj,seed_phase,new_phase,cluster)
        #             inc_not += 1

        X_axis = T
        H_total,total_phase,total_phase2,total_spin,total_spin2 = eval_lattice_new(supercell_obj,Cluster_rules,J_rules,Js,T)
        plt.figure(2)
        plt.plot(X_axis,H_total,lw=3,marker='o',color='b')
        plt.figure(3)
        plt.subplot(311)
        plt.xlabel("Sweeps", fontsize=10)
        plt.ylabel("Magnetic Order Parameter", fontsize=10)
        plt.plot(X_axis,M/supercell_obj.num_sites,lw=3,marker='o',color='g')
        plt.subplot(312)
        plt.xlabel("Sweeps", fontsize=10)
        plt.ylabel("(Average Mag)^2", fontsize=10)
        plt.plot(X_axis,total_spin2,lw=3,marker='o',color='g')###########
        plt.figure(4)
        plt.subplot(411)
        plt.xlabel("Sweeps", fontsize=10)
        plt.ylabel("Average Phase", fontsize=10)
        plt.plot(X_axis,total_phase,lw=3,marker='o',color='r')
        plt.subplot(412)
        plt.xlabel("Sweeps", fontsize=10)
        plt.ylabel("Average Phase^2", fontsize=10)
        plt.plot(X_axis,total_phase2,lw=3,marker='o',color='r')
        print([H_total,H_total])

        T += temp_inc
        print(T)
        print('\n')

    if do_figs is True:
        plt.figure(2)
        plt.xlabel("Sweeps", fontsize=20)
        plt.ylabel("Energy of lattice (eV)", fontsize=20)
        plt.savefig('Enrg.png')
        plt.figure(3)
        plt.savefig('Mag.png')
        plt.figure(4)
        plt.savefig('Phase.png')
    #temp_output = open('Temp_data','a')
    #temp_output.write(str(T)+'  '+str(H_avg/supercell_obj.num_sites)+'  '+str(mag_avg)+'  '+str(mag2_avg)+'  '+str(p_avg)+'  '+str(p2_avg)+'\n')
    #temp_output.close()
    print("\n")
    print(inc_down)
    print(inc_up)
    print(inc_not)
    print("\n")

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
def grow_cluster(site,supercell_obj,seed_phase,new_phase,links,Cluster_rules,J_rules,Js,T): # Recursive function
    Kb = .000086173324
    B = 1/(Kb*T)
    site_phase = supercell_obj.get_site_phase(site)
    J,K = calc_BEG_params(site,supercell_obj,Cluster_rules,J_rules,Js,T)
    BEG_J = B*J
    BEG_K = B*K
    links.append(site)
    # Wolf Algorithm
    if new_phase*seed_phase == -1:
        for neighbor in range(supercell_obj.get_number_of_neighbors(site)):
            if supercell_obj.get_neighbor_order(site,neighbor) == 1:
                if supercell_obj.get_neighbor_phase(site,neighbor) == seed_phase:
                    if supercell_obj.get_neighbor_pos(site,neighbor) not in links:
                        rand = np.random.random()
                        prob = 1-np.exp(2*BEG_J)
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
                            prob = 1-np.exp(BEG_J+BEG_K/3)
                            if rand < prob:
                                new_site = supercell_obj.get_neighbor_pos(site,neighbor)
                                grow_cluster(new_site,supercell_obj,seed_phase,new_phase,links,Cluster_rules,J_rules,Js,T)
                        else:
                            rand = np.random.random()
                            prob = 1-np.exp(BEG_J-BEG_K/3)
                            if rand < prob:
                                links.append(site)
                                new_site = supercell_obj.get_neighbor_pos(site,neighbor)
                                grow_cluster(new_site,supercell_obj,seed_phase,new_phase,links,Cluster_rules,J_rules,Js,T)

    if [seed_phase,new_phase] == [-1,0] or [seed_phase,new_phase] == [0,1]:
        for neighbor in range(supercell_obj.get_number_of_neighbors(site)):
            if supercell_obj.get_neighbor_order(site,neighbor) == 1:
                if supercell_obj.get_neighbor_pos(site,neighbor) not in links:
                        if supercell_obj.get_neighbor_phase(site,neighbor) == -1 or supercell_obj.get_neighbor_phase(site,neighbor) == 0:
                            if supercell_obj.get_neighbor_phase(site,neighbor) == site_phase:
                                rand = np.random.random()
                                prob = 1-np.exp(BEG_J+BEG_K/3)
                                if rand < prob:
                                    new_site = supercell_obj.get_neighbor_pos(site,neighbor)
                                    grow_cluster(new_site,supercell_obj,seed_phase,new_phase,links,Cluster_rules,J_rules,Js,T)
                            else:
                                rand = np.random.random()
                                prob = 1-np.exp(BEG_J-BEG_K/3)
                                if rand < prob:
                                    new_site = supercell_obj.get_neighbor_pos(site,neighbor)
                                    grow_cluster(new_site,supercell_obj,seed_phase,new_phase,links,Cluster_rules,J_rules,Js,T)

#-# Evaluates the total energy of the cluster
def eval_cluster(supercell_obj,seed_phase,new_phase,links,Cluster_rules,J_ruels,Js,T):
    Kb = .000086173324
    total_H = 0
    for i in range(len(links)):
        site = links[i]
        #print('site number is ',i,' and this site has ',supercell_obj.get_number_of_neighbors(site),' number of neighbors.')
        site_phase = supercell_obj.get_site_phase(site)
######### START ELIF COMMENT HERE #############
# This line here in particular, I am concerned that we will evaluate the energy including all neighbors of all atoms in the cluster,
# whether the neighbors are in the cluster or not.  But maybe that is ok.  This is just the site specific J,K calculation
# which I guess will depend on all neighbors.
        BEG_J,BEG_K = calc_BEG_params(site,supercell_obj,Cluster_rules,J_ruels,Js,T)
######### END ELIF COMMENT HERE #############
        for neighbor in range(supercell_obj.get_number_of_neighbors(site)):
            if supercell_obj.get_neighbor_order(site,neighbor) == 1:
                if supercell_obj.get_neighbor_pos(site,neighbor) in links:
                    neighbor_phase = supercell_obj.get_neighbor_phase(site,neighbor)
                    total_H += BEG_J*site_phase*neighbor_phase+BEG_K*(1-site_phase**2)*(1-neighbor_phase**2)
######### START ELIF MODIFICATIONS #############
# Notice that I am not including the constant offset here again.
# But again, if I were it would be divided by two, -Kb*T*np.log(2)/2
        #total_H += -Kb*T*np.log(2)*(1-site_phase**2)
        total_H += Kb*T*np.log(2)*(site_phase**2)
######## END ELIF MODIFICATIONS  ##############################
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

def run_simple_cluster_MC(supercell_obj,numb_passes,temp,BEG_rules,Cluster_rules,J_rules,Js,do_figs=True):
    T = temp
    Kb = .000086173324 #8.6173324(78)×10−5 eV*K^-1
    H_avg = 0
    mag_avg = 0
    mag2_avg = 0
    p_avg = 0
    p2_avg = 0
    inc = 0
    inc_not = 0
    H_total,p,p2,mag,mag2 = eval_lattice(supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
    for passes in range(numb_passes):
        for i in range(supercell_obj.i_length):
            for j in range(supercell_obj.j_length):
                for k in range(supercell_obj.k_length):
                    site = [i,j,k]
                    rand_move_type = np.random.random()
                    if rand_move_type <=.5:
                        #print('single')
                        old_Ham = eval_site(site,supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
                        old_phase = flip_phase(site,supercell_obj)
                        new_Ham = eval_site(site,supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
                        #print([old_Ham,new_Ham])
                        if new_Ham > old_Ham:
                            rand = np.random.random()
                            prob = np.exp(-1/(Kb*T)*(new_Ham-old_Ham))
                            if rand > prob:
                                supercell_obj.set_site_phase(site,old_phase)
                                inc_not += 1
                            else:
                                H_total += new_Ham-old_Ham
                                inc += 1
                        else:
                            H_total += new_Ham-old_Ham
                            inc += 1
                    else:
                        #print('picked')
                        neighbors = []
                        for l in range(supercell_obj.get_number_of_neighbors(site)):
                            if supercell_obj.get_neighbor_order(site,l) == BEG_rules[0].neighbor_order: ######################### fix up this in future
                                neighbors.append(supercell_obj.get_neighbor_pos(site,l))
                        old_Ham = eval_site(site,supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
                        for l in range(len(neighbors)):
                            old_Ham += eval_site(neighbors[l],supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
                        old_phase,old_neighbor_phase = simple_cluster_flip_phase(site,neighbors,supercell_obj)
                        new_Ham = eval_site(site,supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
                        for l in range(len(neighbors)):
                            new_Ham += eval_site(neighbors[l],supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
                        if new_Ham > old_Ham:
                            rand = np.random.random()
                            prob = np.exp(-1/(Kb*T)*(new_Ham-old_Ham))
                            #print(prob)
                            if rand > prob:
                                supercell_obj.set_site_phase(site,old_phase)
                                for m in range(len(neighbors)):
                                    supercell_obj.set_site_phase(neighbors[m],old_neighbor_phase[m])
                                inc_not += 1
                            else:
                                H_total += new_Ham-old_Ham
                                inc += 1
                        else:
                            H_total += new_Ham-old_Ham
                            inc += 1

        H_total2,p,p2,mag,mag2 = eval_lattice(supercell_obj,BEG_rules,Cluster_rules,J_rules,Js,T)
        print([H_total,H_total2])

        if passes >= numb_passes*.9:
            H_avg += H_total/(numb_passes*.1)
            mag_avg += mag/(numb_passes*.1)
            mag2_avg += mag2/(numb_passes*.1)
            p_avg += p/(numb_passes*.1)
            p2_avg += p2/(numb_passes*.1)

        plt.figure(2)
        #plt.plot(passes,H_total/supercell_obj.num_sites,lw=3,marker='o',color='b')
        plt.plot(passes,H_total,lw=3,marker='o',color='b')
        plt.figure(3)
        plt.subplot(311)
        plt.plot(passes,mag,lw=3,marker='o',color='g')
        plt.subplot(312)
        plt.plot(passes,mag2,lw=3,marker='o',color='g')
        plt.figure(4)
        plt.subplot(411)
        plt.plot(passes,p,lw=3,marker='o',color='r')
        plt.subplot(412)
        plt.plot(passes,p2,lw=3,marker='o',color='r')
    if do_figs is True:
        plt.figure(2)
        plt.savefig('Enrg.png')
        plt.figure(3)
        plt.savefig('Mag.png')
        plt.figure(4)
        plt.savefig('Phase.png')
    temp_output = open('Temp_data','a')
    temp_output.write(str(T)+'  '+str(H_avg/supercell_obj.num_sites)+'  '+str(mag_avg)+'  '+str(mag2_avg)+'  '+str(p_avg)+'  '+str(p2_avg)+'\n')
    temp_output.close()
    print(inc)
    print(inc_not)

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