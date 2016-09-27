__author__ = 'brian'
import numpy as np
from copy import deepcopy


def check_phase(home_site):
    home_phase = home_site.phase
    if abs(home_phase) == 1:
        phase = 'mart'
    if home_phase == 0:
        phase = 'aust'
    return phase

def eval_supercell(supercell_obj,beg_rule_list, cluster_rule_list,j_rule_list,js):
    h = 0
    mag = 0
    mag2 = 0
    p = 0
    p2 = 0
    supercell = supercell_obj.supercell
    for i in range(supercell_obj.i_length):
        for j in range(supercell_obj.j_length):
            for k in range(supercell_obj.k_length):
                site = supercell[i,j,k]
                h_site = eval_site(supercell,(i,j,k),beg_rule_list,cluster_rule_list,j_rule_list,js)
                h += float(h_site)
                mag += site.get_spin()/supercell_obj.num_sites
                mag2 += site.get_spin()**2/supercell_obj.num_sites
                p += site.get_phase()/supercell_obj.num_sites
                p2 += site.get_phase()**2/supercell_obj.num_sites

    return float(h),float(mag),float(mag2),float(p),float(p2)

def eval_site(suercell,index,beg_rule_list, cluster_rule_list,j_rule_list,js):
    h = 0
    i_index = index[0]
    j_index = index[1]
    k_index = index[2]
    home_site = suercell[i_index,j_index,k_index]
    inc_BEG = 0
    inc_Clust = 0
    inc_spin = 0
    beg_sum = 0
    clust_sum = 0
    spin_sum = 0
    for i in range(len(home_site.neighbors)):
        neighbor_data = home_site.neighbors[i]
        neighbor_site = suercell[neighbor_data.i_pos,neighbor_data.j_pos,neighbor_data.k_pos]
        # BEG Ham
        for j in range(len(beg_rule_list)):
            if neighbor_data.order == beg_rule_list[j].neighbor_order:
                if home_site.species in beg_rule_list[j].home_atom_list:
                    if neighbor_site.species in beg_rule_list[j].neighbor_atom_list:
                        if check_phase(home_site) == beg_rule_list[j].phase:
                            if beg_rule_list[j].phase == 'mart':
                                h += js[j] * home_site.phase * neighbor_site.phase/6
                                inc_BEG += js[j]/6
                                beg_sum += js[j] * home_site.phase * neighbor_site.phase/6
                            if beg_rule_list[j].phase == 'aust':
                                h += js[j] * (1 - home_site.phase**2) * (1-neighbor_site.phase**2)/6
                                inc_BEG += js[j]/6
                                beg_sum += js[j] * (1 - home_site.phase**2) * (1-neighbor_site.phase**2)/6
        # Cluster Ham
        for j in range(len(cluster_rule_list)):
            if home_site.species in cluster_rule_list[j].home_atom_list:
                if check_phase(home_site) == cluster_rule_list[j].phase:
                    if neighbor_data.plain == cluster_rule_list[j].plane or cluster_rule_list[j].plane == 'ALL':
                        if neighbor_data.order == cluster_rule_list[j].neighbor_order:
                            if cluster_rule_list[j].neighbor_arrangement == 'PERM':
                                if neighbor_site.species != home_site.species:
                                    if neighbor_site.species in cluster_rule_list[j].neighbor_atom_list:
                                        h += js[len(beg_rule_list) + j]
                                        inc_Clust += js[len(beg_rule_list) + j]
                                        clust_sum += js[len(beg_rule_list) + j]
                            if cluster_rule_list[j].neighbor_arrangement == 'COMB':
                                if neighbor_site.species in cluster_rule_list[j].neighbor_atom_list:
                                    h += js[len(beg_rule_list) + j]
                                    inc_Clust += js[len(beg_rule_list) + j]
                                    clust_sum += js[len(beg_rule_list) + j]
        # Mag Ham
        for j in range(len(j_rule_list)):
            if home_site.species in j_rule_list[j].home_atom_list:
                if check_phase(home_site) == j_rule_list[j].phase:
                    if neighbor_data.plain == j_rule_list[j].plane or j_rule_list[j].plane == 'ALL':
                        if neighbor_data.order == j_rule_list[j].neighbor_order:
                            if j_rule_list[j].neighbor_arrangement == 'PERM':
                                if neighbor_site.species != home_site.species:
                                    if neighbor_site.species in j_rule_list[j].neighbor_atom_list:
                                        h += home_site.spin*neighbor_site.spin* float(js[len(beg_rule_list) + len(cluster_rule_list)+j])
                                        inc_spin += js[len(beg_rule_list) + len(cluster_rule_list)+j]
                                        spin_sum += home_site.spin*neighbor_site.spin* float(js[len(beg_rule_list) + len(cluster_rule_list)+j])
                            if j_rule_list[j].neighbor_arrangement == 'COMB':
                                if neighbor_site.species in j_rule_list[j].neighbor_atom_list:
                                    h += home_site.spin*neighbor_site.spin * float(js[len(beg_rule_list) + len(cluster_rule_list)+j])
                                    inc_spin += js[len(beg_rule_list) + len(cluster_rule_list)+j]
                                    spin_sum += home_site.spin*neighbor_site.spin * float(js[len(beg_rule_list) + len(cluster_rule_list)+j])
    return float(beg_sum+clust_sum+spin_sum)

def flip_species(supercell_obj,pos):
    i_pos = pos[0]
    j_pos = pos[1]
    k_pos = pos[2]
    home_site = supercell_obj.supercell[i_pos,j_pos,k_pos]
    old_home_site = deepcopy(supercell_obj.supercell[i_pos,j_pos,k_pos])
    home_species = home_site.species
    site_flipped = False
    while site_flipped == False:
        rand_index = [np.random.randint(0,supercell_obj.i_length),np.random.randint(0,supercell_obj.j_length),np.random.randint(0,supercell_obj.k_length)]
        neighbor_site = supercell_obj.supercell[rand_index[0],rand_index[1],rand_index[2]]
        old_neighbor_site = deepcopy(supercell_obj.supercell[rand_index[0],rand_index[1],rand_index[2]])
        if neighbor_site.species != 0 and neighbor_site.species != home_site.species:
            home_site.species = neighbor_site.species
            neighbor_site.species = home_species
            site_flipped = True
    supercell_obj.supercell[i_pos,j_pos,k_pos] = home_site
    supercell_obj.supercell[rand_index[0],rand_index[1],rand_index[2]] = neighbor_site
    return old_home_site,old_neighbor_site,rand_index

def flip_phase(supercell,pos):
    i_pos = pos[0]
    j_pos = pos[1]
    k_pos = pos[2]
    home_site = supercell[i_pos,j_pos,k_pos]
    old_home_site = deepcopy(supercell[i_pos,j_pos,k_pos])
    old_phase = home_site.get_phase()
    phase_flipped = False
    while phase_flipped == False:
        rand = np.random.random()
        if rand <= 1/2.0:
            phase = 0
        elif rand > 1/2.0 and rand <= 3/4.0:
            phase = -1
        elif rand >3/4.0:
            phase = 1
        if phase != old_phase:
            phase_flipped = True
    home_site.set_phase(phase)
    supercell[i_pos,j_pos,k_pos] = home_site
    return old_home_site

def flip_spin(supercell,pos):
    i_pos = pos[0]
    j_pos = pos[1]
    k_pos = pos[2]
    home_site = supercell[i_pos,j_pos,k_pos]
    old_home_site = deepcopy(supercell[i_pos,j_pos,k_pos])
    if home_site.species == 2:
        home_site.set_spin(0)
    else:
        same_spin = False
        while same_spin == False:
            rand = np.random.random()
            if rand <= 1/3.0:
                spin = -1
            elif rand > 1/3.0 and rand < 2/3.0:
                spin = 0
            elif rand >= 2/3.0:
                spin = 1
            if spin != home_site.get_spin():
                same_spin = True
        home_site.set_spin(spin)
    supercell[i_pos,j_pos,k_pos] = home_site
    return old_home_site
