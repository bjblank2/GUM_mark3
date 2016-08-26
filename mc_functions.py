__author__ = 'brian'
import numpy as np
from copy import deepcopy

def init_supercell(num_pts):
    supercell = np.empty((num_pts,num_pts,num_pts),dtype=list)
    supercell_list = []
    count = 0
    for i in range(num_pts):
        for j in range(num_pts):
            for k in range(num_pts):
                if np.mod(k,2) == 0 :
                    params = [0,i*.5,j*.5,k*.25,1,0,0]
                else:
                    params = [0,i*.5+.25,j*.5+.25,k*.25,0,0,0]

                rand = np.random.random()
                if rand <= 1/3.0:
                    params[5] = -1 #########
                elif rand > 1/3.0 and rand <= 2/3.0:
                    params[5] = 0 #########
                elif rand >2/3.0:
                    params[5] = 1

                # if np.mod(i + j,2) == 0:
                #     params[5] = 1
                # else:
                #     params[5] = -1
                # if params[4] == 0:
                #     params[5] = 0

                rand = np.random.random()
                if rand <= 1/3.0:
                    params[6] = -1 #########
                elif rand > 1/3.0 and rand <= 2/3.0:
                    params[6] = 0 #########
                elif rand >2/3.0:
                    params[6] = 1
                params[0] = count
                count += 1
                supercell[i,j,k] = params
                supercell_list.append(params)
    return supercell_list,supercell

def apply_bc(i,inc,limit):
    if i + inc >= limit:
        new_i = i+inc-limit
    elif i + inc < 0:
        new_i = i + inc + limit
    else:
        new_i = i+inc
    return new_i

def calc_neighbors(supercell):
    length = len(supercell)
    neighbor_list = []
    neighbor_plain_list = []
    for i in range(length):
        for j in range(length):
            for k in range(length):
                nn_list = ['_']*8
                nn_plain = ['-']*8
                nnn_list = ['_']*6
                nnn_plain = ['_']*6
                nnnn_list = ['_']*12
                nnnn_plain = ['_']*12

                nn = supercell[i,j,apply_bc(k,1,length)]
                nn_list[0] = nn[0]
                nn_plain[0] = 'OUT'
                nn = supercell[i,j,apply_bc(k,-1,length)]
                nn_list[1] = nn[0]
                nn_plain[1] = 'OUT'
                nn = supercell[i,apply_bc(j,-1,length),apply_bc(k,1,length)]
                nn_list[2] = nn[0]
                nn_plain[2] = 'OUT'
                nn = supercell[i,apply_bc(j,-1,length),apply_bc(k,-1,length)]
                nn_list[3] = nn[0]
                nn_plain[3] = 'OUT'
                nn = supercell[apply_bc(i,-1,length),apply_bc(j,-1,length),apply_bc(k,1,length)]
                nn_list[4] = nn[0]
                nn_plain[4] = 'OUT'
                nn = supercell[apply_bc(i,-1,length),apply_bc(j,-1,length),apply_bc(k,-1,length)]
                nn_list[5] = nn[0]
                nn_plain[5] = 'OUT'
                nn = supercell[apply_bc(i,-1,length),j,apply_bc(k,1,length)]
                nn_list[6] = nn[0]
                nn_plain[6] = 'OUT'
                nn = supercell[apply_bc(i,-1,length),j,apply_bc(k,-1,length)]
                nn_list[7] = nn[0]
                nn_plain[7] = 'OUT'

                nnn = supercell[apply_bc(i,1,length),j,k]
                nnn_list[0] = nnn[0]
                nnn_plain[0] = 'IN'
                nnn = supercell[apply_bc(i,-1,length),j,k]
                nnn_list[1] = nnn[0]
                nnn_plain[1] = 'IN'
                nnn = supercell[i,apply_bc(j,1,length),k]
                nnn_list[2] = nnn[0]
                nnn_plain[2] = 'IN'
                nnn = supercell[i,apply_bc(j,-1,length),k]
                nnn_list[3] = nnn[0]
                nnn_plain[3] = 'IN'
                nnn = supercell[i,j,apply_bc(k,2,length)]
                nnn_list[4] = nnn[0]
                nnn_plain[4] = 'OUT'
                nnn = supercell[i,j,apply_bc(k,-2,length)]
                nnn_list[5] = nnn[0]
                nnn_plain[5] = 'OUT'

                nnnn = supercell[apply_bc(i,1,length),j,apply_bc(k,2,length)]
                nnnn_list[0] = nnnn[0]
                nnnn_plain[0] = 'OUT'
                nnnn = supercell[apply_bc(i,1,length),j,apply_bc(k,-2,length)]
                nnnn_list[1] = nnnn[0]
                nnnn_plain[1] = 'OUT'
                nnnn = supercell[apply_bc(i,-1,length),j,apply_bc(k,2,length)]
                nnnn_list[2] = nnnn[0]
                nnnn_plain[2] = 'OUT'
                nnnn = supercell[apply_bc(i,-1,length),j,apply_bc(k,-2,length)]
                nnnn_list[3] = nnnn[0]
                nnnn_plain[3] = 'OUT'
                nnnn = supercell[i,apply_bc(j,1,length),apply_bc(k,2,length)]
                nnnn_list[4] = nnnn[0]
                nnnn_plain[4] = 'OUT'
                nnnn = supercell[i,apply_bc(j,1,length),apply_bc(k,-2,length)]
                nnnn_list[5] = nnnn[0]
                nnnn_plain[5] = 'OUT'
                nnnn = supercell[i,apply_bc(j,-1,length),apply_bc(k,2,length)]
                nnnn_list[6] = nnnn[0]
                nnnn_plain[6] = 'OUT'
                nnnn = supercell[i,apply_bc(j,-1,length),apply_bc(k,-2,length)]
                nnnn_list[7] = nnnn[0]
                nnnn_plain[7] = 'OUT'
                nnnn = supercell[apply_bc(i,1,length),apply_bc(j,1,length),k]
                nnnn_list[8] = nnnn[0]
                nnnn_plain[8] = 'IN'
                nnnn = supercell[apply_bc(i,1,length),apply_bc(j,-1,length),k]
                nnnn_list[9] = nnnn[0]
                nnnn_plain[9] = 'IN'
                nnnn = supercell[apply_bc(i,-1,length),apply_bc(j,1,length),k]
                nnnn_list[10] = nnnn[0]
                nnnn_plain[10] = 'IN'
                nnnn = supercell[apply_bc(i,-1,length),apply_bc(j,-1,length),k]
                nnnn_list[11] = nnnn[0]
                nnnn_plain[11] = 'IN'

                neighbor_list.append([nn_list,nnn_list,nnnn_list])
                neighbor_plain_list.append([nn_plain,nnn_plain,nnnn_plain])
    return neighbor_list,neighbor_plain_list

def check_plane(home_site, plain):
        if home_site[6] == 0:
            plane_check = 'ALL'
        else:
            plane_check = plain
        return plane_check

def check_phase(home_site):
    home_phase = home_site[6]
    if abs(home_phase) == 1:
        phase = 'mart'
    if home_phase == 0:
        phase = 'aust'
    return phase

def eval_supercell(supercell_list,neighbors,neighbor_plain,beg_rule_list, cluster_rule_list,j_rule_list,js):
    h = 0
    h_mag = 0
    for i in range(len(supercell_list)):
        home_site = supercell_list[i]
        neighbor_lists = neighbors[i]
        neighbor_plain_lists = neighbor_plain[i]
        h_site,h_mag_site = eval_site(home_site,supercell_list,neighbor_lists,neighbor_plain_lists,beg_rule_list,cluster_rule_list,j_rule_list,js)
        h += float(h_site)
        h_mag += float(h_mag_site)
    return float(h),float(h_mag)

def eval_site(home_site,supercell_list,neighbors,plains,beg_rule_list, cluster_rule_list,j_rule_list,js):
    h = 0
    h_mag = 0
    beg_h = 0
    beg_sum_J = 0
    beg_sum_K = 0

    for i in range(len(neighbors)):
        neighbor_list = neighbors[i]
        plain_list = plains[i]
        for j in range(len(neighbor_list)):
            neighbor_index = neighbor_list[j]
            plain = plain_list[j]
            neighbor_site = supercell_list[neighbor_index]
            # BEG Ham
            for k in range(len(beg_rule_list)):
                if i+1 == beg_rule_list[k].neighbor_order:
                    if home_site[4] in beg_rule_list[k].home_atom_list:
                        if neighbor_site[4] in beg_rule_list[k].neighbor_atom_list:
                            if k == 0:
                                beg_sum_J += float(home_site[6] * neighbor_site[6])
                                beg_h += js[k] * home_site[6] * neighbor_site[6]
                            if k == 1:
                                beg_sum_K += float((1 - float(home_site[6])**2) * (1-float(neighbor_site[6])**2))
                                beg_h += js[k] * (1 - home_site[6]**2) * (1-neighbor_site[6]**2)
            # Cluster Ham
            for k in range(len(cluster_rule_list)):
                if home_site[4] in cluster_rule_list[k].home_atom_list:
                    if check_phase(home_site) == cluster_rule_list[k].phase:
                        if check_plane(home_site,plain) == cluster_rule_list[k].plane or check_plane(home_site,plain) == 'ALL':
                            if i+1 == cluster_rule_list[k].neighbor_order:
                                if cluster_rule_list[k].neighbor_arrangement == 'PERM':
                                    if neighbor_site[4] != home_site[4]:
                                        if neighbor_site[4] in cluster_rule_list[k].neighbor_atom_list:
                                            h += js[len(beg_rule_list) + k]
                                if cluster_rule_list[k].neighbor_arrangement == 'COMB':
                                    if neighbor_site[4] in cluster_rule_list[k].neighbor_atom_list:
                                        h += js[len(beg_rule_list) + k]
            # Mag Ham
            for k in range(len(j_rule_list)):
                if home_site[4] in j_rule_list[k].home_atom_list:
                    if check_phase(home_site) == j_rule_list[k].phase:
                        if check_plane(home_site,plain) == j_rule_list[k].plane or check_plane(home_site,plain) == 'ALL':
                            if i+1 == j_rule_list[k].neighbor_order:
                                if j_rule_list[k].neighbor_arrangement == 'PERM':
                                    if neighbor_site[4] != home_site[4]:
                                        if neighbor_site[4] in j_rule_list[k].neighbor_atom_list:
                                            h += home_site[5]*neighbor_site[5]* float(js[len(beg_rule_list) + len(cluster_rule_list)+k])
                                            h_mag += home_site[5]*neighbor_site[5]* float(js[len(beg_rule_list) + len(cluster_rule_list)+k])
                                if j_rule_list[k].neighbor_arrangement == 'COMB':
                                    if neighbor_site[4] in j_rule_list[k].neighbor_atom_list:
                                        h += home_site[5]*neighbor_site[5] * float(js[len(beg_rule_list) + len(cluster_rule_list)+k])
                                        h_mag += home_site[5]*neighbor_site[5] * float(js[len(beg_rule_list) + len(cluster_rule_list)+k])
    return float(h),float(h_mag)

def flip_species(i,supercell_list):
    home_site = supercell_list[i]
    old_home_site = list(supercell_list[i])
    home_species = home_site[4]
    site_flipped = False
    while site_flipped == False:
        rand = np.random.randint(len(supercell_list)-1)
        neighbor_site = supercell_list[rand]
        old_neighbor_site = supercell_list[rand]
        if neighbor_site[4] != 0 and neighbor_site[4] != home_site[4]:
            home_site[4] = neighbor_site[4]
            neighbor_site[4] = home_species
            site_flipped = True
    supercell_list[i] = home_site
    supercell_list[rand] = neighbor_site
    return old_home_site,old_neighbor_site,rand

def flip_phase(i,supercell_list):
    home_site = supercell_list[i]
    old_home_site = list(supercell_list[i])
    old_phase = home_site[6]
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
    home_site[6] = phase
    supercell_list[i] = home_site
    return old_home_site

def flip_spin(i,supercell_list):
    home_site = supercell_list[i]
    old_home_site = deepcopy(supercell_list[i])
    same_spin = False
    while same_spin == False:
        rand = np.random.random()
        if rand <= 1/3.0:
            spin = -1
        elif rand > 1/3.0 and rand < 2/3.0:
            spin = 0
        elif rand >= 2/3.0:
            spin = 1
        if spin != home_site[5]:
            same_spin = True
    home_site[5] = spin
    supercell_list[i] = home_site
    return old_home_site
