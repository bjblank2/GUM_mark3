_author__ = 'brian'
# This file handles taking the input from vasp and calculating the fit.
# The functions that are important for the actual calculating the fit
# have comments

import numpy as np
import BEG
import clusters
import js
import m_structure
import os
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from sklearn.linear_model import Ridge
from sklearn import linear_model

def import_data(number_of_species, root_dir, output_dir):
    output = open(output_dir, 'w')
    for subdir, dirs, files in os.walk(root_dir):
        contcar_lines = []
        outcar_lines = []
        flag = 0
        for file in files:
            if 'CONTCAR' in files and 'OUTCAR' in files:
                name = subdir.strip(root_dir)
                flag = 1
                if file == "CONTCAR":
                    contcar = open(subdir + '/' + file, 'r')
                    contcar_lines = contcar.readlines()
                    contcar.close()
                if file == "OUTCAR":
                    outcar = open(subdir + '/' + file, 'r')
                    outcar_lines = outcar.readlines()
                    outcar_len = len(outcar_lines)
                    outcar.close()
            else:
                flag = 0
        if len(contcar_lines) > 0 and len(outcar_lines) > 0 and flag == 1:
            lc = float(contcar_lines[1])
            a = contcar_lines[2].split()
            a = float(a[0]) * lc
            b = contcar_lines[3].split()
            b = float(b[1]) * lc
            c = contcar_lines[4].split()
            c = float(c[2]) * lc
            composition = [0] * number_of_species
            comp = contcar_lines[6].split()
            for i in range(len(comp)):
                composition[i] = comp[i]
            lcs = [a, b, c]
            C = max(lcs)
            find_c = [abs(a - b), abs(a - c), abs(b - c)]
            min_ind = find_c.index(min(find_c))
            if min_ind == 0:
                C = c
                A = a
                B = b
                u = c / a
            elif min_ind == 1:
                C = b
                A = a
                B = c
                u = b / a
            elif min_ind == 2:
                C = a
                A = c
                B = b
                u = a / c
            if u > 1.1:
                phase = "mart"
            if u <= 1.1 and u >= 0.9:
                phase = "aust"
            if u < 0.9:
                phase = "pm"
            output.write("#\t")
            total_num = 0
            for i in range(len(composition)):
                output.write(str(composition[i]) + "\t")
                total_num += int(composition[i])
            mag_list = [0] * (total_num)
            for i in range(outcar_len):
                if "TOTEN" in outcar_lines[i]:
                    enrg = outcar_lines[i].split()
                    enrg = float(enrg[4])
                if "magnetization (x)" in outcar_lines[i]:
                    for j in range(total_num):
                        mag = outcar_lines[i + j + 4].split()
                        mag_list[j] = mag[4]
            u = str(u)
            enrg = str(enrg)
            a = str(a)
            b = str(b)
            c = str(c)
            output_line = phase + "\tNA\t" + name + "\t" + u + "\t" + enrg + "\t" + a + "\t" + b + "\t" + c + "\n"
            output.write(output_line)
            index = 0
            for i in range(8, 8 + total_num):
                pos = contcar_lines[i].split()
                output_line = "\t" + str(i - 7) + "\t" + str(mag_list[index]) + "\t" + str(pos[0]) + "\t" + str(
                    pos[1]) + "\t" + str(pos[2]) + "\n"
                output.write(output_line)
                index += 1
    output.close()


def write_beg_rules(rule_file):
    output = open(rule_file, 'w')
    new_rule = True
    while new_rule == True:
        BEG_rule = BEG.BEGObj()
        BEG_rule.create_rule()
        output.write('# ' + str(BEG_rule.name) + '\n')
        output.write(str(BEG_rule.neighbor_order) + '\n')
        output.write(str(BEG_rule.neighbor_arrangement) + '\n')
        output.write(str(BEG_rule.home_atom_list) + '\n')
        output.write(str(BEG_rule.neighbor_atom_list) + '\n')
        output.write(str(BEG_rule.phase) + '\n')
        output.write(str(BEG_rule.plane) + '\n')
        if input('Add another rule? (Y/N):  ') == 'N':
            new_rule = False
    output.close()


def read_beg_rules(rule_file):
    input_file = open(rule_file, 'r')
    lines = input_file.readlines()
    BEG_rule_list = []
    for i in range(len(lines)):
        if '#' in lines[i]:
            BEG_rule = BEG.BEGObj()
            name = lines[i]
            BEG_rule.set_name(name.strip('# '))
            BEG_rule.set_neighbor_order(int(lines[i + 1]))
            BEG_rule.set_neighbor_arrangement(lines[i + 2].strip())
            BEG_rule.set_home_atom_list(lines[i + 3].split())
            BEG_rule.set_neighbor_atom_list(lines[i + 4].split())
            BEG_rule.set_phase(lines[i + 5].strip())
            BEG_rule.set_plane(lines[i + 6].strip())
            BEG_rule.set_composition(lines[i + 7].split())
            BEG_rule_list.append(BEG_rule)
    input_file.close()
    return BEG_rule_list


def write_cluster_rules(rule_file):
    output = open(rule_file, 'w')
    new_rule = True
    while new_rule == True:
        Cluster_rule = clusters.ClusterObj()
        Cluster_rule.create_rule()
        output.write('# ' + str(Cluster_rule.name) + '\n')
        output.write(str(Cluster_rule.neighbor_order) + '\n')
        output.write(str(Cluster_rule.neighbor_arrangement) + '\n')
        output.write(str(Cluster_rule.home_atom_list) + '\n')
        output.write(str(Cluster_rule.neighbor_atom_list) + '\n')
        output.write(str(Cluster_rule.phase) + '\n')
        output.write(str(Cluster_rule.plane) + '\n')
        output.write(str(Cluster_rule.composition) + '\n')
        if input('Add another rule? (Y/N):  ') == 'N':
            new_rule = False
    output.close()


def read_cluster_rules(rule_file):
    input_file = open(rule_file, 'r')
    lines = input_file.readlines()
    Cluster_rule_list = []
    for i in range(len(lines)):
        if '#' in lines[i]:
            Cluster_rule = clusters.ClusterObj()
            name = lines[i]
            Cluster_rule.set_name(name.strip('# '))
            Cluster_rule.set_neighbor_order(int(lines[i + 1]))
            Cluster_rule.set_neighbor_arrangement(lines[i + 2].strip())
            Cluster_rule.set_home_atom_list(lines[i + 3].split())
            Cluster_rule.set_neighbor_atom_list(lines[i + 4].split())
            Cluster_rule.set_phase(lines[i + 5].strip())
            Cluster_rule.set_plane(lines[i + 6].strip())
            Cluster_rule_list.append(Cluster_rule)
    input_file.close()
    return Cluster_rule_list


def write_j_rules(rule_file):
    output = open(rule_file, 'w')
    new_rule = True
    while new_rule == True:
        J_rule = js.JObj()
        J_rule.create_rule()
        output.write('# ' + str(J_rule.name) + '\n')
        output.write(str(J_rule.neighbor_order) + '\n')
        output.write(str(J_rule.neighbor_arrangement) + '\n')
        output.write(str(J_rule.home_atom_list) + '\n')
        output.write(str(J_rule.neighbor_atom_list) + '\n')
        output.write(str(J_rule.phase) + '\n')
        output.write(str(J_rule.plane) + '\n')
        if input('Add another rule? (Y/N):  ') == 'N':
            new_rule = False
    output.close()


def read_j_rules(rule_file):
    input_file = open(rule_file, 'r')
    lines = input_file.readlines()
    J_rule_list = []
    for i in range(len(lines)):
        if '#' in lines[i]:
            J_rule = js.JObj()
            name = lines[i]
            J_rule.set_name(name.strip('# '))
            J_rule.set_neighbor_order(int(lines[i + 1]))
            J_rule.set_neighbor_arrangement(lines[i + 2].strip())
            J_rule.set_home_atom_list(lines[i + 3].split())
            J_rule.set_neighbor_atom_list(lines[i + 4].split())
            J_rule.set_phase(lines[i + 5].strip())
            J_rule.set_plane(lines[i + 6].strip())
            J_rule_list.append(J_rule)
    input_file.close()
    return J_rule_list


def read_Js(num_Js):
    js_file = open('output','r')
    Js = []
    lines = js_file.readlines()
    for i in range(num_Js):
        line = lines[i+1]
        line = line.split()
        Js.append(float(line[2]))
    return Js


def read_m_structure_data(data_file, num_species, num_BEG_rules, num_Cluster_rules, num_J_rules):
    m_struct_list = []
    data = open(data_file, 'r')
    lines = data.readlines()
    for i in range(len(lines)):
        if '#' in lines[i]:
            m_struct = m_structure.MStructureObj(lines[i], num_species, num_BEG_rules, num_Cluster_rules, num_J_rules)
            for j in range(m_struct.num_Atoms):
                atom_data = lines[i + j + 1]
                m_struct.set_atom_properties(j, atom_data)
            m_struct_list.append(m_struct)
    return m_struct_list

# This function sweeps through the VASP data and determines active interactions and clusters
# for each atom site and group of sites
def calculate_sums(m_structure_list, beg_rule_list, cluster_rule_list, j_rule_list):
    for i in range(len(m_structure_list)):
        m_structure_list[i].create_super_cell()
        m_structure_list[i].calculate_distances()
        m_structure_list[i].calculate_minimums()
        for j in range(m_structure_list[i].num_Atoms):
            for k in range(len(m_structure_list[i].basis)):
                # Calc BEG sums
                for l in range(len(beg_rule_list)):
                    if m_structure_list[i].basis[j].species in beg_rule_list[l].home_atom_list:
                        if m_structure_list[i].basis[k].species in beg_rule_list[l].neighbor_atom_list:
                            if m_structure_list[i].distances[j, k] == m_structure_list[i].mins[j, beg_rule_list[l].neighbor_order - 1]:
                                if m_structure_list[i].phase_name == beg_rule_list[l].phase:
                                    if m_structure_list[i].composition == beg_rule_list[l].composition:
                                        if beg_rule_list[l].neighbor_arrangement == 'COMB':
                                            m_structure_list[i].BEG_sums[l] += 1
                                        if beg_rule_list[l].neighbor_arrangement == 'PERM':
                                            if m_structure_list[i].basis[k].species != m_structure_list[i].basis[j].species:
                                                m_structure_list[i].BEG_sums[l] += 1
        ##########################BEG V2#########################
        # for j in range(m_structure_list[i].num_Atoms):
        #     for k in range(len(m_structure_list[i].basis)):
        #         # Calc BEG sums
        #         for l in range(len(beg_rule_list)):
        #             if m_structure_list[i].basis[j].species in beg_rule_list[l].home_atom_list:
        #                 if m_structure_list[i].basis[k].species in beg_rule_list[l].neighbor_atom_list:
        #                     if m_structure_list[i].distances[j, k] == m_structure_list[i].mins[j, beg_rule_list[l].neighbor_order - 1]:
        #                         if m_structure_list[i].phase_name == beg_rule_list[l].phase:
        #                             if m_structure_list[i].composition == beg_rule_list[l].composition:
        #                                 m_structure_list[i].BEG_sums[l] += -1
        #########################################################
                # Calc Cluster sums
                for l in range(len(cluster_rule_list)):
                    if m_structure_list[i].basis[j].species in cluster_rule_list[l].home_atom_list:
                        if m_structure_list[i].distances[j, k] == m_structure_list[i].mins[j, cluster_rule_list[l].neighbor_order - 1]:
                            if m_structure_list[i].check_plane(j, k) == cluster_rule_list[l].plane or m_structure_list[i].check_plane(j, k) == 'ALL':
                                if m_structure_list[i].phase_name == cluster_rule_list[l].phase:
                                    if cluster_rule_list[l].neighbor_arrangement == 'COMB':
                                        if m_structure_list[i].basis[k].species in cluster_rule_list[l].neighbor_atom_list:
                                            m_structure_list[i].Cluster_sums[l] += 1
                                    if cluster_rule_list[l].neighbor_arrangement == 'PERM':
                                        if m_structure_list[i].basis[k].species in cluster_rule_list[l].neighbor_atom_list:
                                            if m_structure_list[i].basis[k].species != m_structure_list[i].basis[j].species:
                                                m_structure_list[i].Cluster_sums[l] += 1
                # Calc J sums
                for l in range(len(j_rule_list)):
                    if m_structure_list[i].basis[j].species in j_rule_list[l].home_atom_list:
                        if m_structure_list[i].distances[j, k] == m_structure_list[i].mins[j, j_rule_list[l].neighbor_order - 1]:
                            if m_structure_list[i].check_plane(j, k) == j_rule_list[l].plane or m_structure_list[i].check_plane(j, k) == 'ALL':
                                if m_structure_list[i].phase_name == j_rule_list[l].phase:
                                    if j_rule_list[l].neighbor_arrangement == 'COMB':
                                        if m_structure_list[i].basis[k].species in j_rule_list[l].neighbor_atom_list:
                                            m_structure_list[i].J_sums[l] += m_structure_list[i].basis[j].spin * m_structure_list[i].basis[k].spin
                                    if j_rule_list[l].neighbor_arrangement == 'PERM':
                                        if m_structure_list[i].basis[k].species in j_rule_list[l].neighbor_atom_list:
                                            if m_structure_list[i].basis[k].species != m_structure_list[i].basis[j].species:
                                                m_structure_list[i].J_sums[l] += m_structure_list[i].basis[j].spin * m_structure_list[i].basis[k].spin


def find_weights(m_structure_list, compositions, tk):
    for i in range(len(compositions)):
        enrgs = []
        for j in range(len(m_structure_list)):
            if compositions[i] == m_structure_list[j].composition[1]:
                enrgs.append(m_structure_list[j].enrg)
        minimum = min(enrgs)
        for j in range(len(m_structure_list)):
            if compositions[i] == m_structure_list[j].composition[1]:
                m_structure_list[j].weight = np.exp(-1.0 * abs(minimum - m_structure_list[j].enrg) / tk) ** (0.5)


def find_weights_2(m_structure_list, compositions,limit):
    for i in range(len(compositions)):
        enrgs = []
        for j in range(len(m_structure_list)):
            if compositions[i] == m_structure_list[j].composition[1]:
                if m_structure_list[j].mag_phase != "pera" and m_structure_list[j].phase_name != "pm" and m_structure_list[j].enrg <= limit:
                    enrgs.append(m_structure_list[j].enrg)
        minimum = min(enrgs)
        maximum = max(enrgs)
        cutoff = (maximum-minimum)*.25
        for j in range(len(m_structure_list)):
            if compositions[i] == m_structure_list[j].composition[1]:
                if m_structure_list[j].enrg > minimum+cutoff:
                    m_structure_list[j].weight = .7


def do_weighted_ls(m_structure_list, limit):
    a = []
    b = []
    for i in range(len(m_structure_list)):
        mat = m_structure_list[i]
        if mat.mag_phase != "pera" and mat.phase_name != "pm" and mat.enrg <= limit:
            # if mat.phase_name != "pm" and mat.enrg <= limit:
            row = mat.BEG_sums + mat.Cluster_sums + mat.J_sums
            for j in range(len(row)):
                row[j] *= mat.weight
            a.append(row)
            b.append(mat.enrg * mat.weight)
    a = np.matrix(a)
    b = np.transpose(np.matrix(b))
    Js = np.linalg.lstsq(a, b)[0]
    return Js


def CV_score(Js,m_structure_list):
    CV2 = 0
    energies = []
    for i in range(len(m_structure_list)):
        mat = m_structure_list[i]
        if mat.mag_phase != "pera" and mat.phase_name != "pm":
            energies.append(m_structure_list[i].enrg)
    x_rows = len(energies)
    x_colm = m_structure_list[0].num_beg_rules+m_structure_list[0].num_cluster_rules+m_structure_list[0].num_j_rules
    x_matrix = np.zeros((x_rows,x_colm))
    x_matrix = np.matrix(x_matrix)
    i = 0
    for itter in range(len(m_structure_list)):
        mat = m_structure_list[itter]
        if mat.mag_phase != "pera" and mat.phase_name != "pm":
            for j in range(m_structure_list[0].num_beg_rules):
                x_matrix[i,j] = m_structure_list[i].BEG_sums[j]
            for k in range(m_structure_list[0].num_cluster_rules):
                x_matrix[i,j+k] = m_structure_list[i].Cluster_sums[k]
            for l in range(m_structure_list[0].num_j_rules):
                x_matrix[i,j+k+l] = m_structure_list[i].J_sums[l]
            i += 1
    en_calc = np.zeros(len(energies))
    for i in range(x_rows):
        for j in range(x_colm):
            en_calc[i] += x_matrix[i,j]*Js[j]
    for i in range(len(energies)):
        x_i = np.matrix(x_matrix[i,:])
        x_trans = x_matrix.getT()
        x_it = x_i.getT()
        CV2 += ((energies[i]-en_calc[i])/(1-x_i*np.linalg.inv(x_trans*x_matrix+np.matrix(np.eye(26,26)*.05))*x_it))**2
    CV2 *= 1/len(m_structure_list)
    return CV2


def CV_score2(m_structure_list):
    energies = []
    for i in range(len(m_structure_list)):
        mat = m_structure_list[i]
        if mat.mag_phase != "pera" and mat.phase_name != "pm":
            energies.append(m_structure_list[i].enrg)
    x_rows = len(energies)
    x_colm = m_structure_list[0].num_beg_rules+m_structure_list[0].num_cluster_rules+m_structure_list[0].num_j_rules
    x_matrix = np.zeros((x_rows,x_colm))
    x_matrix = np.matrix(x_matrix)
    i = 0
    for itter in range(len(m_structure_list)):
        mat = m_structure_list[itter]
        if mat.mag_phase != "pera" and mat.phase_name != "pm":
            for j in range(m_structure_list[0].num_beg_rules):
                x_matrix[i,j] = m_structure_list[i].BEG_sums[j]
            for k in range(m_structure_list[0].num_cluster_rules):
                x_matrix[i,j+k] = m_structure_list[i].Cluster_sums[k]
            for l in range(m_structure_list[0].num_j_rules):
                x_matrix[i,j+k+l] = m_structure_list[i].J_sums[l]
            i += 1
    x_cross = np.zeros((x_rows-1,x_colm))
    CV2 = 0
    for i in range(len(energies)):
        inc = 0
        eng_cross = []
        for j in range(len(energies)):
            if i != j:
                eng_cross.append(energies[j])
                x_cross[inc,:] = x_matrix[j,:]
                inc += 1
        b = np.transpose(np.matrix(eng_cross))
        a = np.matrix(x_cross)
        Js = np.linalg.lstsq(a, b)[0]
        Ei = 0
        for j in range(len(Js)):
            Ei += x_matrix[i,j]*Js[j]
        CV2 += (energies[i]-Ei)**2
    CV2 *= 1/len(m_structure_list)
    return CV2

# Calculates the fitting paramiters using a ridge regression with autocorrelation
def ridge_simple(m_structure_list,alpha):
    a = []
    b = []
    for i in range(len(m_structure_list)):
        mat = m_structure_list[i]
        if mat.mag_phase != "pera" and mat.phase_name != "pm":
            # if mat.phase_name != "pm" and mat.enrg <= limit:
            row = mat.BEG_sums + mat.Cluster_sums + mat.J_sums
            for j in range(len(row)):
                row[j] *= mat.weight
            a.append(row)
            b.append(mat.enrg * mat.weight)
    a = np.matrix(a)
    b = np.transpose(np.matrix(b))
    ridge_fit = linear_model.RidgeCV(alphas=[.01,.05,.1,.15,.2,.25,.3,1,2,3,4,5,10],fit_intercept=False, )
    ridge_fit.fit(a,b)
    Js = ridge_fit.coef_
    JS_list = []
    Js_0 = Js[0]
    for i in range(len(Js_0)):
        Js_0
        JS_list.append(Js_0[i])
    #print(Js)
    #print(ridge_fit.predict(a))
    print(ridge_fit.intercept_)
    return JS_list


def ridge_optimized_fit(m_structure_list,a_min,a_max,step):
    alpha = a_min
    old_CV2_min = 10000000
    alpha_opt = 0
    while alpha < a_max:
        energies = []
        for i in range(len(m_structure_list)):
            mat = m_structure_list[i]
            if mat.mag_phase != "pera" and mat.phase_name != "pm":
                energies.append(m_structure_list[i].enrg)
        x_rows = len(energies)
        x_colm = m_structure_list[0].num_beg_rules+m_structure_list[0].num_cluster_rules+m_structure_list[0].num_j_rules
        x_matrix = np.zeros((x_rows,x_colm))
        x_matrix = np.matrix(x_matrix)
        i = 0
        for itter in range(len(m_structure_list)):
            mat = m_structure_list[itter]
            if mat.mag_phase != "pera" and mat.phase_name != "pm":
                for j in range(m_structure_list[0].num_beg_rules):
                    x_matrix[i,j] = m_structure_list[i].BEG_sums[j]
                for k in range(m_structure_list[0].num_cluster_rules):
                    x_matrix[i,j+k] = m_structure_list[i].Cluster_sums[k]
                for l in range(m_structure_list[0].num_j_rules):
                    x_matrix[i,j+k+l] = m_structure_list[i].J_sums[l]
                i += 1
        x_cross = np.zeros((x_rows-1,x_colm))
        CV2 = 0
        for i in range(len(energies)):
            inc = 0
            eng_cross = []
            for j in range(len(energies)):
                if i != j:
                    eng_cross.append(energies[j])
                    x_cross[inc,:] = x_matrix[j,:]
                    inc += 1
            b = np.transpose(np.matrix(eng_cross))
            a = np.matrix(x_cross)
            r_fit = Ridge(alpha=alpha)
            r_fit.fit(a,b)
            Js = r_fit.get_params()
            Ei = 0
            for j in range(len(Js)):
                Ei += x_matrix[i,j]*Js[j]
            CV2 += (energies[i]-Ei)**2
        CV2 *= 1/len(m_structure_list)
        if CV2 < old_CV2_min:
            old_CV2_min = CV2
            alpha_opt = alpha
    r_fit = Ridge(alpha=alpha_opt)
    r_fit.fit(a,b)
    Js = r_fit.get_params()
    print(alpha_opt)
    print(CV2)
    return Js


def write_data(structures, limit, Js):
    path = 'data'
    file = open(path, 'w')
    file.write("NAME".ljust(15) + "PHASE".ljust(7) + "MAG".ljust(6) + "ENERG".ljust(17) + "SUMS->\n")
    for i in range(len(structures)):
        mat = structures[i]
        if mat.phase_name != "pmmmm":
            out = [mat.enrg, mat.BEG_sums, mat.Cluster_sums, mat.J_sums]
            file.write(mat.name.ljust(15) + mat.phase_name.ljust(7) + mat.mag_phase.ljust(7))
            for j in range(len(out)):
                sums = str(out[j])
                if j == 0:
                    file.write(sums.ljust(17))
                else:
                    file.write(sums.ljust(7))
            file.write("\n")
    file.close()


def write_output(structures, beg_list, clusters_list, j_list, Js, limit):
    label = []
    for i in range(len(beg_list)):
        label.append(beg_list[i].name)
    for i in range(len(clusters_list)):
        label.append(clusters_list[i].name)
    for i in range(len(j_list)):
        label.append(j_list[i].name)
    path = 'output'
    file = open(path, 'w')
    file.write("Fitting Parameters\n")
    for i in range(len(Js)):
        line = label[i].strip() + " = " + str(Js[i]) + "\n"
        line = line.replace("[", "")
        line = line.replace("]", "")
        file.write(line)
    file.write("\n")
    file.write("Original Enrg\tNew Enrg\n")
    for i in range(len(structures)):
        mat = structures[i]
        if mat.phase_name != "pmmmm":
            file.write(str(mat.composition[1]) + "    " + str(mat.enrg).ljust(16) + "    ")
            new_enrg = 0
            for j in range(len(beg_list)):
                new_enrg += Js[j] * structures[i].BEG_sums[j]
            for k in range(len(clusters_list)):
                new_enrg += Js[len(beg_list)+k] * structures[i].Cluster_sums[k]
            for l in range(len(j_list)):
                new_enrg += Js[len(beg_list) + len(clusters_list) + l] * structures[i].J_sums[l]
            line = str(new_enrg)
            line = line.replace("[", "")
            line = line.replace("]", "")
            file.write(line + "    " + mat.mag_phase + "    " + mat.phase_name + "     " + str(mat.weight) + "\n")
    file.close()


def r_function(x,t,y):
    r = np.zeros((1,len(y)))
    r = np.matrix.tolist(r)
    r = r[0]
    for i in range(len(y)):
        sums = t[i]
        for j in range(len(sums)):
            r[i] += x[j]*sums[j]
        r[i] -= y[i]
    return r


def do_robust_ls(M_structures):
    t_train = []
    y_train = []
    for i in range(len(M_structures)):
        line = list(M_structures[i].BEG_sums+M_structures[i].Cluster_sums+M_structures[i].J_sums)
        for j in range(len(line)):
            line[j] *= M_structures[i].weight
        t_train.append(line)
        y_train.append(M_structures[i].enrg * M_structures[i].weight)
    x0 = np.ones((len(line)))*1
    res_lsq = least_squares(r_function, x0, loss='soft_l1', f_scale=0.05, args=(t_train, y_train))
    return res_lsq.x


def ransac(M_structures,error_cutoff,good_fit_cutoff,iterations):
    best_error = 100000000000
    modle_change_count = 0
    for iterations in range (iterations):
        candidate_list = []
        rand_int_list = []
        for i in range(26):
            rand_int = np.random.randint(0,len(M_structures)-1)
            candidate_list.append(M_structures[rand_int])
            rand_int_list.append(rand_int)
        #candidate_model = do_robust_ls(candidate_list)
        candidate_model = do_weighted_ls(candidate_list,500)
        candidate_inliers = []
        for i in range(len(M_structures)):
            if i not in rand_int_list:
                energy = 0
                sums = list(M_structures[i].BEG_sums+M_structures[i].Cluster_sums+M_structures[i].J_sums)
                for j in range(len(candidate_model)):
                    energy += candidate_model[j]*sums[j]
                error = abs(energy-M_structures[i].enrg)
                if error <= error_cutoff:
                    candidate_inliers.append(M_structures[i])
        if len(candidate_inliers)+len(candidate_list) >= good_fit_cutoff:
            new_candidate_list = list(candidate_list+candidate_inliers)
            #new_candidate_model = do_robust_ls(new_candidate_list)
            new_candidate_model = do_weighted_ls(M_structures, 500)
            new_error = 0
            for i in range(len(new_candidate_list)):
                new_energy = 0
                new_sums = list(M_structures[i].BEG_sums+M_structures[i].Cluster_sums+M_structures[i].J_sums)
                for j in range(len(new_candidate_model)):
                    new_energy += new_candidate_model[j]*new_sums[j]
                new_error += abs(new_energy-new_candidate_list[i].enrg)/len(new_candidate_list)
            if new_error < best_error:
                best_error = new_error
                best_model = new_candidate_model
                modle_change_count += 1
    return best_model


def ransacom(M_structures,error_cutoff,good_fit_cutoff,iterations):
    best_error = 100000000000
    modle_change_count = 0
    for iterations in range (iterations):
        candidate_list = []
        rand_int_list = []
        for i in range(26):
            rand_int = np.random.randint(0,len(M_structures)-1)
            if rand_int not in rand_int_list:
                rand_int_list.append(rand_int)
                M_structures[rand_int].weight = 1.0
                candidate_list.append(M_structures[rand_int])
        #candidate_model = do_robust_ls(candidate_list)
        candidate_model = do_weighted_ls(candidate_list,500)
        candidate_inliers = []
        outlire_list = []
        for i in range(len(M_structures)):
            M_structures[i].weight = 1.0
            if i not in rand_int_list:
                energy = 0
                sums = list(M_structures[i].BEG_sums+M_structures[i].Cluster_sums+M_structures[i].J_sums)
                for j in range(len(candidate_model)):
                    energy += candidate_model[j]*sums[j]
                error = abs(energy-M_structures[i].enrg)
                if error <= error_cutoff:
                    candidate_inliers.append(M_structures[i])
                else:
                    M_structures[i].weight = .4
                    outlire_list.append(M_structures[i])

        if len(candidate_inliers)+len(candidate_list) >= good_fit_cutoff:
            new_candidate_list = list(candidate_list+candidate_inliers+outlire_list)
            #new_candidate_model = do_robust_ls(new_candidate_list)
            new_candidate_model = do_weighted_ls(M_structures, 500)
            new_error = 0
            for i in range(len(M_structures)):
                new_energy = 0
                new_sums = list(M_structures[i].BEG_sums+M_structures[i].Cluster_sums+M_structures[i].J_sums)
                for j in range(len(new_candidate_model)):
                    new_energy += new_candidate_model[j]*new_sums[j]
                new_error += abs(new_energy-new_candidate_list[i].enrg)/len(new_candidate_list)
            if new_error < best_error:
                best_outliers = outlire_list
                best_error = new_error
                best_model = new_candidate_model
                modle_change_count += 1
    for i in range(len(M_structures)):
        M_structures[i].weight = 1
    for i in range(len(best_outliers)):
        best_outliers[i].weight = .4
    best_model = do_weighted_ls(M_structures,500)
    print(modle_change_count)
    return best_model


def linearize(M_structures):
    e_comp0 = []
    e_comp50 = []
    for i in range(len(M_structures)):
        structure = M_structures[i]
        comp = structure.composition[2]/structure.composition[0]
        if comp == 0:
            e_comp0.append(structure.enrg)
        if comp == .5:
            e_comp50.append(structure.enrg)
    comp0_min = min(e_comp0)
    comp50_min = min(e_comp50)
    offset = (comp50_min-comp0_min)/.5
    for i in range(len(M_structures)):
        comp = M_structures[i].composition[2]/M_structures[i].composition[0]
        M_structures[i].enrg -= offset*comp + comp0_min


def scale(M_structures):
    min = 1000
    for i in range(len(M_structures)):
        structure = M_structures[i]
        if structure.enrg < min:
            min = structure.enrg
    for i in range(len(M_structures)):
        structure = M_structures[i]
        structure.enrg /= abs(min)


def plot_data():
    plt.rc('lines', linewidth=1)
    path = 'output'
    file = open(path, 'r')
    data = file.readlines()
    length = len(data)
    flag = 0
    for i in range(length):
        if flag == 1:
            enrgy = data[i].split()
            plt.plot([1, 2], [enrgy[1], enrgy[2]])
        if 'Original Enrg' in data[i]:
            flag = 1
    file.close()
    # plt.axis([0.5,2.5,-22.85,-22.4])
    plt.savefig('fit_line.png')


def plot_data2():
    colors = {'mart': 'g', 'aus': 'b', 'pre-mart': 'r'}
    markers = {'FM': 'o', 'AFM': 's', 'none': '^', 'FM/AFM': 'D'}
    labels = ['NiMn', 'Ni4Mn3In1', 'Ni2MnIn']
    actual_labels = ['NiMn', 'Ni$_4$Mn$_3$In', 'Ni$_2$MnIn']
    plt.rc('lines', linewidth=1)
    path = 'output'
    file = open(path, 'r')
    data = file.readlines()
    length = len(data)
    flag = 0
    x = 0
    x2_itter = -.5
    x4_itter = -.5
    x6_itter = -.5
    for i in range(length):
        if flag == 1:
            dat = data[i].split()
            if float(dat[0]) == 8:
                x = 2
                x2_itter += .08
                x_itter = x2_itter
            if float(dat[0]) == 6:
                x = 4
                x4_itter += .08
                x_itter = x4_itter
            if float(dat[0]) == 4:
                x = 8
                x6_itter += .08
                x_itter = x6_itter
            if dat[4] == "aust":
                c = 'b'
            if dat[4] == "mart":
                c = 'g'
            if dat[4] == "pm":
                c = 'r'
            if dat[3] == "fm":
                m = 'o'
            if dat[3] == "afm":
                m = 's'
            if dat[3] == "pera":
                m = '^'
            if dat[3] == "fi":
                m = 'D'
            if dat[3] == "NA":
                m = 'x'
            y = float(dat[1])
            plt.plot([x + x_itter, x + x_itter], [y, float(dat[2])], lw=1, color="k")
            plt.plot(x + x_itter, y, lw=0, markersize=8, marker=m, color=c)
            plt.plot(x + x_itter, float(dat[2]), lw=0, markersize=8, marker=".", color="r")
        if 'Original Enrg' in data[i]:
            flag = 1
    file.close()
    #    plt.xlim(0,8)
    #    plt.ylim(-1,10)
    plt.rc('lines', linewidth=1)
    plt.title("NiMn -- Ni$_2$MnIn Composition Energies", fontsize=24)
    plt.xlabel("Composition", fontsize=24)
    plt.ylabel("Energy above Hull (eV/16 atoms)", fontsize=24)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=12)
    #    plt.xticks([2,4,6],actual_labels, rotation='horizontal',fontsize=18)
    plt.savefig('Fit.png')


def plot_data3(M_structures, beg_list, clusters_list, j_list, Js, limit):
    colors = {'mart': 'g', 'aus': 'b', 'pre-mart': 'r'}
    markers = {'FM': 'o', 'AFM': 's', 'none': '^', 'FM/AFM': 'D'}
    labels = ['NiMn', 'Ni4Mn3In1', 'Ni2MnIn']
    actual_labels = ['Ni$_8$Mn$_8$', 'Ni$_4$Mn$_3$In$_1$', 'Ni$_4$Mn$_2$In$_2$']
    e_comp0 = []
    e_comp50 = []
    for i in range(len(M_structures)):
        if M_structures[i].phase_name != "pmmmm":
            structure = M_structures[i]
            comp = structure.composition[2]/structure.composition[0]
            if comp == 0:
                e_comp0.append(structure.enrg)
            if comp == .5:
                e_comp50.append(structure.enrg)
    comp0_min = min(e_comp0)
    comp50_min = min(e_comp50)
    offset = (comp50_min-comp0_min)/.5
    enrg_list = []
    for i in range(len(M_structures)):
        if M_structures[i].phase_name != "pmmmm":
            comp = M_structures[i].composition[2]/M_structures[i].composition[0]
            enrg_list.append(M_structures[i].enrg-(offset*comp + comp0_min))
    new_enrg_list = []
    for i in range(len(M_structures)):
        mat = M_structures[i]
        comp = mat.composition[2]/mat.composition[0]
        if mat.phase_name != "pmmmm":
            new_enrg = 0
            for j in range(len(beg_list)):
                new_enrg += Js[j] * M_structures[i].BEG_sums[j]
            for k in range(len(clusters_list)):
                new_enrg += Js[len(beg_list)+k] * M_structures[i].Cluster_sums[k]
            for l in range(len(j_list)):
                new_enrg += Js[len(beg_list) + len(clusters_list) + l] * M_structures[i].J_sums[l]
            new_enrg -= offset*comp + comp0_min
            new_enrg_list.append(new_enrg)
    x = 0
    x2_itter = -.5
    x4_itter = -.5
    x6_itter = -.5
    for i in range(len(enrg_list)):
       if M_structures[i].phase_name != "pmmmm":
            if int(M_structures[i].composition[1]) == 8:
                x = 1.5
                x2_itter += .08
                x_itter = x2_itter
            if int(M_structures[i].composition[1]) == 6:
                x = 4.5
                x4_itter += .08
                x_itter = x4_itter
            if int(M_structures[i].composition[1]) == 4:
                x = 9.5
                x6_itter += .08
                x_itter = x6_itter
            if M_structures[i].phase_name == "aust":
                c = 'b'
            if M_structures[i].phase_name == "mart":
                c = 'g'
            if M_structures[i].phase_name == "pm":
                c = 'r'
            if M_structures[i].mag_phase == "fm":
                m = 'o'
            if M_structures[i].mag_phase == "afm":
                m = 's'
            if M_structures[i].mag_phase == "pera":
                m = '^'
            if M_structures[i].mag_phase == "fi":
                m = 'D'
            if M_structures[i].mag_phase == "NA":
                m = 'x'
            y = float(enrg_list[i])
            if M_structures[i].phase_name != 'pm':
                #plt.plot([x + x_itter, x + x_itter], [y, float(new_enrg_list[i])], lw=1, color="k")
                plt.plot(x + x_itter, y, lw=0, markersize=8, marker=m, color=c)
                plt.plot(x + x_itter, float(new_enrg_list[i]), lw=0, markersize=8, marker=".", color="r")
    #plt.xlim(0,8)
    #plt.ylim(-1,9)
    plt.rc('lines', linewidth=1)
    plt.title("NiMn -- Ni$_2$MnIn Composition Energies", fontsize=24)
    plt.xlabel("Composition", fontsize=24)
    plt.ylabel("Energy above Hull (eV/16 atoms)", fontsize=24)
    #plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=12)
    #lt.xticks([2,5.5,10],actual_labels, rotation='horizontal',fontsize=18)
    plt.savefig('Fit.png')
