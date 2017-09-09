_author__ = 'brian'
# This file handles taking the input from vasp and calculating the fit.
# The functions that are important for the actual calculating the fit
# have comments

import numpy as np
#import BEG
import clusters
import js
import m_structure
import os
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from sklearn.linear_model import Ridge
from sklearn import linear_model

def import_vasp(root_dir, output_dir):
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
            output_line = name + "\t" + u + "\t" + enrg + "\t" + a + "\t" + b + "\t" + c + "\n"
            output.write(output_line)
            index = 0
            for i in range(8, 8 + total_num):
                pos = contcar_lines[i].split()
                output_line = "\t" + str(i - 7) + "\t" + str(mag_list[index]) + "\t" + str(pos[0]) + "\t" + str(
                    pos[1]) + "\t" + str(pos[2]) + "\n"
                output.write(output_line)
                index += 1
    output.close()

def generate_m_structure(data_file, num_Cluster_rules, num_J_rules, aust_tol, spin_tol):
    m_struct_list = []
    data = open(data_file, 'r')
    lines = data.readlines()
    for i in range(len(lines)):
        if '#' in lines[i]:
            species = (lines[i].split(' '))
            species.pop(0)
            m_struct = m_structure.MStructureObj(lines[i+1], species, num_Cluster_rules, num_J_rules, aust_tol)
            for j in range(m_struct.num_Atoms):
                atom_data = lines[i + j + 2]
                m_struct.set_atom_properties(j, atom_data, spin_tol)
            if m_struct.phase_name != 'prem':
                m_struct_list.append(m_struct)
    return m_struct_list

def write_structures_processedvasp(structures,data_file_pp):
    file = open(data_file_pp, 'w')
    for i in range(len(structures)):
        file.write('# Ni Mn In \n')            # hard coding!
        mat = structures[i]
        out = [mat.composition, mat.name, mat.enrg, mat.phase_name, mat.LCs]
        for j in range(len(out)):
            sums = str(out[j])
            file.write(sums.ljust(20))
        file.write("\n")
        for j in range(mat.num_Atoms):
            file.write("\t")
            out = [mat.basis[j].atom_index,mat.basis[j].species,mat.basis[j].spin,mat.basis[j].a_pos,mat.basis[j].b_pos,mat.basis[j].c_pos]
            for k in range(len(out)):
                sums = str(out[k])
                file.write(sums.ljust(17))
            file.write("\n")
    file.close()

# This function sweeps through the VASP data and determines active interactions and clusters
# for each atom site and group of sites
def calculate_sums(m_structure_list, cluster_rule_list, j_rule_list, spin_tol):
    for i in range(len(m_structure_list)):          # maybe this loop should be outside, just have this calc sums for a given structure
        m_structure_list[i].create_supercell(spin_tol)
        m_structure_list[i].calculate_distances()
        m_structure_list[i].calculate_minimums()
        for j in range(m_structure_list[i].num_Atoms):
            for k in range(len(m_structure_list[i].basis)):
#                # Calc BEG sums
#                for l in range(len(beg_rule_list)):
#                    if m_structure_list[i].basis[j].species in beg_rule_list[l].home_atom_list:
#                        if m_structure_list[i].basis[k].species in beg_rule_list[l].neighbor_atom_list:
#                            if m_structure_list[i].distances[j, k] == m_structure_list[i].mins[j, beg_rule_list[l].neighbor_order - 1]:
#                                if m_structure_list[i].phase_name == beg_rule_list[l].phase:
#                                    if m_structure_list[i].composition == beg_rule_list[l].composition:
#                                        if beg_rule_list[l].neighbor_arrangement == 'COMB':
#                                            m_structure_list[i].BEG_sums[l] += 1
#                                        if beg_rule_list[l].neighbor_arrangement == 'PERM':
#                                            if m_structure_list[i].basis[k].species != m_structure_list[i].basis[j].species:
#                                                m_structure_list[i].BEG_sums[l] += 1
#        ##########################BEG V2#########################
#        # for j in range(m_structure_list[i].num_Atoms):
#        #     for k in range(len(m_structure_list[i].basis)):
#        #         # Calc BEG sums
#        #         for l in range(len(beg_rule_list)):
#        #             if m_structure_list[i].basis[j].species in beg_rule_list[l].home_atom_list:
#        #                 if m_structure_list[i].basis[k].species in beg_rule_list[l].neighbor_atom_list:
#        #                     if m_structure_list[i].distances[j, k] == m_structure_list[i].mins[j, beg_rule_list[l].neighbor_order - 1]:
#        #                         if m_structure_list[i].phase_name == beg_rule_list[l].phase:
#        #                             if m_structure_list[i].composition == beg_rule_list[l].composition:
#        #                                 m_structure_list[i].BEG_sums[l] += -1
#        #########################################################
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

def summarize_fitting_structures(structures):
    path = 'summary_fitting_structures'
    file = open(path, 'w')
    file.write("NAME".ljust(15) + "PHASE".ljust(7) + "LCONST".ljust(15) + "MAG".ljust(6) + "ENERG".ljust(17) + "SUMS->\n")
    for i in range(len(structures)):
        mat = structures[i]
        out = [mat.enrg, mat.Cluster_sums, mat.J_sums]
        file.write(mat.name.ljust(15) + mat.phase_name.ljust(7) + str(round(mat.LCs[0],2)).ljust(5) + str(round(mat.LCs[1],2)).ljust(5) + str(round(mat.LCs[2],2)).ljust(5) + mat.mag_phase.ljust(7))
        for j in range(len(out)):
            sums = str(out[j])
            if j == 0:
                file.write(sums.ljust(17))
            else:
                file.write(sums.ljust(7))
        file.write("\n")
    file.close()

def check_duplicate_structures(structures):
    dupl = 'False';
    for i in range(len(structures)):
        for j in range(i+1,len(structures)):
            if (structures[i].Cluster_sums == structures[j].Cluster_sums):
                if (structures[i].J_sums == structures[j].J_sums):
                    dupl = 'True';
                    print('Duplicate fitting structure found: ',structures[i].name,',',structures[j].name,'\n')
    return dupl



