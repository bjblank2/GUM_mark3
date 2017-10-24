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
#from scipy.optimize import least_squares
#from sklearn.linear_model import Ridge
#from sklearn import linear_model

def generate_m_structure(data_file, num_Cluster_rules, num_J_rules, aust_tol, spin_style, spin_tol, Cluster_rules, J_rules):
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
                m_struct.set_atom_properties(j, atom_data, spin_style, spin_tol)
            calculate_sums(m_struct, Cluster_rules, J_rules, spin_style, spin_tol)
            # NEED TO ASSIGN SPIN TYPE HERE BASED ON CALCULATED SUMS - need a function for this.
            if m_struct.phase_name != 'prem':
                if check_duplicate_structures(m_struct,m_struct_list)=='False':
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
def calculate_sums(m_struct, cluster_rule_list, j_rule_list, spin_style, spin_tol):
        m_struct.create_supercell(spin_style, spin_tol)
        m_struct.calculate_distances()
        m_struct.calculate_minimums()
        for j in range(m_struct.num_Atoms):
            for k in range(len(m_struct.basis)):
                # Calc Cluster sums
                for l in range(len(cluster_rule_list)):
                    if m_struct.basis[j].species in cluster_rule_list[l].home_atom_list:
                        if m_struct.distances[j, k] == m_struct.mins[j, cluster_rule_list[l].neighbor_order - 1]:
                            if m_struct.check_plane(j, k) == cluster_rule_list[l].plane or m_struct.check_plane(j, k) == 'ALL':
                                if m_struct.phase_name == cluster_rule_list[l].phase:
                                    if cluster_rule_list[l].neighbor_arrangement == 'COMB':
                                        if m_struct.basis[k].species in cluster_rule_list[l].neighbor_atom_list:
                                            m_struct.Cluster_sums[l] += 1
                                    if cluster_rule_list[l].neighbor_arrangement == 'PERM':
                                        if m_struct.basis[k].species in cluster_rule_list[l].neighbor_atom_list:
                                            if m_struct.basis[k].species != m_struct.basis[j].species:
                                                m_struct.Cluster_sums[l] += 1
                # Calc J sums
                for l in range(len(j_rule_list)):
                    if m_struct.basis[j].species in j_rule_list[l].home_atom_list:
                        if m_struct.distances[j, k] == m_struct.mins[j, j_rule_list[l].neighbor_order - 1]:
                            if m_struct.check_plane(j, k) == j_rule_list[l].plane or m_struct.check_plane(j, k) == 'ALL':
                                if m_struct.phase_name == j_rule_list[l].phase:
                                    if j_rule_list[l].neighbor_arrangement == 'COMB':
                                        if m_struct.basis[k].species in j_rule_list[l].neighbor_atom_list:
                                            m_struct.J_sums[l] += m_struct.basis[j].spin * m_struct.basis[k].spin
                                    if j_rule_list[l].neighbor_arrangement == 'PERM':
                                        if m_struct.basis[k].species in j_rule_list[l].neighbor_atom_list:
                                            if m_struct.basis[k].species != m_struct.basis[j].species:
                                                m_struct.J_sums[l] += m_struct.basis[j].spin * m_struct.basis[k].spin

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

def check_duplicate_structures(structure,structure_list):
    dupl = 'False';
    for i in range(len(structure_list)):
        if (structure.Cluster_sums == structure_list[i].Cluster_sums):
            if (structure.J_sums == structure_list[i].J_sums):
                dupl = 'True';
                print('Duplicate fitting structure found: ',structure.name,'(energy =',structure.enrg,'eV), ',structure_list[i].name,'(energy =',structure_list[i].enrg,'eV)')
    return dupl



