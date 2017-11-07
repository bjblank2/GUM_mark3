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
    norms = []
    for i in range(len(lines)):
        if '#' in lines[i]:
            species = (lines[i].split(' '))
            species.pop(0)
            m_struct = m_structure.MStructureObj(lines[i+1], species, num_Cluster_rules, num_J_rules, aust_tol)
            for j in range(m_struct.num_Atoms):
                atom_data = lines[i + j + 2]
                m_struct.set_atom_properties(j, atom_data, spin_style, spin_tol)
            #calculate_sums(m_struct, Cluster_rules, J_rules, spin_style, spin_tol)
            # NEED TO ASSIGN SPIN TYPE HERE BASED ON CALCULATED SUMS - need a function for this.
            j_count = calculate_sums_scaled(m_struct, Cluster_rules, J_rules, spin_style, spin_tol)
            if m_struct.phase_name != 'prem':
                if check_duplicate_structures(m_struct,m_struct_list)=='False':
                    m_struct_list.append(m_struct)
                    norms.append(j_count)
    return m_struct_list,norms

def write_structures_processedvasp(structures,data_file_pp,norms):
    file = open(data_file_pp, 'w')
    for i in range(len(structures)):
        file.write('# Ni Mn In \n')            # hard coding!
        mat = structures[i]
        out = [mat.composition, mat.name, mat.enrg, mat.phase_name, mat.LCs]
        for j in range(len(out)):
            sums = str(out[j])
            file.write(sums.ljust(20))
        file.write("".ljust(10))
        file.write(str(norms[i]).ljust(20))
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

# This function calculates the sum scaled by the number of Mg-Mg bonds
def calculate_sums_scaled(m_structure_list, cluster_rule_list, j_rule_list, spin_style, spin_tol):
    #for i in range(len(m_structure_list)):          # maybe this loop should be outside, just have this calc sums for a given structure
    j_count = [0]*len(j_rule_list)
    m_structure_list.create_supercell(spin_style, spin_tol)
    #print(i)

    m_structure_list.calculate_distances()
    m_structure_list.calculate_minimums()
    for j in range(m_structure_list.num_Atoms):
        for k in range(len(m_structure_list.basis)):
            # Calc Cluster sums
            for l in range(len(cluster_rule_list)):
                if m_structure_list.basis[j].species in cluster_rule_list[l].home_atom_list:
                    if m_structure_list.distances[j, k] == m_structure_list.mins[j, cluster_rule_list[l].neighbor_order - 1]:
                        if m_structure_list.check_plane(j, k) == cluster_rule_list[l].plane or m_structure_list.check_plane(j, k) == 'ALL':
                            if m_structure_list.phase_name == cluster_rule_list[l].phase:
                                if cluster_rule_list[l].neighbor_arrangement == 'COMB':
                                    if m_structure_list.basis[k].species in cluster_rule_list[l].neighbor_atom_list:
                                        m_structure_list.Cluster_sums[l] += 1
                                if cluster_rule_list[l].neighbor_arrangement == 'PERM':
                                    if m_structure_list.basis[k].species in cluster_rule_list[l].neighbor_atom_list:
                                        if m_structure_list.basis[k].species != m_structure_list.basis[j].species:
                                            m_structure_list.Cluster_sums[l] += 1
            # Calc J sums
            for l in range(len(j_rule_list)):
                if m_structure_list.basis[j].species in j_rule_list[l].home_atom_list:
                    if m_structure_list.distances[j, k] == m_structure_list.mins[j, j_rule_list[l].neighbor_order - 1]:
                        if m_structure_list.check_plane(j, k) == j_rule_list[l].plane or m_structure_list.check_plane(j, k) == 'ALL':
                            if m_structure_list.phase_name == j_rule_list[l].phase:
                                if j_rule_list[l].neighbor_arrangement == 'COMB':
                                    if m_structure_list.basis[k].species in j_rule_list[l].neighbor_atom_list:
                                        m_structure_list.J_sums[l] += m_structure_list.basis[j].spin * m_structure_list.basis[k].spin
                                        # will change in the future determine Mg - Mg only
                                        if m_structure_list.basis[j].species != 2:
                                            if m_structure_list.basis[k].species != 2:
                                                j_count[l]+=1
                                if j_rule_list[l].neighbor_arrangement == 'PERM':
                                    if m_structure_list.basis[k].species in j_rule_list[l].neighbor_atom_list:
                                        if m_structure_list.basis[k].species != m_structure_list.basis[j].species:
                                            m_structure_list.J_sums[l] += m_structure_list.basis[j].spin * m_structure_list.basis[k].spin
                                            if m_structure_list.basis[j].species != 2:
                                                if m_structure_list.basis[k].species != 2:
                                                    j_count[l]+=1

    return j_count

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


def phase_determine(structures,norms,J_rules):
    # find out the corresponding index for each J_rule we want to look at
    for i in range(len(J_rules)):
        J_rule = J_rules[i]
        name = ''.join(J_rule.name.split())
        if name == "J1b0":
            J1b0 = i
        elif name == "J1bU":
            J1bU = i
        elif name == "J2bU":
            J2bU = i
        elif name == "J1c0":
            J1c0 = i
        elif name == "J1cU":
            J1cU = i
        elif name == "J2cU":
            J2cU = i

    for i in range(len(structures)):
        mat = structures[i]
        Jsums = mat.J_sums
        Jcounts = norms[i]
        Jscale = []
        # for each m_structure, get the j_sum scaled by number of Mg-Mg bonds
        for j in range(len(Jcounts)):
            Jcount = Jcounts[j]
            if Jcount == 0:
                Jscale.append(0.0)
            else:
                Jscale.append(Jsums[j] / Jcount)
        # print(Jscale)
        # step 1: check if all 0
        if (all( v >= -0.1 and v <= 0.1 for v in Jscale)):
            mat.mag_phase = "NA"
        else:
            #step 2: if it is austinite
            if mat.phase_name == "aust":
                # check J1b0 first
                phase =  check_aust(Jscale,J1b0)
                if phase != 0:
                    mat.mag_phase = phase
                else:
                    #check J1c0 then
                    phase = check_aust(Jscale,J1c0)
                    if phase != 0:
                        mat.mag_phase = phase
                    else:
                        mat.mag_phase = "sd"
            #step 3: if it is martensite
            elif mat.phase_name == "mart":
                # check J2bU and J1bU first
                phase = check_mart(Jscale,J1bU,J2bU)
                if phase != 0:
                    mat.mag_phase = phase
                else:
                    #check J1cU and J2cU
                    phase = check_mart(Jscale,J1cU,J2cU)
                    if phase!= 0:
                        mat.mag_phase = phase
                    else:
                        mat.mag_phase = "sd"
            else:
                # when it is not aust or mart? can be pm i guess?
                mat.mag_phase = "sd"

def check_aust (Jscale,J_rule):
    determine = Jscale[J_rule]
    if determine <= 1 and determine >= 0.9:
        return "fm"
    elif determine >= -1 and determine <= -0.9:
        return "afm"
    else:
        return 0

def check_mart (Jscale,J1,J2):
    first = Jscale[J1]
    second = Jscale[J2]
    # if one of them is -1, afm
    if (first <= -0.9 and first >= -1) or (second <= -0.9 and second >= -1):
        return "afm"
    # if both are 1, fm
    elif (first >= 0.9 and first <= 1) and (second >= 0.9 and second <= 1):
        return  "fm"
    else:
        return 0

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

def summarize_classification(structures,norms):
    path = 'summary_classification'
    file = open(path, 'w')
    file.write("NAME".ljust(15) + "SCALED_SUMS->\n")

    for i in range(len(structures)):
        mat = structures[i]
        Jsums = mat.J_sums
        Jcounts = norms[i]
        Jscale = []
        for j in range(len(Jcounts)):
            Jcount = Jcounts[j]
            if Jcount == 0:
                Jscale.append(0.0)
            else:
                Jscale.append(Jsums[j]/Jcount)
        file.write(mat.name.ljust(15) + str(Jscale).ljust(17))
        file.write("\n")
    file.close()

