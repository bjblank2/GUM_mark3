_author__ = 'brian'
# This file handles taking the input from vasp and calculating the fit.
# The functions that are important for the actual calculating the fit
# have comments

import numpy as np
cimport numpy as np
import clusters
import js
import m_structure
import os
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from sklearn.linear_model import Ridge
from sklearn import linear_model


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





