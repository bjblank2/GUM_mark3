__author__ = 'brian'
import atom
import numpy as np


class MStructureObj:
    def __init__(self, data, num_species,num_beg_rules, num_cluster_rules, num_j_rules):
        data = data.split()
        itter = 0
        self.composition = [0] * num_species
        self.num_Atoms = 0
        for i in range(num_species):
            self.num_Atoms += int(data[itter + 1])
            self.composition[i] = int(data[itter + 1])
            itter += 1
        self.phase_name = data[itter + 1]
        self.mag_phase = data[itter + 2]
        self.name = data[itter + 3]
        self.u = float(data[itter + 4])
        self.enrg = float(data[itter + 5])
        a = float(data[itter + 6])
        b = float(data[itter + 7])
        c = float(data[itter + 8])
        self.LCs = [a, b, c]
        find_c = [abs(b - c), abs(a - c), abs(a - b)]
        self.Cindex = find_c.index(min(find_c))
        self.C = self.LCs[self.Cindex]
        self.LCs[self.Cindex] = self.LCs[2]
        self.LCs[2] = self.C
        self.weight = 1.0
        self.BEG_sums = [0.0] * num_beg_rules
        self.Cluster_sums = [0.0] * num_cluster_rules
        self.J_sums = [0.0] * num_j_rules
        self.basis = []
        self.distances = np.ones([self.num_Atoms, self.num_Atoms * 27]) * 100
        self.mins = np.ones([self.num_Atoms, 10]) * 100

    def set_atom_properties(self, index, atom_data):
        atom_data = atom_data.split()
        mag = float(atom_data[1])
        pos = [round(float(atom_data[2]), 5), round(float(atom_data[3]), 5), round(float(atom_data[4]), 5)]
        self.basis.append(atom.AtomObj(index, self.composition, mag, pos, self.Cindex))

    def create_super_cell(self):
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    if abs(i) + abs(j) + abs(k) != 0:
                        for l in range(self.num_Atoms):
                            atom_copy = self.basis[l]
                            mag = atom_copy.mag
                            pos = [atom_copy.a_pos + i, atom_copy.b_pos + j, atom_copy.c_pos + k]
                            self.basis.append(atom.AtomObj(l, self.composition, mag, pos))

    def calculate_distances(self):
        for i in range(self.num_Atoms):
            for j in range(len(self.basis)):
                dist_a = round(self.basis[i].a_pos - self.basis[j].a_pos, 5)
                dist_b = round(self.basis[i].b_pos - self.basis[j].b_pos, 5)
                dist_c = round(self.basis[i].c_pos - self.basis[j].c_pos, 5)
                self.distances[i, j] = round((dist_a ** 2 + dist_b ** 2 + dist_c ** 2) ** (0.5), 5)
                if self.distances[i, j] == 0:
                    self.distances[i, j] = 100

    def calculate_minimums(self):
        dists = [0] * len(self.basis)
        for i in range(self.num_Atoms):
            old_min = 0.0
            for j in range(len(self.basis)):
                dists[j] = self.distances[i, j]
            for j in range(10):
                self.mins[i, j] = self.next_min(dists, old_min)
                old_min = self.mins[i, j]

    def next_min(self, dists, old_min):
        for i in range(len(dists)):
            if round(dists[i], 5) <= old_min:
                dists[i] = 100
        new_min = np.min(dists)
        return new_min

    def check_plane(self, home_atom_index, neighbor_atom_index):
        if round(self.basis[home_atom_index].c_pos, 5) == round(self.basis[neighbor_atom_index].c_pos, 5):
            plane = 'IN'
        else:
            plane = 'OUT'
        if self.phase_name == 'aust':
            plane = 'ALL'
        return plane
