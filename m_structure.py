__author__ = 'brian'
import atom
import numpy as np

def determine_phaseu(LCs,aust_tol):
    a = LCs[0];
    b = LCs[1];
    c = LCs[2];
    find_c = [abs(a - b), abs(a - c), abs(b - c)]
    min_ind = find_c.index(min(find_c))
    if min_ind == 0:
        C = c
        A = a
        B = b
        u = c / ((a+b)/2)
    elif min_ind == 1:
        C = b
        A = a
        B = c
        u = b / ((a+c)/2)
    elif min_ind == 2:
        C = a
        A = c
        B = b
        u = a / ((c+b)/2)
    if u > (1+aust_tol):
        phase = "mart"
    if u <= (1+aust_tol) and u >= (1-aust_tol):
        phase = "aust"
    if u < (1-aust_tol):
        phase = "prem"
    return(phase,u)

def determine_spin():
    return('sd')

class MStructureObj:
    def __init__(self, data, species, num_cluster_rules, num_j_rules, aust_tol):
        data = data.split()
        self.num_cluster_rules = num_cluster_rules
        self.num_j_rules = num_j_rules
        self.species = species
        itter = 0
        self.composition = [0] * (len(species)-1)
        self.num_Atoms = 0
        for i in range(len(species)-1):
            self.num_Atoms += int(data[itter])
            self.composition[i] = int(data[itter])
            itter += 1
        self.name = data[itter]
        self.enrg = float(data[itter + 1])
        a = float(data[itter + 2])
        b = float(data[itter + 3])
        c = float(data[itter + 4])
        self.LCs = [a, b, c]
        self.phase_name,self.u = determine_phaseu(self.LCs,aust_tol)
        self.mag_phase = determine_spin()   # does this ever get used anywhere?
        find_c = [abs(b - c), abs(a - c), abs(a - b)]
        self.Cindex = find_c.index(min(find_c))
        self.C = self.LCs[self.Cindex]
        self.LCs[self.Cindex] = self.LCs[2]
        self.LCs[2] = self.C
        self.weight = 1.0
        self.Cluster_sums = [0.0] * num_cluster_rules
        self.J_sums = [0.0] * num_j_rules
        self.basis = []
        self.distances = np.ones([self.num_Atoms, self.num_Atoms * 27]) * 100
        self.mins = np.ones([self.num_Atoms, 10]) * 100

    def set_atom_properties(self, index, atom_data, spin_style, spin_tol):
        atom_data = atom_data.split()
        mag = float(atom_data[1])
        #print("atom number ",index," , magnetization is ",mag)
        pos = [round(float(atom_data[2]), 5), round(float(atom_data[3]), 5), round(float(atom_data[4]), 5)]
        self.basis.append(atom.AtomObj(index, self.composition, mag, pos, spin_style, spin_tol, self.Cindex))

    def create_supercell(self,spin_style, spin_tol):
        for i in range(-1, 2):
            for j in range(-1, 2):
                for k in range(-1, 2):
                    if abs(i) + abs(j) + abs(k) != 0:
                        for l in range(self.num_Atoms):
                            atom_copy = self.basis[l]
                            mag = atom_copy.mag
                            pos = [atom_copy.a_pos + i, atom_copy.b_pos + j, atom_copy.c_pos + k]
                            self.basis.append(atom.AtomObj(l, self.composition, mag, pos, spin_style, spin_tol))
## wait - should append line above include self.Cindex??

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
        if (abs(self.basis[home_atom_index].c_pos-self.basis[neighbor_atom_index].c_pos))<0.03:
            plane = 'IN'        # has been rotated so that c_pos is the one skew axis
                #            if(neighbor_atom_index<17) :
                #print( 'now j = ',home_atom_index,', k = ',neighbor_atom_index,': these are in plane \n')
        else:
            plane = 'OUT'
        if self.phase_name == 'aust':
            plane = 'ALL'
        return plane
