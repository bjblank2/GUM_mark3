__author__ = 'brian'


class BEGObj:
    def __init__(self):
        self.name = ''
        self.neighbor_order = 0
        self.neighbor_arrangement = ''
        self.home_atom_list = []
        self.neighbor_atom_list = []
        self.phase = ''
        self.plane = ''

    def create_rule(self):
        self.name = input('Input J name:  ')
        self.neighbor_order = int(input('Input neighbor order:  '))
        self.neighbor_arrangement = input('Permutation (PERM) or Combination (COMB) of home/neighbor pairs:  ')
        r_list = input('Input all home atom species (separated by a space):  ')
        self.home_atom_list = r_list.split()
        r_list = input('Input all neighbor atom species (separated by a space):  ')
        self.neighbor_atom_list = r_list.split()
        self.phase = input('Input phase (aust / mart):  ')
        self.plane = input('In-plane (IN) or Out-of-plane (OUT):  ')

    def set_name(self, name):
        self.name = name

    def set_neighbor_order(self, order):
        self.neighbor_order = order

    def set_neighbor_arrangement(self, arrangement):
        self.neighbor_arrangement = arrangement

    def set_home_atom_list(self, ha_list):
        self.home_atom_list = ha_list
        for i in range(len(ha_list)):
            self.home_atom_list[i] = float(self.home_atom_list[i])

    def set_neighbor_atom_list(self, na_list):
        self.neighbor_atom_list = na_list
        for i in range(len(na_list)):
            self.neighbor_atom_list[i] = float(self.neighbor_atom_list[i])

    def set_phase(self, phase):
        self.phase = phase

    def set_plane(self, plane):
        self.plane = plane
__author__ = 'brian'
