__author__ = 'brian'

import numpy as np
cimport numpy as np
from mc_neighborObjCY cimport mc_neighborObj

cdef class mc_siteObj:
    def __init__(self, int index, list pos, int species, int spin, int phase):
        self.index = index
        self.pos = pos
        self.species = species
        self.spin = spin
        self.phase = phase
        self.neighbors = []

    cdef int get_index(self):
        return self.index

    cdef list get_pos(self):
        return self.pos

    cdef int get_species(self):
        return self.species

    cdef int get_spin(self):
        return self.spin

    cdef int get_phase(self):
        return self.phase

    cdef void set_species(self,int new_value):
        self.species = new_value

    cdef void set_spin(self,int new_value):
        self.spin = new_value

    cdef void set_phase(self, int new_value):
        self.phase = new_value

    cdef void add_neighbor(self, mc_neighborObj neighbor_obj):
        self.neighbors.append(neighbor_obj)