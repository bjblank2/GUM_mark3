__author__ = 'brian'
import numpy as np
cimport numpy as np

cdef class mc_neighborObj:
    def __init__(self,index,pos,neighbor_order,plain):
        self.index = index
        self.pos = pos
        self.i_pos = pos[0]
        self.j_pos = pos[1]
        self.k_pos = pos[2]
        self.order = neighbor_order
        self.plain = plain

    cdef int get_index(self):
        return self.index

    cdef list get_pos(self):
        return self.pos

    cdef float get_i_pos(self):
        return self.i_pos

    cdef float get_j_pos(self):
        return self.j_pos

    cdef float get_k_pos(self):
        return self.k_pos

    cdef int get_order(self):
        return self.order

    cdef str get_plain(self):
        return self.plain