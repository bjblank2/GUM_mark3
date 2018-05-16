__author__ = 'brian'

cdef class mc_neighborObj:

    cdef int index
    cdef list pos
    cdef float i_pos
    cdef float j_pos
    cdef float k_pos
    cdef int order
    cdef str plain

    cdef int get_index(self)
    cdef list get_pos(self)
    cdef float get_i_pos(self)
    cdef float get_j_pos(self)
    cdef float get_k_pos(self)
    cdef int get_order(self)
    cdef str get_plain(self)