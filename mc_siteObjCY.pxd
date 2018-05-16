__author__ = 'brian'

from mc_neighborObjCY cimport mc_neighborObj

cdef class mc_siteObj:

    cdef int index
    cdef list pos
    cdef int species
    cdef int spin
    cdef int phase
    cdef list neighbors

    cdef int get_index(self)
    cdef list get_pos(self)
    cdef int get_species(self)
    cdef int get_spin(self)
    cdef int get_phase(self)
    cdef void set_species(self,int new_value)
    cdef void set_spin(self,int new_value)
    cdef void set_phase(self,int new_value)
    cdef void add_neighbor(self,mc_neighborObj neighbor_obj)