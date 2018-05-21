__author__ = 'brian'
import numpy as np
cimport numpy as np
from mc_siteObjCY cimport mc_siteObj
from mc_neighborObjCY cimport mc_neighborObj

cdef class mc_supercellObj:

    cdef int i_length
    cdef int j_length
    cdef int k_length
    cdef list composition
    cdef str phase_init
    cdef str spin_init
    cdef str species_init
    cdef int num_sites
    cdef np.ndarray supercell
    cdef mc_siteObj new_site

    cdef int apply_bc(self,int i,int inc,int limit)
    cdef void find_neighbors(self)
    cdef list get_composition(self)
    cdef int get_site_index(self,list site)
    cdef list get_site_pos(self,list site)
    cdef int get_site_species(self,list site)
    cdef int get_site_spin(self,list site)
    cdef int get_site_phase(self,list site)
    cdef void set_site_species(self,list site,int species)
    cdef void set_site_spin(self,list site,int spin)
    cdef void set_site_phase(self,list site,int phase)
    cdef void set_neighbor_phase(self,list site,int neighbor,int phase)
    cdef int get_number_of_neighbors(self,list site)
    cdef list get_neighbor_pos(self,list site,int neighbor)
    cdef int get_neighbor_order(self,list site,int neighbor)
    cdef str get_neighbor_plain(self,list site,int neighbor)
    cdef int get_neighbor_phase(self,list site,int neighbor)
    cdef int get_neighbor_species(self,list site,int neighbor)
    cdef int get_neighbor_spin(self,list site,int neighbor)
    cdef str check_site_phase(self,list site)