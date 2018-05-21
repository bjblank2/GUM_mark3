__author__ = 'brian'
import numpy as np
cimport numpy as np
from mc_siteObjCY cimport mc_siteObj
from mc_neighborObjCY cimport mc_neighborObj

# mc_supercellObj is the class that holds all the information on the MC simulation as well as the array of atoms.
# Each atom is a mc_siteObj (lattice site) and holds information like spin,position,phase ect...

cdef class mc_supercellObj:
    # The most important function in this class is the initialization function. This is the one I modify to
    # change the starting phase and magnetic structure for each simulation.
    # The lines to modify are followed by #############
    def __init__(self,list size, list species, list composition,str phase_init,str spin_init,str species_init):
        self.i_length = size[0]
        self.j_length = size[1]
        self.k_length = size[2]
        self.composition = composition
        self.num_sites = size[0]*size[1]*size[2]
        self.supercell = np.empty((self.i_length,self.j_length,self.k_length),dtype=mc_siteObj)
        cdef int index = 0
        cdef int i
        cdef int j
        cdef int k
        cdef int spin_rand
        cdef int spin
        cdef int phase
        cdef int site_species
        cdef int species_count
        cdef list rand_index_list

        for i in range(self.i_length):
            for j in range(self.j_length):
                for k in range(self.k_length):
                    spin_rand = np.random.random()
                    phase_rand = np.random.random()
                    if spin_init == 'FM':
                        spin = 1
                    elif spin_init == 'rand':
                        if spin_rand <= 1/3: ##########
                            spin = -1 #################
                        elif spin_rand <= 2/3: ########
                            spin = 0 ##################
                        else: #########################
                            spin = 1 ##################
                    elif spin_init =='AFM':
                        if np.mod(k,2) == 0:
                            if np.mod(i+j,2) == 0:
                                spin = 1
                            else:
                                spin = -1
                        else:
                            spin = 0
                    else: spin = 1
                    if phase_init == 'aust':
                        phase = 0
                    elif phase_init =='mart':
                        phase = 1
                    elif phase_init == 'rand':
                        if phase_rand <= 1/3: ##########
                            phase = 0 #################
                        elif phase_rand <= 2/3: ########
                            phase = 1 ##################
                        else: ##########################
                            phase = -1 ##################
                    else: phase = 1
                    if np.mod(k,2) == 0:
                        site_species = species[1]
                    else:
                        site_species = species[0]
                    self.supercell[i,j,k] = <mc_siteObj>mc_siteObj(index,[i,j,k],site_species,spin,phase)
                    index += 1
        if composition[2] != 0:
            if composition[1]/composition[2] == 1:
                for i in range(self.i_length):
                    for j in range(self.j_length):
                        for k in range(self.k_length):
                            if np.mod(k,2) == 0:
                                if np.mod(i+j+k*.5,2) == 0:
                                    self.set_site_species([i,j,k],species[2])
                                else: self.set_site_species([i,j,k],species[1])
                            else:
                                self.set_site_species([i,j,k],species[0])
            else:
                if species_init == "rand":
                    species_count = 0
                    rand_index_list = []
                    while species_count < composition[2]:
                        species_not_0 = False
                        while species_not_0 == False:
                            rand_index = [np.random.randint(0,size[0]),np.random.randint(0,size[1]),np.random.randint(0,size[2])]
                            if self.supercell[rand_index[0],rand_index[1],rand_index[2]].species != species[0]:
                                species_not_0 = True
                                if rand_index not in rand_index_list:
                                    self.supercell[rand_index[0],rand_index[1],rand_index[2]].species = species[2]
                                    self.supercell[rand_index[0],rand_index[1],rand_index[2]].spin = 0
                                    rand_index_list.append(rand_index)
                                    species_count += 1
                if species_init == "ordered":
                    if (composition[0]+composition[1]+composition[2])/composition[2] == 8:
                        for i in range(self.i_length):
                            for j in range(self.j_length):
                                for k in range(self.k_length):
                                    if np.mod(k,2) == 0:
                                        if np.mod(i*.5+(j+1)*.5+k*.5,2) == 0:
                                            self.set_site_species([i,j,k],species[2])
                                        else: self.set_site_species([i,j,k],species[1])
                                    else: self.set_site_species([i,j,k],species[0])
        self.find_neighbors()

    cdef int apply_bc(self,int i, int inc, int limit):
        cdef int new_i
        if i + inc >= limit:
            new_i = i+inc-limit
        elif i + inc < 0:
            new_i = i + inc + limit
        else:
            new_i = i+inc
        return new_i

    cdef void find_neighbors(self):

        cdef int i
        cdef int j
        cdef int k
        cdef mc_siteObj neighbor_site
        cdef mc_siteObj home_site
        cdef mc_neighborObj neighbor_Obj

        for i in range(self.i_length):
            for j in range(self.j_length):
                for k in range(self.k_length):

                    home_site = <mc_siteObj>self.supercell[i,j,k]

                    neighbor_site = <mc_siteObj>self.supercell[i,j,self.apply_bc(k,1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,1,'OUT')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[i,j,self.apply_bc(k,-1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,1,'OUT')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[self.apply_bc(i,-1,self.i_length),j,self.apply_bc(k,1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,1,'OUT')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[self.apply_bc(i,-1,self.i_length),j,self.apply_bc(k,-1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,1,'OUT')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[i,self.apply_bc(j,-1,self.j_length),self.apply_bc(k,1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,1,'OUT')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[i,self.apply_bc(j,-1,self.j_length),self.apply_bc(k,-1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,1,'OUT')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[self.apply_bc(i,-1,self.i_length),self.apply_bc(j,-1,self.j_length),self.apply_bc(k,1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,1,'OUT')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[self.apply_bc(i,-1,self.i_length),self.apply_bc(j,-1,self.j_length),self.apply_bc(k,-1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,1,'OUT')
                    home_site.add_neighbor(neighbor_obj)


                    neighbor_site = <mc_siteObj>self.supercell[self.apply_bc(i,1,self.i_length),j,k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,2,'IN')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[self.apply_bc(i,-1,self.i_length),j,k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,2,'IN')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[i,self.apply_bc(j,1,self.j_length),k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,2,'IN')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[i,self.apply_bc(j,-1,self.j_length),k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,2,'IN')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[i,j,self.apply_bc(k,2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,2,'OUT')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[i,j,self.apply_bc(k,-2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,2,'OUT')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[self.apply_bc(i,1,self.i_length),j,self.apply_bc(k,2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,3,'OUT')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[self.apply_bc(i,1,self.i_length),j,self.apply_bc(k,-2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,3,'OUT')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[self.apply_bc(i,-1,self.i_length),j,self.apply_bc(k,2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,3,'OUT')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[self.apply_bc(i,-1,self.i_length),j,self.apply_bc(k,-2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,3,'OUT')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[i,self.apply_bc(j,1,self.j_length),self.apply_bc(k,2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,3,'OUT')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[i,self.apply_bc(j,1,self.j_length),self.apply_bc(k,-2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,3,'OUT')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[i,self.apply_bc(j,-1,self.j_length),self.apply_bc(k,2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,3,'OUT')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[i,self.apply_bc(j,-1,self.j_length),self.apply_bc(k,-2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,3,'OUT')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[self.apply_bc(i,1,self.i_length),self.apply_bc(j,1,self.j_length),k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,3,'IN')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[self.apply_bc(i,1,self.i_length),self.apply_bc(j,-1,self.j_length),k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,3,'IN')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[self.apply_bc(i,-1,self.i_length),self.apply_bc(j,1,self.j_length),k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,3,'IN')
                    home_site.add_neighbor(neighbor_obj)

                    neighbor_site = <mc_siteObj>self.supercell[self.apply_bc(i,-1,self.i_length),self.apply_bc(j,-1,self.j_length),k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    neighbor_obj = <mc_neighborObj>mc_neighborObj(site_index,site_pos,3,'IN')
                    home_site.add_neighbor(neighbor_obj)


    cdef list get_composition(self):
        return <mc_siteObj>self.composition

    cdef int get_site_index(self,list site):
        return <mc_siteObj>self.supercell[site[0],site[1],site[2]].get_index()

    cdef list get_site_pos(self,list site):
        return <mc_siteObj>self.supercell[site[0],site[1],site[2]].get_pos()

    cdef int get_site_species(self, list site):
        return <mc_siteObj>self.supercell[site[0],site[1],site[2]].get_species()

    cdef int get_site_spin(self, list site):
        return <mc_siteObj>self.supercell[site[0],site[1],site[2]].get_spin()

    cdef int get_site_phase(self, list site):
        return <mc_siteObj>self.supercell[site[0],site[1],site[2]].get_phase()

    cdef void set_site_species(self, list site, int species):
        <mc_siteObj>self.supercell[site[0],site[1],site[2]].set_species(species)

    cdef void set_site_spin(self, list site, int spin):
        <mc_siteObj>self.supercell[site[0],site[1],site[2]].set_spin(spin)

    cdef void set_site_phase(self, list site, int phase):
        <mc_siteObj>self.supercell[site[0],site[1],site[2]].set_phase(phase)

    cdef void set_neighbor_phase(self, list site, int neighbor, int phase):
        cdef list neighbor_site = <mc_siteObj>self.get_neighbor_pos(site,neighbor)
        <mc_siteObj>self.set_site_phase(neighbor_site,phase)

    cdef int get_number_of_neighbors(self, list site):
        return len(<mc_siteObj>self.supercell[site[0],site[1],site[2]].neighbors)

    cdef list get_neighbor_pos(self, list site, int neighbor):
        return <mc_siteObj>self.supercell[site[0],site[1],site[2]].neighbors[neighbor].get_pos()

    cdef int get_neighbor_order(self, list site, int neighbor):
        return <mc_siteObj>self.supercell[site[0],site[1],site[2]].neighbors[neighbor].get_order()

    cdef str get_neighbor_plain(self, list site, int neighbor):
        return <mc_siteObj>self.supercell[site[0],site[1],site[2]].neighbors[neighbor].get_plain()

    cdef int get_neighbor_phase(self, list site, int neighbor):
        cdef list neighbor_pos = <mc_siteObj>self.supercell[site[0],site[1],site[2]].neighbors[neighbor].get_pos()
        return <mc_siteObj>self.supercell[neighbor_pos[0],neighbor_pos[1],neighbor_pos[2]].get_phase()

    cdef int get_neighbor_species(self,list site, int neighbor):
        cdef list neighbor_pos = <mc_siteObj>self.supercell[site[0],site[1],site[2]].neighbors[neighbor].get_pos()
        return <mc_siteObj>self.supercell[neighbor_pos[0],neighbor_pos[1],neighbor_pos[2]].get_species()

    cdef int get_neighbor_spin(self, list site, int neighbor):
        cdef list neighbor_pos = <mc_siteObj>self.supercell[site[0],site[1],site[2]].neighbors[neighbor].get_pos()
        return <mc_siteObj>self.supercell[neighbor_pos[0],neighbor_pos[1],neighbor_pos[2]].get_spin()

    cdef str check_site_phase(self, list site):
        phase = <mc_siteObj>self.supercell[site[0],site[1],site[2]].get_phase()
        if abs(phase) == 1:
            phase_string = 'mart'
        if phase == 0:
            phase_string = 'aust'
        return phase_string
