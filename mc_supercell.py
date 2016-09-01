__author__ = 'brian'
import numpy as np
from copy import deepcopy

class mc_siteObj:
    def __init__(self, index, pos, species, spin, phase):
        self.index = index
        self.pos = pos
        self.species = species
        self.spin = spin
        self.phase = phase
        self.neighbors = []

    def get_index(self):
        return self.index

    def get_pos(self):
        return self.pos

    def get_species(self):
        return self.species

    def get_spin(self):
        return self.spin

    def get_phase(self):
        self.phase

    def set_species(self,new_value):
        self.species = new_value

    def set_spin(self,new_value):
        self.spin = new_value

    def set_phase(self,new_value):
        self.phase = new_value

    def add_neighbor(self,neighbor_obj):
        self.neighbors.append(neighbor_obj)



class mc_neighborObj:
    def __init__(self,index,pos,neighbor_order,plain):
        self.index = index
        self.pos = pos
        self.i_pos = pos[0]
        self.j_pos = pos[1]
        self.k_pos = pos[2]
        self.order = neighbor_order
        self.plain = plain


class mc_supercellObj:
    def __init__(self, size, species, composition):
        self.i_length = size[0]
        self.j_length = size[1]
        self.k_length = size[2]
        self.supercell = np.empty((self.i_length,self.j_length,self.k_length),dtype=mc_siteObj)
        self.supercell_list = []
        index = 0
        for i in range(self.i_length):
            for j in range(self.j_length):
                for k in range(self.k_length):
                    spin_rand = np.random.random()
                    phase_rand = np.random.random()
                    if spin_rand <= 1/3:
                        spin = -1
                    elif spin_rand <= 2/3:
                        spin = 0
                    else:
                        spin = 1
                    if phase_rand <= 1/3:
                        phase = -1
                    elif phase_rand <= 2/3:
                        phase = 0
                    else:
                        phase = 1
                    if np.mod(k,2) == 0:
                        site_species = 1
                    else:
                        site_species = 0
                    self.supercell[i,j,k] = mc_siteObj(index,(i,j,k),site_species,spin,phase)
                    self.supercell_list.append(self.supercell[i,j,k])
                    index += 1

    def apply_bc(self,i,inc,limit):
        if i + inc >= limit:
            new_i = i+inc-limit
        elif i + inc < 0:
            new_i = i + inc + limit
        else:
            new_i = i+inc
        return new_i

    def find_neighbors(self):
        for i in range(self.i_length):
            for j in range(self.j_length):
                for k in range(self.k_length):

                    neighbor_site = self.supercell[i,j,self.apply_bc(k,1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nn','OUT'))

                    neighbor_site = self.supercell[i,j,self.apply_bc(k,-1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nn','OUT'))

                    neighbor_site = self.supercell[i,j,self.apply_bc(k,-1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nn','OUT'))

                    neighbor_site = self.supercell[i,self.apply_bc(j,-1,self.j_length),self.apply_bc(k,1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nn','OUT'))

                    neighbor_site = self.supercell[i,self.apply_bc(j,-1,self.j_length),self.apply_bc(k,-1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nn','OUT'))

                    neighbor_site = self.supercell[self.apply_bc(i,-1,self.i_length),self.apply_bc(j,-1,self.j_length),self.apply_bc(k,1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nn','OUT'))

                    neighbor_site = self.supercell[self.apply_bc(i,-1,self.i_length),self.apply_bc(j,-1,self.j_length),self.apply_bc(k,-1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nn','OUT'))

                    neighbor_site = self.supercell[self.apply_bc(i,-1,self.i_length),j,self.apply_bc(k,1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nn','OUT'))

                    neighbor_site = self.supercell[self.apply_bc(i,-1,self.i_length),j,self.apply_bc(k,-1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nn','OUT'))


                    neighbor_site = self.supercell[self.apply_bc(i,1,self.i_length),j,k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nnn','IN'))

                    neighbor_site = self.supercell[self.apply_bc(i,-1,self.i_length),j,k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nnn','IN'))

                    neighbor_site = self.supercell[i,self.apply_bc(j,1,self.j_length),k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nnn','IN'))

                    neighbor_site = self.supercell[i,self.apply_bc(j,-1,self.j_length),k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nnn','IN'))

                    neighbor_site = self.supercell[i,j,self.apply_bc(k,2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nnn','OUT'))

                    neighbor_site = self.supercell[i,j,self.apply_bc(k,-2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nnn','OUT'))


                    neighbor_site = self.supercell[self.apply_bc(i,1,self.i_length),j,self.apply_bc(k,2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nnnn','OUT'))

                    neighbor_site = self.supercell[self.apply_bc(i,1,self.i_length),j,self.apply_bc(k,-2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nnnn','OUT'))

                    neighbor_site = self.supercell[self.apply_bc(i,-1,self.i_length),j,self.apply_bc(k,2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nnnn','OUT'))

                    neighbor_site = self.supercell[self.apply_bc(i,-1,self.i_length),j,self.apply_bc(k,-2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nnnn','OUT'))

                    neighbor_site = self.supercell[i,self.apply_bc(j,1,self.j_length),self.apply_bc(k,2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nnnn','OUT'))

                    neighbor_site = self.supercell[i,self.apply_bc(j,1,self.j_length),self.apply_bc(k,-2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nnnn','OUT'))

                    neighbor_site = self.supercell[i,self.apply_bc(j,-1,self.j_length),self.apply_bc(k,2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nnnn','OUT'))

                    neighbor_site = self.supercell[i,self.apply_bc(j,-1,self.j_length),self.apply_bc(k,-2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nnnn','OUT'))

                    neighbor_site = self.supercell[self.apply_bc(i,1,self.i_length),self.apply_bc(j,1,self.j_length),k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nnnn','IN'))

                    neighbor_site = self.supercell[self.apply_bc(i,1,self.i_length),self.apply_bc(j,-1,self.j_length),k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nnnn','IN'))

                    neighbor_site = self.supercell[self.apply_bc(i,-1,self.i_length),self.apply_bc(j,1,self.j_length),k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nnnn','IN'))

                    neighbor_site = self.supercell[self.apply_bc(i,-1,self.i_length),self.apply_bc(j,-1,self.j_length),k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,'nnnn','IN'))