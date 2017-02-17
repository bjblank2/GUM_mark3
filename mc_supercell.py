__author__ = 'brian'
import numpy as np

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
        return self.phase

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

    def get_index(self):
        return self.index

    def get_pos(self):
        return self.pos

    def get_i_pos(self):
        return self.i_pos

    def get_j_pos(self):
        return self.j_pos

    def get_k_pos(self):
        return self.k_pos

    def get_order(self):
        return self.order

    def get_plain(self):
        return self.plain

# mc_supercellObj is the class that holds all the information on the MC simulation as well as the array of atoms.
# Each atom is a mc_siteObj (lattice site) and holds information like spin,position,phase ect...
class mc_supercellObj:
    # The most important function in this class is the initialization function. This is the one I modify to
    # change the starting phase and magnetic structure for each simulation.
    # The lines to modify are followed by #############
    def __init__(self, size, species, composition):
        self.i_length = size[0]
        self.j_length = size[1]
        self.k_length = size[2]
        self.composition = composition
        self.num_sites = size[0]*size[1]*size[2]
        self.supercell = np.empty((self.i_length,self.j_length,self.k_length),dtype=mc_siteObj)
        index = 0
        for i in range(self.i_length):
            for j in range(self.j_length):
                for k in range(self.k_length):
                    spin_rand = np.random.random()
                    phase_rand = np.random.random()
                    if spin_rand <= 1/3: ##########
                        spin = -1 #################
                    elif spin_rand <= 2/3: ########
                        spin = 0 ##################
                    else: #########################
                        spin = 1 ##################
                    # if np.mod(k,2) == 0:
                    #     if np.mod(i+j,2) == 0:
                    #         spin = 1
                    #     else:
                    #         spin = -1
                    # else:
                    #     spin = 0
                    if phase_rand <= 1/3: ##########
                        phase = 1 #################
                    elif phase_rand <= 2/3: ########
                        phase = 1 ##################
                    else: ##########################
                        phase = 1 ##################
                    if np.mod(k,2) == 0: ###########
                        site_species = species[1]
                    else:
                        site_species = species[0]
                    self.supercell[i,j,k] = mc_siteObj(index,(i,j,k),site_species,spin,phase)
                    index += 1
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
        self.find_neighbors()

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
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,1,'OUT'))

                    neighbor_site = self.supercell[i,j,self.apply_bc(k,-1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,1,'OUT'))

                    neighbor_site = self.supercell[self.apply_bc(i,-1,self.i_length),j,self.apply_bc(k,1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,1,'OUT'))

                    neighbor_site = self.supercell[self.apply_bc(i,-1,self.i_length),j,self.apply_bc(k,-1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,1,'OUT'))

                    neighbor_site = self.supercell[i,self.apply_bc(j,-1,self.j_length),self.apply_bc(k,1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,1,'OUT'))

                    neighbor_site = self.supercell[i,self.apply_bc(j,-1,self.j_length),self.apply_bc(k,-1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,1,'OUT'))

                    neighbor_site = self.supercell[self.apply_bc(i,-1,self.i_length),self.apply_bc(j,-1,self.j_length),self.apply_bc(k,1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,1,'OUT'))

                    neighbor_site = self.supercell[self.apply_bc(i,-1,self.i_length),self.apply_bc(j,-1,self.j_length),self.apply_bc(k,-1,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,1,'OUT'))


                    neighbor_site = self.supercell[self.apply_bc(i,1,self.i_length),j,k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,2,'IN'))

                    neighbor_site = self.supercell[self.apply_bc(i,-1,self.i_length),j,k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,2,'IN'))

                    neighbor_site = self.supercell[i,self.apply_bc(j,1,self.j_length),k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,2,'IN'))

                    neighbor_site = self.supercell[i,self.apply_bc(j,-1,self.j_length),k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,2,'IN'))

                    neighbor_site = self.supercell[i,j,self.apply_bc(k,2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,2,'OUT'))

                    neighbor_site = self.supercell[i,j,self.apply_bc(k,-2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,2,'OUT'))


                    neighbor_site = self.supercell[self.apply_bc(i,1,self.i_length),j,self.apply_bc(k,2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,3,'OUT'))

                    neighbor_site = self.supercell[self.apply_bc(i,1,self.i_length),j,self.apply_bc(k,-2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,3,'OUT'))

                    neighbor_site = self.supercell[self.apply_bc(i,-1,self.i_length),j,self.apply_bc(k,2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,3,'OUT'))

                    neighbor_site = self.supercell[self.apply_bc(i,-1,self.i_length),j,self.apply_bc(k,-2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,3,'OUT'))

                    neighbor_site = self.supercell[i,self.apply_bc(j,1,self.j_length),self.apply_bc(k,2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,3,'OUT'))

                    neighbor_site = self.supercell[i,self.apply_bc(j,1,self.j_length),self.apply_bc(k,-2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,3,'OUT'))

                    neighbor_site = self.supercell[i,self.apply_bc(j,-1,self.j_length),self.apply_bc(k,2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,3,'OUT'))

                    neighbor_site = self.supercell[i,self.apply_bc(j,-1,self.j_length),self.apply_bc(k,-2,self.k_length)]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,3,'OUT'))

                    neighbor_site = self.supercell[self.apply_bc(i,1,self.i_length),self.apply_bc(j,1,self.j_length),k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,3,'IN'))

                    neighbor_site = self.supercell[self.apply_bc(i,1,self.i_length),self.apply_bc(j,-1,self.j_length),k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,3,'IN'))

                    neighbor_site = self.supercell[self.apply_bc(i,-1,self.i_length),self.apply_bc(j,1,self.j_length),k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,3,'IN'))

                    neighbor_site = self.supercell[self.apply_bc(i,-1,self.i_length),self.apply_bc(j,-1,self.j_length),k]
                    site_index = neighbor_site.get_index()
                    site_pos = neighbor_site.get_pos()
                    self.supercell[i,j,k].add_neighbor(mc_neighborObj(site_index,site_pos,3,'IN'))

    def get_composition(self):
        return self.composition

    def get_site_index(self,site):
        return self.supercell[site[0],site[1],site[2]].get_index()

    def get_site_pos(self,site):
        return self.supercell[site[0],site[1],site[2]].get_pos()

    def get_site_species(self,site):
        return self.supercell[site[0],site[1],site[2]].get_species()

    def get_site_spin(self,site):
        return self.supercell[site[0],site[1],site[2]].get_spin()

    def get_site_phase(self,site):
        return self.supercell[site[0],site[1],site[2]].get_phase()

    def set_site_species(self,site,species):
        self.supercell[site[0],site[1],site[2]].set_species(species)

    def set_site_spin(self,site,spin):
        self.supercell[site[0],site[1],site[2]].set_spin(spin)

    def set_site_phase(self,site,phase):
        self.supercell[site[0],site[1],site[2]].set_phase(phase)

    def set_neighbor_phase(self,site,neighbor,phase):
        neighbor_site = self.get_neighbor_pos(site,neighbor)
        self.set_site_phase(neighbor_site,phase)

    def get_number_of_neighbors(self,site):
        return len(self.supercell[site[0],site[1],site[2]].neighbors)

    def get_neighbor_pos(self,site,neighbor):
        return self.supercell[site[0],site[1],site[2]].neighbors[neighbor].get_pos()

    def get_neighbor_order(self,site,neighbor):
        return self.supercell[site[0],site[1],site[2]].neighbors[neighbor].get_order()

    def get_neighbor_plain(self,site,neighbor):
        return self.supercell[site[0],site[1],site[2]].neighbors[neighbor].get_plain()

    def get_neighbor_phase(self,site,neighbor):
        neighbor_pos = self.supercell[site[0],site[1],site[2]].neighbors[neighbor].get_pos()
        return self.supercell[neighbor_pos[0],neighbor_pos[1],neighbor_pos[2]].get_phase()

    def get_neighbor_species(self,site,neighbor):
        neighbor_pos = self.supercell[site[0],site[1],site[2]].neighbors[neighbor].get_pos()
        return self.supercell[neighbor_pos[0],neighbor_pos[1],neighbor_pos[2]].get_species()

    def get_neighbor_spin(self,site,neighbor):
        neighbor_pos = self.supercell[site[0],site[1],site[2]].neighbors[neighbor].get_pos()
        return self.supercell[neighbor_pos[0],neighbor_pos[1],neighbor_pos[2]].get_spin()

    def check_site_phase(self,site):
        phase = self.supercell[site[0],site[1],site[2]].get_phase()
        if abs(phase) == 1:
            phase_string = 'mart'
        if phase == 0:
            phase_string = 'aust'
        return phase_string
