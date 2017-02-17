__author__ = 'brian'

class AtomObj:
    def __init__(self, index, species_list, mag, pos, c_index=None):
        self.atom_index = index
        self.species = 0
        self.set_species(index, species_list)
        self.mag = mag
        self.spin = 0
        self.set_spin(mag)
        self.set_pos(pos[0], pos[1], pos[2])
        if c_index is not None:
            self.rotate(c_index)

    def set_spin(self, mag):
        if self.species == 0:
            if abs(round(mag, 2)) >= .1:
                self.spin = round(abs(mag) / mag, 5)
            else:
                self.spin = 0
        if self.species == 1:
            #if abs(round(mag, 2)) >= 2:
            if abs(mag) >= 2:
                self.spin = round(abs(mag) / mag, 5)
            else:
                self.spin = 0
        if self.species == 2:
            self.spin = 0

    def set_pos(self, a_pos, b_pos, c_pos):
        if a_pos >= .2 and a_pos <= .3:
            a_pos = .25
        if a_pos >= .7 and a_pos <= .8:
            a_pos = .75
        if b_pos >= .2 and b_pos <= .3:
            b_pos = .25
        if b_pos >= .7 and b_pos <= .8:
            b_pos = .75
        if c_pos >= .2 and c_pos <= .3:
            c_pos = .25
        if c_pos >= .7 and c_pos <= .8:
            c_pos = .75
        self.a_pos = a_pos
        self.b_pos = b_pos
        self.c_pos = c_pos

    def rotate(self, cindex):
        a_pos = self.a_pos
        b_pos = self.b_pos
        c_pos = self.c_pos
        pos = [a_pos, b_pos, c_pos]
        C = pos[cindex]
        pos[cindex] = pos[2]
        pos[2] = C
        self.a_pos = pos[0]
        self.b_pos = pos[1]
        self.c_pos = pos[2]

    def set_species(self, index, species_list):
        species_index = 0
        species_index_old = 0
        for i in range(len(species_list)):
            species_index += species_list[i]
            if index < species_index and index >= species_index_old:
                self.species = i
            species_index_old = species_index
