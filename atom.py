__author__ = 'brian'
import numpy as np

class AtomObj:                      # does this need mag in here?  The object should only contain what is           needed for fitting and Monte Carlo from here on out.
    def __init__(self, index, species_list, mag, pos, spin_style, spin_tol, c_index=None):
        self.atom_index = index
        self.species = 0
        self.set_species(index, species_list)
        self.mag = mag
        self.spin = 0
        self.set_spin(mag,spin_style, spin_tol)
        self.set_pos(pos[0], pos[1], pos[2])
        if c_index is not None:
            self.rotate(c_index)

    def set_spin(self, mag, spin_style, spin_tol):
            if self.species == 2:                   # again somehow assuming that 2 = In
                self.spin = 0
            elif round(mag,2) == 0:
                self.spin = 0
            else:
                if spin_style[self.species] == 'factor':
                    self.spin = np.floor(abs(mag)/spin_tol[self.species])*round(abs(mag) / mag, 5)
                else:
                    if abs(mag) > spin_tol[self.species]:
                        self.spin = round(abs(mag) / mag, 5)
                    else:
                        self.spin = 0

    def set_pos(self, a_pos, b_pos, c_pos):        # again very specific assumptions about positions
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
