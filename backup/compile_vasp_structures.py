_author__ = 'brian'
# This file handles taking the input from vasp and calculating the fit.
# The functions that are important for the actual calculating the fit
# have comments

import numpy as np
import clusters
import js
import m_structure
import os
import matplotlib.pyplot as plt

def import_vasp(root_dir, output_dir,species):
    output = open(output_dir, 'w')
    for subdir, dirs, files in os.walk(root_dir):
        contcar_lines = []
        outcar_lines = []
        flag = 0
        for file in files:
            if 'CONTCAR' in files and 'OUTCAR' in files:
                name = subdir.replace(root_dir,"")
                name = name.replace("/","")
                flag = 1
                if file == "CONTCAR":
                    contcar = open(subdir + '/' + file, 'r')
                    contcar_lines = contcar.readlines()
                    contcar.close()
                if file == "OUTCAR":
                    outcar = open(subdir + '/' + file, 'r')
                    outcar_lines = outcar.readlines()
                    outcar_len = len(outcar_lines)
                    outcar.close()
            else:
                flag = 0
        if len(contcar_lines) > 0 and len(outcar_lines) > 0 and flag == 1:
            # species order is in contcar_line[0]
            order = contcar_lines[5].split()
            #print(order)
            lc = float(contcar_lines[1])                # ARE ALL OF THESE IN ORDER A, B, C???
            a = contcar_lines[2].split()
            a = float(a[0]) * lc
            b = contcar_lines[3].split()
            b = float(b[1]) * lc
            c = contcar_lines[4].split()
            c = float(c[2]) * lc
            #composition = [0] * len(species)
            comp = contcar_lines[6].split()
            composition = [[0,species[-1]]] *len(species)
            for i in range(len(comp)):
                composition[i]=[comp[i],order[i]]

            composition_ordered = sorted(composition, key=lambda d: species.index(d[1]))


            output.write("# ")
            for i in range(len(species)):
                output.write(str(species[i])+ " ")
            output.write('\n')
            total_num = 0
            for i in range(len(composition_ordered)):
                output.write(str(composition_ordered[i][0]) + "\t")
                total_num += int(composition_ordered[i][0])
            mag_list = [[0,species[-1]]] * (total_num)
            for i in range(outcar_len):
                if "TOTEN" in outcar_lines[i]:
                    enrg = outcar_lines[i].split()
                    enrg = float(enrg[4])
                if "magnetization (x)" in outcar_lines[i]:
                    for j in range(total_num):
                        mag = outcar_lines[i + j + 4].split()
                        if j < int(composition[0][0]):
                            kind = order[0]
                        elif j< int(composition[0][0]) + int(composition[1][0]):
                            kind = order[1]
                        else:
                            kind = order[2]
                        mag_list[j] = [mag[4],kind]
            #print(name)
            #print(mag_list)
            mag_list_ordered = sorted(mag_list, key=lambda d: species.index(d[1]))
            enrg = str(enrg)
            a = str(a)
            b = str(b)
            c = str(c)
            output_line = name + "\t" + enrg + "\t" + a + "\t" + b + "\t" + c + "\n"
            output.write(output_line)
            index = 0
            pos_list = []
            for i in range(8, 8 + total_num):
                pos = contcar_lines[i].split()
                if i-8 < int(composition[0][0]):
                    kind = order[0]
                elif i-8 < int(composition[0][0]) + int(composition[1][0]):
                    kind = order[1]
                else:
                    kind = order[2]
                pos_list.append([pos,kind])

            # resort pos according to the order of species
            pos_list_ordered = sorted(pos_list, key=lambda d: species.index(d[1]))
            for i in range(total_num):
                output_line = "\t" + str(i ) + "\t" + str(mag_list_ordered[index][0]) + "\t" + str(pos_list_ordered[i][0][0]) + "\t" + str(
                    pos_list_ordered[i][0][1]) + "\t" + str(pos_list_ordered[i][0][2]) + "\n"
                output.write(output_line)
                index += 1
    output.close()


