_author__ = 'brian'
# This file handles taking the input from vasp and calculating the fit.
# The functions that are important for the actual calculating the fit
# have comments

import numpy as np
#import BEG
import clusters
import js
import m_structure
import os
import matplotlib.pyplot as plt
#from scipy.optimize import least_squares
#from sklearn.linear_model import Ridge
#from sklearn import linear_model

def import_vasp(root_dir, output_dir):
    output = open(output_dir, 'w')
    for subdir, dirs, files in os.walk(root_dir):
        contcar_lines = []
        outcar_lines = []
        flag = 0
        for file in files:
            if 'CONTCAR' in files and 'OUTCAR' in files:
                name = subdir.strip(root_dir)
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
            lc = float(contcar_lines[1])
            a = contcar_lines[2].split()
            a = float(a[0]) * lc
            b = contcar_lines[3].split()
            b = float(b[1]) * lc
            c = contcar_lines[4].split()
            c = float(c[2]) * lc
            composition = [0] * number_of_species
            comp = contcar_lines[6].split()
            for i in range(len(comp)):
                composition[i] = comp[i]
            output.write("#\t")
            total_num = 0
            for i in range(len(composition)):
                output.write(str(composition[i]) + "\t")
                total_num += int(composition[i])
            mag_list = [0] * (total_num)
            for i in range(outcar_len):
                if "TOTEN" in outcar_lines[i]:
                    enrg = outcar_lines[i].split()
                    enrg = float(enrg[4])
                if "magnetization (x)" in outcar_lines[i]:
                    for j in range(total_num):
                        mag = outcar_lines[i + j + 4].split()
                        mag_list[j] = mag[4]
            u = str(u)
            enrg = str(enrg)
            a = str(a)
            b = str(b)
            c = str(c)
            output_line = name + "\t" + u + "\t" + enrg + "\t" + a + "\t" + b + "\t" + c + "\n"
            output.write(output_line)
            index = 0
            for i in range(8, 8 + total_num):
                pos = contcar_lines[i].split()
                output_line = "\t" + str(i - 7) + "\t" + str(mag_list[index]) + "\t" + str(pos[0]) + "\t" + str(
                    pos[1]) + "\t" + str(pos[2]) + "\n"
                output.write(output_line)
                index += 1
    output.close()


