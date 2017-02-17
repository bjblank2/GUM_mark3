__author__ = 'chendi'

import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import calc_params
from scipy.optimize import lsq_linear
# with open('Temp_data') as f:
#     a, b,c,d,e,g = [str(x) for x in next(f).split()] # read first line
#     print(a,b,c,d,e,g)
#     array = []
#   #  print([str(x) for x in next(f).split()])
#     for line in f: # read rest of lines
#         array.append([float(x) for x in line.split()])
#
# array = np.array(array)
# print(array)
# temp = array[:,0]
# H_ave = array[:,1]
# mag_ave = array[:,2]
# mag2_ave = array[:,3]
# phase_ave = array[:,4]
# phase2_ave = array[:,5]
#
# fig1 = plt.plot(temp,H_ave,'-')
# plt.xlabel('temp')
# plt.ylabel('H_ave')
# plt.savefig('H_ave vs temp')
#
# plt.show()
# plt.figure()
# plt.plot(temp,mag_ave)
# plt.xlabel('temp')
# plt.ylabel('mag_ave')
# plt.savefig('mag_ave vs temp')
#
# plt.show()
#
# plt.figure()
# plt.plot(temp,mag2_ave)
# plt.xlabel('temp')
# plt.ylabel('mag2_ave')
# plt.savefig('mag2_ave vs temp')
#
# plt.show()
#
# plt.figure()
# plt.plot(temp,phase_ave)
# plt.xlabel('temp')
# plt.ylabel('phase_ave')
# plt.savefig('phase_ave vs temp')
#
# plt.show()
#
# plt.figure()
# plt.plot(temp,phase2_ave)
# plt.xlabel('temp')
# plt.ylabel('phase2_ave')
# plt.savefig('phase2_ave vs temp')

# plt.show()

# fitting part

# should define the function we want to fit here
# data should be different I will leave it here\
# ave_mag = np.array([0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.4975, 0.495, 0.5, 0.4975, 0.485, 0.4575, 0.4525, 0.3575, 0.0675, 0.0525, 0.0375, 0.0525, 0.02, 0.015, 0.0175, 0.025, 0.04, 0.0025, 0.0475, 0.0075, 0.005, 0.02, 0.045, 0.0425, 0.005, 0.0075, 0.0175, 0.0075, 0.025, 0.0025, 0.0125, 0.0025, 0.0175, 0.03, 0.0475])
# kbt = np.array([-5.0,  -2.0, -1.0, -0.95860731484177486, -0.91721462968354983, -0.8758219445253248, -0.83442925936709966, -0.79303657420887463, -0.75164388905064961, -0.71025120389242458, -0.66885851873419944,  -0.58607314841774938, -0.54468046325952435, -0.50328777810129921, -0.46189509294307418, -0.4205024077848491, -0.37910972262662401, -0.33771703746839893, -0.2963243523101739, -0.25493166715194882, -0.2135389819937237, -0.17214629683549859, -0.13075361167727353, -0.089360926519048478, -0.047968241360823373, -0.0065755562025982955, 0.03481712895562681, 0.076209814113851887, 0.11760249927207696, 0.15899518443030203, 0.20038786958852711, 0.24178055474675214, 0.28317323990497723, 0.32456592506320231, 0.36595861022142739, 0.40735129537965253, 0.44874398053787756, 0.49013666569610265, 0.53152935085432773, 0.57292203601255276, 0.6143147211707779, 0.65570740632900293, 0.69710009148722796, 0.7384927766454531, 0.77988546180367813, 0.82127814696190327, 0.86267083212012829, 0.90406351727835343, 0.94545620243657846, 0.98684888759480349, 1.0282415727530287, 2.0282415727530285, 3.0282415727530285, 4.0282415727530285, 5.0282415727530285])
#magwewant = ave_mag[4:len(ave_mag)-4]
#kbtwewant = kbt[4:len(ave_mag)-4]
#magwewant = ave_mag[16:len(ave_mag)-16]
#kbtwewant = kbt[16:len(ave_mag)-16]

def setup_matrix (M_structures):
    matrix = []
    enrg = []
    for i in range (len(M_structures)):
        if M_structures[i].mag_phase!= 'pera':
            if M_structures[i].phase_name!='pm':
                row = M_structures[i].BEG_sums + M_structures[i].Cluster_sums + M_structures[i].J_sums
                matrix.append(row)
                enrg.append(M_structures[i].enrg)
    return(matrix,enrg)

def cal_params (M_structures):
    matrix,enrg = setup_matrix(M_structures)
    params = linalg.lstsq(matrix,enrg)
    return params
def func (matrix,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25)  :
    y = []
    for i in range(len(matrix)):
        row = matrix[i]
        value = row[0]*a0 + row[1]*a1 + row[2]*a2 + row[3]*a3+ row[4]*a4+ row[5]*a5+ row[6]*a6 + row[7]*a7 + row[8]*a8 + row[9]*a9 + row[10]*a10
        value = value +   row[11]*a11 + row[12]*a12 + row[13]*a13+ row[14]*a14+ row[15]*a15+ row[16]*a16 + row[17]*a17 + row[18]*a18 + row[19]*a19 + row[20]*a20
        value = value +row[21]*a1 + row[22]*a22 + row[23]*a23+ row[24]*a24+ row[25]*a25
        y.append(value)
    return(y)

# not working well
def cal_params2(M_structures):
#     initial_guess = np.array([0.001,0.0001,-0.0001,0.005,0.01,-0.005,-0.001,0.02,0.03,0.04,0.05,-0.02,-0.03,-0.04,-0.05,0.1,-0.1,0.075,-0.075])
#     bestguess = 0
     matrix,engr = setup_matrix(M_structures)
#     bestcount = 0
#     for i in range (len(initial_guess)):
#         count = 0
#         p0 = initial_guess[i]
#         print(p0)
# # use computational method to find out best initial guess
# # try to figure out initial guess
#
#         params = curve_fit(func, matrix, engr,(p0+0.5,p0+0.5,p0+0.5,p0+0.5,p0+0.5,p0+0.5,p0+0.5,p0+1,p0+1,p0+1,p0+1,p0+1,p0+1,p0+1,p0+1,p0,p0,p0,p0,p0,p0,p0,p0,p0,p0,p0),maxfev = 100000)
#
#         [a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25] = params[0]
# #params = curve_fit(func, kbt, ave_mag)
#
#
#
#         y = func(matrix,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,)
#         for j in range (len(y)):
#             if (abs((y[j]-engr[j])/engr[j])<0.01):
#                 count+=1
#
#         if count>bestcount:
#             bestguess = i
#
#     print(i)
#     now_guess = initial_guess[i]
#     print(now_guess)
#     step = now_guess/100
#     kth_step = now_guess-50*step
#     bestcount = 0
#     print(kth_step)
#     print(now_guess+50*step)
#     while (abs(kth_step)<=abs(now_guess+50*step)):
#         count = 0
#         p0 = kth_step
#     # initial guess for J can be different
#         params = curve_fit(func, matrix, engr,(p0+0.5,p0+0.5,p0+0.5,p0+0.5,p0+0.5,p0+0.5,p0+0.5,p0+1,p0+1,p0+1,p0+1,p0+1,p0+1,p0+1,p0+1,p0,p0,p0,p0,p0,p0,p0,p0,p0,p0,p0),maxfev = 100000)
#
#         [a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25] = params[0]
#     # params = curve_fit(func, kbt, ave_mag)
#
#
#         y = func(matrix,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25)
#         for j in range (len(y)):
#              if (abs((y[j]-engr[j])/engr[j])<0.01):
#                  count+=1
#
#         if count > bestcount:
#             bestguess = kth_step
#         print(kth_step)
#         kth_step+=step
#     print(bestguess)
#     p0 = initial_guess[bestguess]
#     params = curve_fit(func, matrix, engr,(p0+0.5,p0+0.5,p0+0.5,p0+0.5,p0+0.5,p0+0.5,p0+0.5,p0+1,p0+1,p0+1,p0+1,p0+1,p0+1,p0+1,p0+1,p0,p0,p0,p0,p0,p0,p0,p0,p0,p0,p0),maxfev = 100000)
     p0= 0
     #params_bounds = ([0,0.01715002,10],[0,0.01402785,10],[0,0.01449548,10],[0,0.01244507,10],[0,0.00768251,10],[0,0.00544701,10],[-np.inf,0.01312224,np.inf],[-np.inf,0.00583003,np.inf],[-np.inf,-0.01917217,np.inf],[-np.inf,-0.0066806,np.inf],[-np.inf,0.00427139,np.inf],[-np.inf,0.04631006,np.inf],[-np.inf,0.00552834,np.inf],[-np.inf,-0.00026941,np.inf],[-np.inf,-0.02222189,np.inf],[-np.inf, 0.00165294,np.inf],[-np.inf,-0.00190895,np.inf],[-np.inf,-0.0006273,np.inf],[-np.inf,0.03032933,np.inf],[-np.inf,-0.0266784,np.inf],[-np.inf,0.0006747,np.inf],[-np.inf,-0.02276243,np.inf],[-np.inf,-0.00108391,np.inf],[-np.inf,-0.00956571,np.inf],[-np.inf,0.00917091,np.inf],[-np.inf,0.004585,np.inf])
     params_bounds = ([0,0,0,0,0,0,-50,-50,-50,-50,-50,-50,-50,-50,-50,-50,-50,-50,-50,-50,-50,-50,-50,-50,-50,-50],[10,10,10,10,10,10,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50])
     params = curve_fit(func, matrix, engr,bounds = params_bounds)
     [a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25] = params[0]
     return(params[0])



def cal_params3 (M_structures):
    # least square with constraint
    matrix, engr = setup_matrix(M_structures)
    lb = np.ones(26)*(-np.inf)
    ub = np.zeros(26)
    for i in range (20):
       # lb[i+6] = -np.inf
        ub[i+6] = np.inf
    print(lb,ub)
    params = lsq_linear(matrix,engr,bounds=(lb, ub), lsmr_tol='auto', verbose=1)
    print(params.x)
    return (params.x)
        # intial_guess = [0.01715002,0.01402785,0.01449548,0.01244507,0.00768251,0.00544701,0.01312224,0.00583003,-0.01917217,-0.0066806,0.00427139,0.04631006,0.00552834,-0.00026941,-0.02222189, 0.00165294,-0.00190895,-0.0006273,0.03032933,-0.0266784,0.0006747,-0.02276243,-0.00108391,-0.00956571,0.00917091,0.004585]

