__author__ = 'chendi'

import numpy as np
import matplotlib.pyplot as plt


with open('Temp_data') as f:
    a, b,c,d,e,g = [str(x) for x in next(f).split()] # read first line
    print(a,b,c,d,e,g)
    array = []
  #  print([str(x) for x in next(f).split()])
    for line in f: # read rest of lines
        array.append([float(x) for x in line.split()])

array = np.array(array)
print(array)
temp = array[:,0]
H_ave = array[:,1]
mag_ave = array[:,2]
mag2_ave = array[:,3]
phase_ave = array[:,4]
phase2_ave = array[:,5]

fig1 = plt.plot(temp,H_ave,'-')
plt.xlabel('temp')
plt.ylabel('H_ave')
plt.savefig('H_ave vs temp')

plt.show()
plt.figure()
plt.plot(temp,mag_ave)
plt.xlabel('temp')
plt.ylabel('mag_ave')
plt.savefig('mag_ave vs temp')

plt.show()

plt.figure()
plt.plot(temp,mag2_ave)
plt.xlabel('temp')
plt.ylabel('mag2_ave')
plt.savefig('mag2_ave vs temp')

plt.show()

plt.figure()
plt.plot(temp,phase_ave)
plt.xlabel('temp')
plt.ylabel('phase_ave')
plt.savefig('phase_ave vs temp')

plt.show()

plt.figure()
plt.plot(temp,phase2_ave)
plt.xlabel('temp')
plt.ylabel('phase2_ave')
plt.savefig('phase2_ave vs temp')

plt.show()

# fitting part

# should define the function we want to fit here
# data should be different I will leave it here\
ave_mag = np.array([0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.4975, 0.495, 0.5, 0.4975, 0.485, 0.4575, 0.4525, 0.3575, 0.0675, 0.0525, 0.0375, 0.0525, 0.02, 0.015, 0.0175, 0.025, 0.04, 0.0025, 0.0475, 0.0075, 0.005, 0.02, 0.045, 0.0425, 0.005, 0.0075, 0.0175, 0.0075, 0.025, 0.0025, 0.0125, 0.0025, 0.0175, 0.03, 0.0475])
kbt = np.array([-5.0,  -2.0, -1.0, -0.95860731484177486, -0.91721462968354983, -0.8758219445253248, -0.83442925936709966, -0.79303657420887463, -0.75164388905064961, -0.71025120389242458, -0.66885851873419944,  -0.58607314841774938, -0.54468046325952435, -0.50328777810129921, -0.46189509294307418, -0.4205024077848491, -0.37910972262662401, -0.33771703746839893, -0.2963243523101739, -0.25493166715194882, -0.2135389819937237, -0.17214629683549859, -0.13075361167727353, -0.089360926519048478, -0.047968241360823373, -0.0065755562025982955, 0.03481712895562681, 0.076209814113851887, 0.11760249927207696, 0.15899518443030203, 0.20038786958852711, 0.24178055474675214, 0.28317323990497723, 0.32456592506320231, 0.36595861022142739, 0.40735129537965253, 0.44874398053787756, 0.49013666569610265, 0.53152935085432773, 0.57292203601255276, 0.6143147211707779, 0.65570740632900293, 0.69710009148722796, 0.7384927766454531, 0.77988546180367813, 0.82127814696190327, 0.86267083212012829, 0.90406351727835343, 0.94545620243657846, 0.98684888759480349, 1.0282415727530287, 2.0282415727530285, 3.0282415727530285, 4.0282415727530285, 5.0282415727530285])
#magwewant = ave_mag[4:len(ave_mag)-4]
#kbtwewant = kbt[4:len(ave_mag)-4]
#magwewant = ave_mag[16:len(ave_mag)-16]
#kbtwewant = kbt[16:len(ave_mag)-16]
initial_guess = np.array([0.001,0.0001,-0.0001,0.005,0.01,-0.005,-0.001,0.02,0.03,0.04,0.05,-0.02,-0.03,-0.04,-0.05,0.1,-0.1,0.075,-0.075])
bestguess = 0

bestcount = 0
for i in range (len(initial_guess)):
    count = 0
    p0 = initial_guess[i]
    print(p0)
# use computational method to find out best initial guess
# try to figure out initial guess
    params = curve_fit(func, kbt, ave_mag,(p0,p0,p0,p0,p0),maxfev = 100000)

    [a,b,c,e,f] = params[0]
#params = curve_fit(func, kbt, ave_mag)


    print a,b,c,e,f
    y = func(kbt,a,b,c,e,f)
    for j in range (len(y)):
        if (abs((y[j]-ave_mag[j])/ave_mag[j])<0.01):
            count+=1

    if count>bestcount:
        bestguess = i

print(i)
now_guess = initial_guess[i]
print(now_guess)
step = now_guess/100
kth_step = now_guess-50*step
bestcount = 0
print(kth_step)
print(now_guess+50*step)
while (abs(kth_step)<=abs(now_guess+50*step)):
    count = 0
    p0 = kth_step
    params = curve_fit(func, kbt, ave_mag, (p0, p0, p0, p0, p0),maxfev = 100000)

    [a, b, c, e, f] = params[0]
    # params = curve_fit(func, kbt, ave_mag)


    print a, b, c, e, f
    y = func(kbt, a, b, c, e, f)
    for j in range(len(y)):
        if (abs((y[j] - ave_mag[j]) / ave_mag[j]) < 0.01):
            count += 1

    if count > bestcount:
        bestguess = kth_step
    print(kth_step)
    kth_step+=step
print(bestguess)
p0 = initial_guess[bestguess]
params = curve_fit(func, kbt, ave_mag,(p0,p0,p0,p0,p0),maxfev = 100000)

[a,b,c,e,f] = params[0]


print a,b,c,e,f
y = func(kbt,a,b,c,e,f)
print(zip(kbt,y))
plt.plot(kbt,y,'ro-')
plt.plot(kbt,ave_mag,'o')
plt.show()


