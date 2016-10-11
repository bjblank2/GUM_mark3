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


