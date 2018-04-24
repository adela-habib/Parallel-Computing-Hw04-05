#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt


IOtime = np.loadtxt("IOtime.txt")
time = np.loadtxt("Time.txt", delimiter =', ')
lab = [0, 4, 8, 32, 64]

ranks = IOtime[:, 0]
threads = IOtime[:, 1]
exct=IOtime[:, 2]

plt.figure(3)
plt.plot(np.log2(ranks), exct, 'o-')
plt.xlabel("Number of ranks")
plt.ylabel("Execution time [s]")
plt.xticks([6,7,8,9,10,11,12,13,14])

plt.figure(1)
c = ["red", "orange", "green", "blue"]
for i in range(4):
	x = []
	y = []
	k = i
	for j in range(4):
		x.append(time[k, 1]*time[k, 2])
		y.append(time[k, 3])
        	k += 5
	plt.plot(np.log2(x),y, 'o-', color = c[i], label = str(lab[i]) +"-pthreads")
plt.legend(frameon=False)
plt.xlabel(r"log$_{2}$(ranks*pthreads)")
plt.ylabel("execution time [s]")
plt.xticks([7,8,9,10,11,12,13,14])

plt.figure(2)
for i in range(4):
        x = []
        y = []
        k = i+5  
        for j in range(3):
                print k, time[k, 1],time[i, 3]
                x.append(time[k, 1]*time[k, 2])
                y.append(time[i, 3]/time[k, 3])
                k += 5
        plt.plot(np.log2(x),y, 'o-', color = c[i], label = str(lab[i]) +"-pthreads")
plt.legend(loc = "upper left", frameon=False)
plt.xlabel(r"log$_{2}$(ranks*pthreads)")
plt.ylabel("relative speedUp")
plt.xticks([9,10,11,12,13,14])

plt.show()


