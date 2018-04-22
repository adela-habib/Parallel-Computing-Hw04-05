#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

time = np.loadtxt("Time.txt", usecols = (0, 1, 2, 3), delimiter=', ')
print time.shape

plt.figure()

k = 0
for i in range(4):
	x = []
	y = []
	for i in range(4):
		x.append(time[i+k,1]*time[i+k,2])
		y.append(time[i+k, 3])
	
	k += 4
	plt.plot(x, y, 'o')

plt.show()
exit(0)



time = np.reshape(time, (4, 20))
print time
c = ["red", "orange", "green", "blue"]
k = 3
for i in range(4):
	x = []
	y = []
	for j in range(4):
		print  k, time[j, k-2],time[j, k-1], time[j,k]   
		x.append(time[j, k-1]*time[j, k-2])
		y.append(time[j,k])
        plt.plot(x, y, 'o', color = c[i], label = str(int(time[j, k-1])) +"-pthreads")
	k = k+3+1
plt.legend()
plt.show()


