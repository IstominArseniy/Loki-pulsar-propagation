import numpy as np
from matplotlib import pyplot as plt

xs = []
th1 = []
th2 = []

file_name = "./bin/my_output/my_run_theta.dat"
with open(file_name) as f:
    lines = f.readlines()
    for line in lines:
        tmp = line.split(', ')
        xs.append(float(tmp[0]))
        th1.append(float(tmp[1]))
        th2.append(float(tmp[2]))
f1 = plt.figure()
f2 = plt.figure()
ax1 = f1.add_subplot(111)
ax2 = f2.add_subplot(111)

ax1.plot(xs, th1)
f1.savefig('x_th1.png')
#ax2.scatter(xs[::10], th2[::10], marker='.')
ax2.plot(xs, th2)
f2.savefig('x_th2.png')



