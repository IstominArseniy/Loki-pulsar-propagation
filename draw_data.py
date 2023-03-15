import numpy as np
from matplotlib import pyplot as plt
print("print 0 to see initial values and 1 to see calculated values")
request = input()
if request == 0:
    file_name = "./bin/my_output/my_run_0.dat"
else:
    file_name = "./bin/my_output/my_run_1.dat"
phases = []
Is = []
Vs = []
PAs = []
with open(file_name) as f:
    lines = f.readlines()
    for line in lines:
        tmp = line.split(' ')
        if not np.isnan(float(tmp[1])) and not np.isnan(float(tmp[2])) and not np.isnan(float(tmp[3])):
            Is.append(float(tmp[1]))
            Vs.append(float(tmp[2]))
            PAs.append(float(tmp[3]))
            phases.append(float(tmp[0]))
f1 = plt.figure()
f2 = plt.figure()
ax1 = f1.add_subplot(111)
ax2 = f2.add_subplot(111)
ax1.plot(phases, Is, label='I')
ax1.plot(phases, Vs, label='V')
ax1.legend()
f1.savefig('I_V_output.png')
ax2.plot(phases, PAs, label='PA')
ax2.legend()
f2.savefig('PA_output.png')



