import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

def shift_angle(angle):
    angle = np.asarray(angle)
    scalar_in = (angle.ndim==0)
    angle_vals = np.array(angle, copy=False, ndmin=1, dtype=float)
    for i in range(angle_vals.shape[0]):
        while angle_vals[i] > np.pi / 2:
            angle_vals[i] -= np.pi
        while angle_vals[i] < -np.pi / 2:
            angle_vals[i] += np.pi
    if scalar_in:
        return angle_vals[0]
    else:
        return angle_vals

print("print 0 to see initial values and 1 to see calculated values")
request = int(input())
if request == 0:
    file_name = "./bin/my_output/my_run_0.dat"
else:
    file_name = "./bin/my_output/my_run_1.dat"
phases = []
Is = []
Vs = []
PAs = []
Fs = []
with open(file_name) as f:
    lines = f.readlines()
    for line in lines:
        tmp = line.split(' ')
        if not np.isnan(float(tmp[1])) and not np.isnan(float(tmp[2])) and not np.isnan(float(tmp[3])):
            Is.append(float(tmp[1]))
            Vs.append(float(tmp[2]))
            PAs.append(shift_angle(float(tmp[3]) * np.pi / 180) * 180 / np.pi)
            phases.append(float(tmp[0]))

fig, axs = plt.subplots(2, height_ratios=[1, 4])
Is = np.array(Is)
Vs = np.array(Vs)
Ls = np.sqrt(Is**2 - Vs**2)
data = pd.DataFrame({'phase':phases, 'I':Is.tolist(), 'V':Vs.tolist(), 'L':Ls.tolist(), 'PA':PAs})
data = data.sort_values('phase')
axs[0].scatter(data['phase'], data['PA'], c='black', s=3)
axs[1].plot(data['phase'], data['I'], c='black', label='I')
axs[1].plot(data['phase'], data['V'], c='blue', label='V')
axs[1].plot(data['phase'], data['L'], c='red', label='L')
fig.legend()
fig.show()
fig.savefig('IVLPA_output.png', dpi=400, bbox_inches='tight')
# f1 = plt.figure()
# f2 = plt.figure()
# ax1 = f1.add_subplot(111)
# ax2 = f2.add_subplot(111)
# ax1.plot(phases, Is, label='I')
# ax1.plot(phases, Vs, label='V')
# ax1.legend()
# f1.savefig('I_V_output.png')
# ax2.scatter(phases, PAs, label='PA')
# ax2.legend()
# f2.savefig('PA_output.png')



