import numpy as np
from matplotlib import pyplot as plt


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

xs = []
th1 = []
th2 = []
freq = []
Fs = []
print("Print phase")
phase = float(input())
phase_str = f"{phase:.6f}"
file_name = "./bin/my_output/my_run_theta_data/my_run_theta_" + phase_str + ".dat"
with open(file_name) as f:
    lines = f.readlines()
    for line in lines:
        tmp = line.split(', ')
        xs.append(float(tmp[0]))
        th1.append(float(tmp[1]))
        th2.append(float(tmp[2]))
        freq.append(float(tmp[3]))
        Fs.append(float(tmp[4]))
xs = np.array(xs)
th1 = np.array(th1)
th2 = np.array(th2)
freq = np.array(freq)
Fs = np.array(Fs)
f1 = plt.figure()
f2 = plt.figure()
f3 = plt.figure()
ax1 = f1.add_subplot(111)
ax2 = f2.add_subplot(111)
ax3 = f3.add_subplot(111)
ax1.plot(xs, shift_angle(th1), label=r'$\theta$')
ax1.plot(xs, shift_angle(Fs), label=r'$\beta + \delta$')
ax1.plot(xs, shift_angle(np.array(Fs) - np.array(th1)), label=r'$\beta + \delta - \theta$')
ax1.plot()
ax1.legend()
f1.savefig('x_th1.png')
#ax2.scatter(xs[::10], th2[::10], marker='.')
ax2.plot(xs, th2)
#ax2.plot(xs, 0.05 + 0.01 * np.sin(freq * xs))
f2.savefig('x_th2.png')
ax3.plot(xs, freq)
f3.savefig('tests.png')



