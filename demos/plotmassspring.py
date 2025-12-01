import numpy as np
data = np.loadtxt('output_test_ode.txt', usecols=(0, 1, 2))
# print (data)

import matplotlib.pyplot as plt


import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--tend_relative",type=str, required=True)
parser.add_argument("--steps",type=str, required=True)
parser.add_argument("--algorithm",type=str, required=True)
args = parser.parse_args()


plt.plot(data[:,0], data[:,1], label='position')
plt.plot(data[:,0], data[:,2], label='velocity')
plt.xlabel('time')
plt.ylabel('value')
plt.title('Mass-Spring System Time Evolution, tend = ' + args.tend_relative + "pi, steps = " +args.steps + " , " + args.algorithm)
plt.legend()
plt.grid()
plt.savefig("Time_evolution.png",dpi=300)
plt.show()
plt.close()

plt.plot(data[:,1], data[:,2], label='phase plot')
plt.xlabel('position')
plt.ylabel('velocity')
plt.title('Mass-Spring System Phase Plot, tend = ' + args.tend_relative + "pi, steps = " +args.steps + " , " + args.algorithm)
plt.legend()
plt.grid()
plt.savefig("Phase_plot.png",dpi=300)
plt.show()
plt.close()