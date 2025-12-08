# include/plot_autodiff_legendre.py
import numpy as np
import matplotlib.pyplot as plt
data = np.genfromtxt("legendre.csv", delimiter=',', names=True)
x = data['x']
plt.figure(figsize=(9,6))
for k in range(6):
    plt.plot(x, data[f'P{k}'], label=f'P{k}')
plt.legend(); plt.xlabel('x'); plt.title('P0..P5'); plt.grid(True); plt.savefig('legendre_Ps.png')
plt.figure(figsize=(9,6))
for k in range(6):
    plt.plot(x, data[f'dP{k}'], label=f'dP{k}')
plt.legend(); plt.xlabel('x'); plt.title('dP0..dP5'); plt.grid(True); plt.savefig('legendre_dPs.png')
print("Wrote legendre_Ps.png and legendre_dPs.png")
