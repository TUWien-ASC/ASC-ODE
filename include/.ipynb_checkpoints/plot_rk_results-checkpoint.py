# include/plot_rk_results.py
import numpy as np
import matplotlib.pyplot as plt
def read(fname): return np.genfromtxt(fname, delimiter=',', names=True)
fe = read("pendulum_euler.csv")
f2 = read("pendulum_rk2.csv")
f4 = read("pendulum_rk4.csv")
plt.figure(figsize=(10,5))
plt.plot(fe['t'], fe['alpha'], label='Euler')
plt.plot(f2['t'], f2['alpha'], label='RK2')
plt.plot(f4['t'], f4['alpha'], label='RK4')
plt.legend(); plt.xlabel('t'); plt.title('Pendulum alpha(t)'); plt.grid(True); plt.savefig('pendulum_alpha_compare.png')
plt.figure(figsize=(6,6))
plt.plot(fe['alpha'], fe['alpha_dot'], label='Euler')
plt.plot(f2['alpha'], f2['alpha_dot'], label='RK2')
plt.plot(f4['alpha'], f4['alpha_dot'], label='RK4')
plt.legend(); plt.xlabel('alpha'); plt.ylabel('alpha_dot'); plt.title('Pendulum phase'); plt.grid(True); plt.savefig('pendulum_phase_compare.png')
print("Wrote pendulum plots")
# error vs tau
le = read("linear_error_vs_tau.csv")
taus = le['tau']; plt.figure(); plt.loglog(taus, le['err_euler'], 'o-', label='Euler'); plt.loglog(taus, le['err_rk2'], 's-', label='RK2'); plt.loglog(taus, le['err_rk4'], '^-', label='RK4'); plt.legend(); plt.xlabel('tau'); plt.ylabel('error'); plt.grid(True); plt.savefig('linear_error_vs_tau.png'); print("Wrote linear_error_vs_tau.png")
