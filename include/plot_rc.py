# include/plot_rc.py
import csv
import numpy as np
import matplotlib.pyplot as plt

def read_csv(fname):
    t, v, th = [], [], []
    with open(fname) as f:
        r = csv.DictReader(f)
        for row in r:
            t.append(float(row['t']))
            v.append(float(row['v']))
            th.append(float(row['theta']))
    return np.array(t), np.array(v), np.array(th)

files = [("rc_explicit.csv","Explicit Euler"),
         ("rc_improved.csv","Improved Euler"),
         ("rc_crank.csv","Crank-Nicolson")]

plt.figure(figsize=(9,5))
for fname,label in files:
    t,v,th = read_csv(fname)
    plt.plot(t, v, label=label)
plt.xlabel("t")
plt.ylabel("capacitor voltage v(t)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("rc_v_compare.png")
print("Wrote rc_v_compare.png")
