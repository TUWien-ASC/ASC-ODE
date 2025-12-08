# python/plot_mass_spring.py
import csv
import matplotlib.pyplot as plt
import numpy as np

def read_csv(filename):
    t, x, v = [], [], []
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            t.append(float(row['t']))
            x.append(float(row['x']))
            v.append(float(row['v']))
    return np.array(t), np.array(x), np.array(v)

files = [
    ("mass_explicit.csv", "Explicit Euler"),
    ("mass_improved.csv", "Improved Euler"),
    ("mass_crank.csv", "Crank-Nicolson (Picard)"),
]

plt.figure(figsize=(9,5))
for fname, label in files:
    t, x, v = read_csv(fname)
    plt.plot(t, x, label=label)
plt.xlabel("t")
plt.ylabel("position x(t)")
plt.legend()
plt.title("Mass-Spring: position vs time")
plt.grid(True)
plt.tight_layout()
plt.savefig("mass_position_compare.png")
print("Wrote mass_position_compare.png")

# phase plot
plt.figure(figsize=(6,6))
for fname, label in files:
    t, x, v = read_csv(fname)
    plt.plot(x, v, label=label)
plt.xlabel("position x")
plt.ylabel("velocity v")
plt.legend()
plt.title("Mass-Spring: phase plot")
plt.grid(True)
plt.tight_layout()
plt.savefig("mass_phase_compare.png")
print("Wrote mass_phase_compare.png")
