import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

# We will run this script from the build directory:
#   python3 ../demos/plot_rc.py --tend 0.01 --steps 2000

parser = argparse.ArgumentParser()
parser.add_argument("--tend", type=float, required=True)
parser.add_argument("--steps", type=int, required=True)
args = parser.parse_args()

methods = ["explicit", "improved", "implicit", "CN"]
colors = {
    "explicit": "tab:blue",
    "improved": "tab:orange",
    "implicit": "tab:green",
    "CN": "tab:red",
}

data = {}

# load data from the text files in the current (build) directory
for m in methods:
    fname = f"rc_{m}.txt"

    if os.path.exists(fname):
        t, uc, _ = np.loadtxt(fname, unpack=True)
        data[m] = (t, uc)

if not data:
    print("No RC data files found (explicit_rc.txt, ...). Run test_rc first.")
    raise SystemExit

# where to save plots (relative to this script location)
script_dir = os.path.dirname(os.path.abspath(__file__))
plots_dir = os.path.join(script_dir, "plots", "RC")
os.makedirs(plots_dir, exist_ok=True)

# --------- Plot 1: capacitor voltage vs time (all methods) ----------
plt.figure()
for m, (t, uc) in data.items():
    plt.plot(t, uc, label=m, linewidth=1.0)
plt.xlabel("time t [s]")
plt.ylabel("capacitor voltage U_C(t)")
plt.title(f"RC circuit, tend = {args.tend}, steps = {args.steps}")
plt.grid(True)
plt.legend()

out1 = os.path.join(plots_dir, f"RC_Time_evolution_tend_{args.tend}_steps_{args.steps}.png")
plt.savefig(out1, dpi=300)
print("Saved:", out1)

# (optional) zoomed view of the beginning
plt.figure()
for m, (t, uc) in data.items():
    plt.plot(t, uc, label=m, linewidth=1.0)
plt.xlabel("time t [s]")
plt.ylabel("U_C(t)")
plt.title("RC circuit (zoomed, first 2 ms)")
plt.xlim(0.0, min(args.tend, 0.002))
plt.grid(True)
plt.legend()

out2 = os.path.join(plots_dir, f"RC_Time_evolution_zoom_tend_{args.tend}_steps_{args.steps}.png")
plt.savefig(out2, dpi=300)
print("Saved:", out2)

# do not call plt.show() in WSL (no GUI)
