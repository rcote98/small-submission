import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

# read simulation data
file = Path("sims/prod/rdf.xvg")
df = pd.read_csv(file, sep="\s+", names=["r", "g(r)"])
r_sim = df["r"].values*10 # convert to angstroms
g_sim = df["g(r)"].values

# read experimental data
file = Path("figures/rdf_exp.csv")
df = pd.read_csv(file).sort_values(by="x")
r_exp = df["x"].values
g_exp = df["y"].values

# generate plot
fig, ax = plt.subplots()
ax.plot(r_sim, g_sim, label="Simulation")
ax.plot(r_exp, g_exp, label="Experiment")
ax.set_ylim([0, 3])
ax.set_xlim([0, 10])
ax.set_xlabel("r (Å)")
ax.set_ylabel("g(r)")
ax.grid(True)
ax.legend()


fig.tight_layout()
plt.savefig("figures/rdf.png")
plt.show()