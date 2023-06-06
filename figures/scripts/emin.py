import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

file = Path("sims/emin/emin.log")

# read the output file from the emin run
with open(file) as f:
    df = pd.DataFrame(columns=["step", "dmax", "epot", "fmax"])
    line = f.readline()
    while line:
        line = f.readline()
        if line.strip().startswith("Step="):
            step = int(line.split()[1][:-1])
            dmax = float(line.split()[3])
            epot = float(line.split()[6])
            fmax = float(line.split()[8][:-1])
            df.loc[len(df)] = [step, dmax, epot, fmax]

# ~ apply conversion factors ~

# from KJ/mol to eV
df["epot"] = df["epot"] * 0.010364

# force, from KJ/mol/nm to eV/A
df["fmax"] = df["fmax"] * 0.010364 / 10

# print minimum energy and force
print("Minimum pot_energy: ", df["epot"].min())
print("Minimum max_force: ", df["fmax"].min())

# generate plot

fig, ax1 = plt.subplots()

color = 'royalblue'
ax1.set_xlabel('Step')
ax1.set_ylabel('Maximum Force (eV/A)', color=color)
ax1.plot(df["step"], df["fmax"], color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim([0, 10])

ax2 = ax1.twinx()  
color = 'tab:red'
ax2.set_ylabel('Potential Energy (eV)', color=color)
ax2.plot(df["step"], df["epot"], color=color, linewidth=3)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim([-2140, -1800])

fig.tight_layout()
plt.savefig("figures/emin.png")
plt.show()