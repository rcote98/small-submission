import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

file = Path("sims/equi/data.xvg")
df = pd.read_csv(file, sep="\s+", skiprows=1, names=["step", "etot", "temp"])

# ~ apply conversion factors ~

# from KJ/mol to eV
df["etot"] = df["etot"] * 0.010364

# generate plot

fig, ax1 = plt.subplots()

color = 'royalblue'
ax1.set_xlabel('Step')
ax1.set_ylabel('Temperature', color=color)
ax1.plot(df["step"], df["temp"], color=color)
ax1.tick_params(axis='y', labelcolor=color)
#ax1.set_ylim([290, 310])

ax2 = ax1.twinx()  
color = 'tab:red'
ax2.set_ylabel('Total Energy (eV)', color=color)
ax2.plot(df["step"], df["etot"], color=color, linewidth=3)
ax2.tick_params(axis='y', labelcolor=color)
#ax2.set_ylim([-2140, -1800])

fig.tight_layout()
plt.savefig("figures/equi.png")
plt.show()