import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

file = Path("benchmark/output/bench.csv")
df = pd.read_csv(file)

# calculate speedup and efficiency
dfs = df.groupby("num_threads")[["cpu_time", "wall_time"]].std().reset_index()
df = df.groupby("num_threads")[["cpu_time", "wall_time"]].mean().reset_index()
df["speedup"] = df["wall_time"].max() / df["wall_time"]
df["efficiency"] = df["speedup"] / df["num_threads"]

# generate plot
fig, ax1 = plt.subplots()

color = 'royalblue'
ax1.set_xlabel('Number of Threads')
ax1.set_ylabel('Speedup', color=color)
ax1.plot(df["num_threads"], df["speedup"], color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  
color = 'tab:red'
ax2.set_ylabel('Wall Time', color=color)
ax2.errorbar(df["num_threads"], df["wall_time"], yerr=dfs["wall_time"], color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()
plt.savefig("figures/bench.png")
plt.show()