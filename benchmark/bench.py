from pathlib import Path
from shutil import copy
import pandas as pd
import subprocess
import os

PRESET_FOLD = Path("benchmark/preset")
NUM_THREADS = [1,2,4,8,16]
CV_REPEATS  = 5

RUN_SIMS = 0
GAT_DATA = 1

# save current directory
cwd = Path.cwd()
df = pd.DataFrame(columns=["num_threads", "cv", "cpu_time", "wall_time"])

for nt in NUM_THREADS:
    for cv in range(CV_REPEATS):

        # create output folder
        out = Path(f"benchmark/output/runs/{nt}c_{cv}cv")
        out.mkdir(parents=True, exist_ok=True)


        if RUN_SIMS:
            # copy over all the files from preset
            for file in PRESET_FOLD.iterdir():
                copy(file, out / file.name)

            # change directory to output folder
            os.chdir(out)

            # run grompp using subprocess
            subprocess.run(["gmx", "grompp", "-f", "bench.mdp", "-c", "waterbox-equi.gro", "-p", "topol.top", "-t", "state.cpt"])

            # run mdrun using subprocess
            subprocess.run(["gmx", "mdrun", "-nt", str(nt), "-pin", "on"])
            
            # change directory back to original
            os.chdir(cwd)

        if GAT_DATA:

            with open(out / "md.log", "r") as f:
                lines = f.readlines()
                df.loc[len(df)] = [nt, cv, float(lines[-5].split()[1]), float(lines[-5].split()[2])]

if GAT_DATA:
    df.to_csv("benchmark/output/bench.csv", index=False)