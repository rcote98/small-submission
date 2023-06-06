# create a base directory
mkdir emin

# switch to the base directory
cd emin

# copy the geometry and topology from the init directory
cp ../init/waterbox.gro .
cp ../init/topol.top .

# copy the mdp corresponding mdp file
cp ../mdps/emin.mdp .

# run the energy minimization
gmx grompp -f emin.mdp -c waterbox.gro -p topol.top
gmx mdrun -v &> emin.log

# Go back to the main directory
cd ..

