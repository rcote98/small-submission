# create a base directory
mkdir prod

# switch to the base directory
cd prod

# copy the geometry and topology from the init directory
cp ../equi/confout.gro waterbox-equi.gro
cp ../equi/state.cpt .
cp ../emin/topol.top . 

# copy the mdp corresponding mdp file
cp ../mdps/prod.mdp .

# run the energy minimization
gmx grompp -f prod.mdp -c waterbox-equi.gro -p topol.top -t state.cpt
gmx mdrun 

# calculate the dielectric constant
gmx dipoles -f traj.trr -xvg none << EOF
0
EOF &> dipoles.log

# generate an index file for the oxygen atoms
gmx make_ndx -f topol.tpr << EOF
a OW
q
EOF

# calculate the radial distribution function
gmx rdf -f traj.trr -s topol.tpr -n index.ndx -xvg none -rmax 1 << EOF
3 3
EOF &> rdf.log

# Go back to the main directory
cd ..

