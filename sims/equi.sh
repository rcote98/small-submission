# create a base directory
mkdir equi

# switch to the base directory
cd equi

# copy the geometry and topology from the init directory
cp ../emin/confout.gro waterbox-emin.gro
cp ../emin/topol.top . 

# copy the mdp corresponding mdp file
cp ../mdps/equi.mdp .

# run the energy minimization
gmx grompp -f equi.mdp -c waterbox-emin.gro -p topol.top
gmx mdrun 

# get the energy and temperature for plotting
gmx energy -f ener.edr -o data -xvg none << EOF
7 9 0
EOF

# Go back to the main directory
cd ..

