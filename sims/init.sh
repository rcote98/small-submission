# create a base directory
mkdir init

# switch to the base directory
cd init

# generate a water box of 5nm
gmx solvate -box 5 5 5 -scale 0.52 -o waterbox.gro -cs spc216.gro 

# generate a topology for the water box
gmx pdb2gmx -f waterbox.gro << EOF
1 1
EOF

# Go back to the main directory
cd ..
