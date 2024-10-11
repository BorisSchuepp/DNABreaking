# Workflow for Simulations

Start with the nicked or non-nicked structures generated in StructureGeneration.
1. Adjust the box size to account for stretching of DNA \
gmx editconf -f {Sequence}\_{Nick-type}.gro -o {Sequence}{Nick-type}_newbox.gro -box 7 7 65
2. Solvate the system \
solvate: gmx solvate -cp {Sequence}\_{Nick-type}\_newbox.gro -cs spc216.gro -o {Sequence}_{Nick-type}\_solvated.gro -p {Sequence}\_{Nick-type}.top
3. Add ions \
gmx grompp -f ions.mdp -c {Sequence}\_{Nick-type}\_solvated.gro -p {Sequence}_{Nick-type}.top -o ions.tpr \
gmx genion -s ions.tpr -o {Sequence}\_{Nick-type}\_ions.gro -p {Sequence}\_{Nick-type}.top -pname NA -nname CL -neutral
4. Perform energy minimization
gmx grompp -f minim.mdp -c {Sequence}\_{Nick-type}\_ions.gro -p {Sequence}\_{Nick-type}.top -o {Sequence}\_{Nick-type}\_EM.tpr
gmx mdrun -v -deffnm {Sequence}\_{Nick-type}\_EM
5. Select pulling group (here: terminal residues from both strands on both sides)\
gmx make_ndx -f {Sequence}\_{Nick-type}\_EM.gro -o {Sequence}\_{Nick-type}\_EM.ndx \
use ri 1 | ri 200 for the group PullLeft and ri 100 | ri 101 for the group PullRight
6. Rename the topology file to {Sequence}\_{Nick-Type}\_EM.top
7. Submit simulation 
resubmit_nvt.sh {Sequence}\_{Nick-Type}\_EM {Sequence}\_{Nick-Type}\_NVT {Sequence}\_{Nick-Type}\_NPT {Sequence}\_{Nick-Type}_Pull
8. After NVT, NPT is finished submit reruns of the pulling simulation
resubmit_pull{i}.sh {Sequence}\_{Nick-Type}\_EM {Sequence}\_{Nick-Type}\_NVT {Sequence}\_{Nick-Type}\_NPT {Sequence}\_{Nick-Type}\_Pull_{i}
where i is a counter indicating the i-th rerun of the simulation

The .mdp files have been adapted from http://www.mdtutorials.com/gmx/lysozyme/index.html and are provided in this directory\
The scripts to run the simulations are also provided here. Run input files for the pulling simulations can be found in StructureGeneration/InputStructrues/ as .tpr files. \
The resulting trajectories can be obtained from Zenodo as .xtc files with the names {Sequence}\_{Nick-Type}_{Rerun counter}.xtc. For the simulations under different forces
the pull.mdp has been changed in the pull-coord1-k and pull-coord2-k field.

