#!/bin/bash
#SBATCH -N 1
#SBATCH -p cascade.p
#SBATCH --gres=gpu:1
#SBATCH --mincpus=20
#SBATCH --time=24:00:00
#SBATCH --job-name=EnterJobName

module use /hits/sw/its/doserbd/cascade/modules/all
module load GROMACS/2021-fosscuda-2019b
echo "1. input_name, 2. nvt_name, 3. npt_name, 4. pull_name"
input_name=$1
nvt_name=$2
npt_name=$3
pull_name=$4

checkpoint=./$4.cpt

export CYCLE=24
START=$(date +"%s")

if ! test -f "$checkpoint";then
	gmx grompp -f pull.mdp -c $npt_name.gro -r $npt_name.gro -o $pull_name.tpr -p $input_name.top -maxwarn 2 -n $input_name.ndx
fi 
ext_x=_pullx.xvg
ext_f=_pullf.xvg
gmx mdrun -v -deffnm $pull_name -px $pull_name$ext_x -pf $pull_name$ext_f  -maxh 24 -cpi $pull_name.cpt -cpo $pull_name.cpt -dlb yes -ntomp 2 -ntmpi 10 -pme gpu -bonded gpu -npme 1

END=$(date +"%s")
echo "$(((END-START))) seconds ran"
echo "$(((END-START)/3600)) full hours ran"
let "CYCLE--"
if [ $(((END-START)/3600)) -lt $CYCLE ]
then
        echo "last cycle was just $(((END-START)/3600))h long and finished"
else
        sbatch resubmit_pull1.sh $1 $2 $3 $4
fi

