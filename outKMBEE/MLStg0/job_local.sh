#!/bin/bash
#SBATCH  -J SIM
#SBATCH -p fe13
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --time=168:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

echo "starting $SLURM_JOBID at `date` on `hostname`"
echo "my home dir is $HOME"

python /home/icsefs01/spec1084/software/MPTS/modsKM/simbruteKM.py --inp_file /home/icsefs01/spec1084/software/MPTS/outKMBEE/MLStg0/siminp.yaml --batch_file /home/icsefs01/spec1084/software/MPTS/outKMBEE/prpsl_0.txt --run_style recursive

echo "ended at `date` on `hostname`"
exit 0
