#!/bin/bash

#SBATCH --job-name=trendProfile
#SBATCH --mail-user=kevtan@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
 
#SBATCH --account=stats_dept1
#SBATCH --partition=standard

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1

## 5GB/cpu is the basic share
#SBATCH --mem-per-cpu=1GB

## wall time hours:minutes:seconds
#SBATCH --time=3-00:00:00

###   Load software modules

module load R/4.1.0
module list

####  Commands your job should run follow this line

echo "Running on $SLURM_JOB_NODELIST"
echo "Running in $(pwd)"

Rscript trendProfileNoRecal.R
