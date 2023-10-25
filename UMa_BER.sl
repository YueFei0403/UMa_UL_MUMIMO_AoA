#!/bin/bash -l
#SBATCH --job-name=UMa_BER
#SBATCH --account=def-rsadve
#SBATCH --time=10:10:30          # adjust this to match the walltime of your job
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --output=UMa_BER-%j.out
#SBATCH --cpus-per-task=1      # adjust this if you are using parallel commands
#SBATCH --mem-per-cpu=2G             # adjust this according to the memory requirement per node you need
#SBATCH --mail-user=y.fei@mail.utoronto.ca # adjust this to match your email address
#SBATCH --mail-type=END


######  Run Single-thread Job
module load matlab
matlab -nodisplay -singleCompThread -r UMa_UL_MUMIMO_AoA.m

######  Run Parallel Job
# MAIN="AoA_MultiP"
# NWORKERS=${SLURM_NTASKS}
# NTHREADS=${SLURM_CPUS_PER_TASK}
# ARGS="($NWORKERS,$NTHREADS)"
# MAIN_WITH_ARGS=${MAIN}${ARGS}

# module load matlab
# matlab -batch "${MAIN_WITH_ARGS}"