#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 
#SBATCH --cpus-per-task=1
#SBATCH --time=65:00:00
#SBATCH --mem=2GB
#SBATCH --export=ALL
#SBATCH --job-name=M2x4_0.1_v3

export JOB_WORK_DIR="${SLURM_JOB_NAME}.${SLURM_JOB_ID%%.*}.$(date +%y%m%d_%H%M%S)"
export TMP_WORK_DIR=${SCRATCH}/${JOB_WORK_DIR}
mkdir -pv ${TMP_WORK_DIR}
cd ${TMP_WORK_DIR}

cp ${SLURM_SUBMIT_DIR}/input0M.jl ${TMP_WORK_DIR}
cp ${SLURM_SUBMIT_DIR}/init_X_traj0.dat ${TMP_WORK_DIR}
cp ${SLURM_SUBMIT_DIR}/init_P_traj0.dat ${TMP_WORK_DIR}
$HOME/PhD/code/julia-1.6.1/bin/julia input0M.jl 

mkdir -vp "$SLURM_SUBMIT_DIR"/
tar czvf "$SLURM_SUBMIT_DIR"/data0M.tgz data0M
cd ${SCRATCH}
rm -rvf "${TMP_WORK_DIR}"
exit
