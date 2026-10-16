#!/bin/bash
#SBATCH --job-name=mc_stage2_fsed
#SBATCH --output=logs/mc_stage2_%A_%a.out
#SBATCH --error=logs/mc_stage2_%A_%a.err
#SBATCH --time=06:40:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from 2_HPRC_Monte_Carlo}"
# Activate your Python environment before submitting; adapt resources for your cluster.
python run_second_stage_one_job.py
