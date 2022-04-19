#!/bin/bash
#
# Tutorial 1: Basic Slurm example batch script
#
# This script simply runs a program (xthi) to
# do a pointless calculation for 60 seconds.
# We provide xthi on our Ubuntu Linux platform.
#
#SBATCH --job-name=tut1_basic
#SBATCH --partition=short
#SBATCH --mem=100G

srun lmp -in run2.ligation.lam 
srun lmp -in run3.ligation.lam 
srun lmp -in run5.ligation.lam 
srun lmp -in run6.ligation.lam 
srun lmp -in run11.ligation.lam 
srun lmp -in run21.ligation.lam 
srun lmp -in run26.ligation.lam 

