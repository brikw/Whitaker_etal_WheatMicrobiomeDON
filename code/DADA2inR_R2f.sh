#!/bin/bash
#SBATCH --job-name="DADA2inR_fungi"
#SBATCH --time=16:00:00   # days-HH:MM:SS walltime 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96   # 48 processor core(s) per node X 2 threads per core
#SBATCH --mem=400G
#SBATCH --partition=ceres    # partition type
#SBATCH --account=mpm_fusarium_genomics
#SBATCH --mail-user=briana.whitaker@usda.gov
#SBATCH --mail-type=BEGIN,END,FAIL

cd /project/mpm_fusarium_genomics/whitaker/niko19/fungi
date
pwd
module load r/4.4.1
module load gdal/3.6.4 #one version older than avail on scinet

Rscript DADA2_hpcc_R2fungi.R


date


#legacy partition? partition=short
