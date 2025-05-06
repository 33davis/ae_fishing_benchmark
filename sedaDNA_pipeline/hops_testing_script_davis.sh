#!/usr/bin/bash -l
#PBS -N hops_fullprocessing_synbaseiterative__euks
#PBS -M emily.davis@utas.edu.au
# Send mail for (a)borted, (b)eginning, (e)nd job
#PBS -m abe
# Select 1 node with 128 cpus per node
#PBS -l select=1:ncpus=128:mem=990G
# Allow the job to run for up to 24 hours
#PBS -l walltime=72:00:00

## BEGIN SCRIPT
cd ~
date 

### load necessary modules and conda environment
#module load rosalind/1.0 
module load Anaconda3/2024.02-1

source activate hops_env


### script to run HOPS - this will run MALTExtract and HOPS together

### to run full pipeline from MALT thru HOPS, need to specify "-m full" and have -input "/path/to/postdedupe/*.fastq"


### to run just MaltExtract and  HOPS postprocessing specify "-m me_po" and have input "/path/to/post_maltrun/*rma6"
 
hops -Xmx990G -input /u/davisee/Scratch/synthetic_files/fq_input/single_fq/20_09_2024-MARES_BAR_PR2v5_20240216/MALT/*.rma6 -output /u/davisee/Scratch/synthetic_files/fq_input/single_fq/20_09_2024-MARES_BAR_PR2v5_20240216/MALT/mepo_results -m me_po -c /u/davisee/Scratch/scripts/configFile.txt


module purge

#END OF SCRIPT
