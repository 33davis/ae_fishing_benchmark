#!/usr/bin/bash -l
#PBS -N rosalind_malt-index-MARESPR2
#PBS -M emily.davis@utas.edu.au
# Send mail for (a)borted, (b)eginning and (e)nd job
#PBS -m abe
# Select 2 nodes with 28 cpus per node
#PBS -l select=2:ncpus=28
# Allow the job to run for up to 24  hours
#PBS -l walltime=24:00:00
# force the job to a specific queue
####PBS -q [destination queue]

## BEGIN SCRIPT
#need to activate malt first? need to put this into the documentation
#check that malt is in miniconda2/bin?
date # Print the current time and date to standard output (so you can see when the job started)

#load necessary modules
module load rosalind/1.0
module load gcc-env/6.4.0
module load rosalind adapterremoval komplexity bbtools megan


# Malt-index

# MALT
# make indices directory first
cd /u/davisee/Scratch/
mkdir MARES_BAR_20240216

cd /u/davisee/miniconda2/pkgs/malt-0.62-hdfd78af_0/bin
export PATH="/u/davisee/miniconda2/pkgs/malt-0.62-hdfd78af_0/bin:$PATH"
./malt-build -i /data/imas_projects/ancient/share/Databases/MARES_BAR_16012023.fasta.gz --sequenceType DNA --index /u/davisee/Scratch/MARES_BAR_20240216 --threads 16 --verbose

module purge

#END OF SCRIPT
