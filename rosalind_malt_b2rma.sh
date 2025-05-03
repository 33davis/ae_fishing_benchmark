#!/usr/bin/bash -l
#PBS -N fullrun_definedpath
#PBS -M emily.davis@utas.edu.au
# Send mail for (a)borted, (b)eginning and (e)nd job
#PBS -m abe
# Select 2 nodes with 28 cpus per node
#PBS -l select=2:ncpus=28
# Allow the job to run for up to 24 hours
#PBS -l walltime=24:00:00
# force the job to a specific queue
####PBS -q [destination queue]

## BEGIN SCRIPT
cd ~
date # Print the current time and date to standard output (so you can see when the job started)

module load rosalind/1.0
module load GCCcore/6.4.0 
module load megan/6.24.12


# MALT
## have to use a local copy of the database (CANNOT use path to shared folder databases, as I do not have writing perissions
cd /u/davisee/miniconda2/pkgs/malt-0.62-hdfd78af_0/bin
export PATH="/u/davisee/miniconda2/pkgs/malt-0.62-hdfd78af_0/bin:$PATH"


./malt-run -i /u/davisee/Scratch/synthetic_files/diff_deams/*-dd.gz -a /u/davisee/Scratch/synthetic_files/diff_deams/ -ou /u/davisee/Scratch/synthetic_files/diff_deams/ -at SemiGlobal -f Text -d /u/davisee/Scratch/MARES_BAR_PR2v5_20240216/ -mem load --mode BlastN -t 8 -v


## blast2rma

module load parallel/20210622-GCCcore-10.3.0

input=/u/davisee/Scratch/synthetic_files/diff_deams

ms=50
me=0.01
supp=0
mpi=90

for i in $input/*.kx-dd.gz; do b=${i/.kx-dd.gz/.blastn.gz}; echo blast2rma -r $i -i $b -o $input -v -ms $ms -me $me -f BlastText -supp $supp -mpi $mpi -alg weighted -lcp 80; done | parallel -j 2 --delay 6


module purge

#END OF SCRIPT

