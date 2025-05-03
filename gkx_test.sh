#!/usr/bin/bash -l
#PBS -N rosalind_Kdmb_test
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
module load rosalind komplexity bbtools megan

# Komplexity


# Navigate to directory containing collapsed sample files
cd /u/davisee/Scratch/Exp382-Deep_collapsed/U1538test

# Create new directory to input all results
## d == date testing script
d=$(date +"%d_%m_%Y")
mkdir "/u/davisee/Scratch/Exp382-Deep_collapsed/U1538test/$d"

##### Komplexity ######
# Run Komplexity
for sample in *.collapsed; do kz --filter --threshold 0.55 < $sample > $sample.kx; done

# Move Komplexity results to new directory
mv *.kx "/u/davisee/Scratch/Exp382-Deep_collapsed/U1538test/$d/"

####### Dedupe #########
cd "/u/davisee/Scratch/Exp382-Deep_collapsed/U1538test/$d"
for sample in *.kx; do dedupe.sh in=$sample out=$sample-dd.gz outd=$sample-duplicates.fa ac=f; done

####### MALT ###########
## Use local copy of the database
cd /u/davisee/miniconda2/pkgs/malt-0.62-hdfd78af_0/bin
export PATH="/u/davisee/miniconda2/pkgs/malt-0.62-hdfd78af_0/bin:$PATH"
./malt-run -i /u/davisee/Scratch/Exp382-Deep_collapsed/U1538test/$d/*-dd.gz -a /u/davisee/Scratch/Exp382-Deep_collapsed/U1538test/$d/ -ou /u/davisee/Scratch/Exp382-Deep_collapsed/U1538test/$d/ -at SemiGlobal -f Text -d /u/davisee/Scratch/MARES_BAR_PR2v5_20240216/ -mem load --mode BlastN -t 8 -v

####### blast2rma ###################	
module load parallel/20210622-GCCcore-10.3.0

input=/u/davisee/Scratch/Exp382-Deep_collapsed/U1538test/$d/
ms=50
me=0.01
supp=0
mpi=95

for i in $input/*.kx-dd.gz; do b=${i/.kx-dd.gz/.blastn.gz}; echo blast2rma -r $i -i $b -o $input -v -ms $ms -me $me -f BlastText -supp $supp -mpi $mpi -alg weighted -lcp 80; done | parallel -j 2 --delay 6

# Organize results into subdirectory
mkdir "/u/davisee/Scratch/Exp382-Deep_collapsed/U1538test/$d/all_MARES2021_PR2v5"
mv /u/davisee/Scratch/Exp382-Deep_collapsed/U1538test/$d/*.rma6 /u/davisee/Scratch/Exp382-Deep_collapsed/U1538test/$d/all_MARES2021_PR2v5/
mv /u/davisee/Scratch/Exp382-Deep_collapsed/U1538test/$d/*.blastn.gz /u/davisee/Scratch/Exp382-Deep_collapsed/U1538test/$d/all_MARES2021_PR2v5/

module purge 

#END OF SCRIPT 
