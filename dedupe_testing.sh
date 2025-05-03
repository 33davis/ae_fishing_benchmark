#!/usr/bin/bash -l
#PBS -N komp_dedudpe_only_testing
#PBS -M emily.davis@utas.edu.au
# Send mail for (a)borted, (b)eginning and (e)nd job
#PBS -m abe
# Select 1 nodes with 28 cpus per node
#PBS -l select=1:ncpus=28
# Allow the job to run for up to 24 hours
#PBS -l walltime=24:00:00
# force the job to a specific queue
####PBS -q [destination queue]

## BEGIN SCRIPT

date # Print the current time and date to standard output (so you can see when the job started)

module load rosalind/1.0
module load gcc-env/6.4.0
module load rosalind adapterremoval komplexity bbtools megan

#add in memory stipulations
export JAVA_TOOL_OPTIONS="-Xmx600g -Xms5000m"

# Komplexity
cd /u/davisee/Scratch/Exp382-Deep_collapsed/U1538	

for sample in *.collapsed; do kz --filter --threshold 0.55 < $sample > $sample.kx; done

cd /u/username/Scratch/
mkdir dedupe_testing
mv /u/davisee/Scratch/Exp382-Deep_collapsed/U1538/*.kx /u/davisee/Scrath/dedupe_testing

# dedupe

cd /u/davisee/Scrath/dedupe_testing
for sample in *.kx; do dedupe.sh in=$sample out=$sample-dd.gz outd=$sample-duplicates.fa ac=f; done


module purge

#END OF SCRIPT


