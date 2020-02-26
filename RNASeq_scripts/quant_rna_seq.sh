#!/bin/bash

#############################
# les directives Slurm vont ici:

# Your job name (displayed by the queue)
#SBATCH -J RNA-Seq analysis

# walltime (hh:mm::ss)
#SBATCH -t 02:00:00

# Specify the number of nodes(nodes=) and the number of cores per nodes(tasks-pernode=) to be used
#SBATCH -N 4
#SBATCH --tasks-per-node=32

# change working directory
# SBATCH --chdir=.

# fin des directives PBS
#############################



STAR --runThreadN 10 --genomeLoad NoSharedMemory --genomeDir tempstargenomedir --readFilesIn /data/dnb02/galaxy_db/files/013/734/dataset_13734635.dat /data/dnb02/galaxy_db/files/013/734/dataset_13734637.dat --outSAMtype BAM SortedByCoordinate --twopassMode None --quantMode - --outSAMstrandField intronMotif --outSAMattrIHstart 1 --outSAMattributes NH HI AS nM NM MD MC jM jI ch --outSAMprimaryFlag OneBestScore --outSAMmapqUnique 60 --outSAMunmapped Within --outFilterType Normal --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 10 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.3 --outFilterMismatchNoverReadLmax 1.0 --outFilterScoreMin 0 --outFilterScoreMinOverLread 0.66 --outFilterMatchNmin 0 --outFilterMatchNminOverLread 0.66 --outSAMmultNmax -1 --outSAMtlen 1 --outBAMsortingThreadN 10 --outBAMsortingBinsN 50 --limitBAMsortRAM 143360000000