#!/bin/bash

#############################
# les directives Slurm vont ici:

# Your job name (displayed by the queue)
#SBATCH -J TPM

# walltime (hh:mm::ss)
#SBATCH -t 00:20:00

# Specify the number of nodes(nodes=) and the number of cores per nodes(tasks-pernode=) to be used
#SBATCH -N 4
#SBATCH --tasks-per-node=32

# change working directory
# SBATCH --chdir=.

# fin des directives PBS
#############################


# useful informations to print
echo "#############################" 
echo "User:" $USER
echo "Date:" `date`
echo "Host:" `hostname`
echo "Directory:" `pwd`
echo "SLURM_JOBID:" $SLURM_JOBID
echo "SLURM_SUBMIT_DIR:" $SLURM_SUBMIT_DIR
echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST
echo "#############################" 

#############################

## Load necessary modules
module load jdk1.8/8u22
module load python/3.7.2 

export LD_LIBRARY_PATH=/gpfs/softs/contrib/apps/gcc/7.3.0/lib64/:/gpfs/softs/contrib/apps/gcc/7.3.0/lib
export LD_LIBRARY_PATH=/gpfs/softs/contrib/apps/python/3.7.2/bin/python3.7/lib/libpython3.7m.so:$LD_LIBRARY_PATH
# list=$(ls subdata/*fq.gz | xargs -n 1 basename | sed 's/\(.*\)_.*/\1/' | sort -u)

# # mkdir -p -m 755 trimm_data 
# # mkdir -p -m 755 trimm_data/quality

# gtf=$(find ./ -not -path '*/\R*' -name "*.gtf")
# gff=$(find ./ -not -path '*/\R*' -name "*gff[3]*")
# fasta=$(find ./ -name "*.fasta")

# if [ ! -d "genome_ind" ] 
#  then  ## Genome indexes are generated if necessary 
#     mkdir -m 755 -p genome_ind 
#     ./STAR \
#     --runThreadN 16 \
#     --runMode genomeGenerate \
#     --genomeDir /gpfs/home/juagarcia/genome_ind \
#     --genomeFastaFiles "$fasta" \
#     --sjdbOverhang 149 \
#     --sjdbGTFfile "$gtf" \
#     --genomeChrBinNbits 12 
# fi  

# mkdir -p -m 755 STAR_Align  
# mkdir -p -m 755 counts
# for I in $list
# do
 
#     # cp subdata/"$I"_1.fq.gz subdata/"$I"_1.fastqsanger 
#     # cp subdata/"$I"_2.fq.gz subdata/"$I"_2.fastqsanger

#     #FastQC/fastqc subdata/"$I"_1.fq subdata/"$I"_2.fq --outdir=trimm_data/quality

#     # java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 128 -phred33 ./subdata/"$I"_1.fq.gz ./subdata/"$I"_2.fq.gz trimm_data/"$I"_1.par.fq.gz "$I"_1.unp.fq.gz trimm_data/"$I"_2.par.fq.gz "$I"_2.unp.fq.gz ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 LEADING:25 TRAILING:25 SLIDINGWINDOW:5:20 MINLEN:50

#     # rm -rf "$I"_1.unp.fq.gz "$I"_2.unp.fq.gz
#./scripts/STAR --genomeLoad LoadAndExit --genomeDir ./bin/genome_ind
## Load genome just once to save RAM memory 
#     ./STAR \ --runThreadN 128 --genomeLoad NoSharedMemory --readFilesIn  --outSAMtype BAM SortedByCoordinate --twopassMode None --quantMode - --outSAMstrandField intronMotif --outSAMattrIHstart 1 --outSAMattributes NH HI AS nM NM MD MC jM jI ch --outSAMprimaryFlag OneBestScore --outSAMmapqUnique 60 --outSAMunmapped Within --outFilterType Normal --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 10 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.3 --outFilterMismatchNoverReadLmax 1.0 --outFilterScoreMin 0 --outFilterScoreMinOverLread 0.66 --outFilterMatchNmin 0 --outFilterMatchNminOverLread 0.66 --outSAMmultNmax -1 --outSAMtlen 1 --outBAMsortingThreadN 10 --outBAMsortingBinsN 50 --limitBAMsortRAM 143360000000
#     --runThreadN 128 \
#     --runMode alignReads \
#     --outFilterMatchNmin 16 \
#     --readFilesCommand zcat \
#     --outSAMtype BAM Unsorted SortedByCoordinate \
#     --genomeDir /gpfs/home/juagarcia/genome_ind \
#     --outFileNamePrefix STAR_Align/"$I" \
#     --readFilesIn trimm_data/"$I"_1.par.fq.gz  trimm_data/"$I"_2.par.fq.gz

#     ## In featureCounts, paired-end reads must have -p option!!! 
#     ./featureCounts \
#     -p \
#     -t exon \
#     -g gene_id \
#     -F GTF \
#     -Q 32 \
#     -T 16 \
#     -a "$gtf" \
#     -o /gpfs/home/juagarcia/counts/"$I" \
#     STAR_Align/"$I"Aligned.out.bam 
# done

#STAR

## Check for installed modules 

 
##python3.6 fpkm.py -d "$1"/fpkm/

##python3.6 tpm_to_C.py -d "$1"/tpm/ -f "$2"""

cd "$1"; 

files=$(find Conc/ -type f | sort)
files=$(readlink -f $files)

count=0
touch Conc.csv

for J in $files 
do 
    count=$((count+1))
    iden=$(echo "$J" | rev | cut -d '/' -f 2 | rev | cut -d '_' -f 1)
    cat "$J" > temp 
    sed "1d" temp > tempfile ; mv tempfile temp
    sed "1 i\Gene,$iden" temp > tempfile ; mv tempfile temp 
    if (($count == 1))
    then 
        cat temp > Conc.csv
    else 
        cut -d ',' -f 2 temp  > count_ind
        paste -d ',' Conc.csv count_ind > temp2 && mv temp2 Conc.csv 
        rm -f temp2 
    fi 
done 

rm count_ind temp 