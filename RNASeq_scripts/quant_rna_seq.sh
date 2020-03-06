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
# module load boost/1.69.0
# module load samtools/1.9
# module load htslib/1.9

check=$(python3.7 -c "import HTSeq" | echo $?)
if (($check==1))
    then
    pip3.7 install --user HTSeq 
fi

export LD_LIBRARY_PATH=/gpfs/softs/contrib/apps/gcc/7.3.0/lib64/:/gpfs/softs/contrib/apps/gcc/7.3.0/lib
export LD_LIBRARY_PATH=/gpfs/softs/contrib/apps/python/3.7.2/bin/python3.7/lib/libpython3.7m.so:$LD_LIBRARY_PATH

cd /gpfs/home/juagarcia/ 
list=$(ls data/*fq | xargs -n 1 basename | sed 's/\(.*\)_.*/\1/' | sort -u)

mkdir -p -m 755 bin/trimm_data 
mkdir -p -m 755 bin/trimm_data/quality

#gtf=$(find data/ -not -path '*/\R*' -name "*.gtf")
gff=$(find data/ -not -path '*/\R*' -name "*gff*")
fasta=$(find data/ -name "*.fa")

if [ ! -d "bin/genome_ind" ] 
 then  ## Genome indexes are generated if necessary 
    mkdir -m 755 -p bin/genome_ind 
    ./scripts/STAR \
    --runThreadN 32 \
    --runMode genomeGenerate \
    --genomeDir ./bin/genome_ind \
    --genomeFastaFiles "$fasta" \
    --genomeSAindexNbases 12 \
    --sjdbOverhang 149 \
    --sjdbGTFfile "$gff" \
    --genomeChrBinNbits 12 
fi  

mkdir -p -m 755 bin/STAR_Align  
mkdir -p -m 755 bin/counts


# cp subdata/"$I"_R1.fq.gz subdata/"$I"_R1.fastqsanger 
# cp subdata/"$I"_R2.fq.gz subdata/"$I"_R2.fastqsanger

parallel -j 3 "./scripts/FastQC/fastqc ./data/{}_R1.fq ./data/{}_R2.fq --outdir=./bin/trimm_data/quality" ::: $list
parallel -j 3 "java -jar ./scripts/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 32 -phred33 ./data/{}_R1.fq ./data/{}_R2.fq ./bin/trimm_data/{}_R1.par.fq {}_R1.unp.fq ./bin/trimm_data/{}_R2.par.fq {}_R2.unp.fq ILLUMINACLIP:./scripts/Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25" ::: $list

parallel -j 3 "rm -rf {}_R1.unp.fq {}_R2.unp.fq" ::: $list

./scripts/STAR --genomeLoad LoadAndExit --genomeDir ./bin/genome_ind
# Load genome just once to save RAM memory 
parallel --compress -j 3 "./scripts/STAR --runThreadN 32 --genomeDir ./bin/genome_ind --outFileNamePrefix ./bin/STAR_Align/{} --runMode alignReads --genomeLoad LoadAndKeep --readFilesIn ./bin/trimm_data/{}_R1.par.fq  ./bin/trimm_data/{}_R2.par.fq  --outSAMtype BAM SortedByCoordinate --twopassMode None --quantMode - --outSAMstrandField intronMotif --outSAMattrIHstart 1 --outSAMattributes NH HI AS nM NM MD MC jM jI ch --outSAMprimaryFlag OneBestScore --outSAMmapqUnique 60 --outSAMunmapped Within --outFilterType Normal --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 10 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.3 --outFilterMismatchNoverReadLmax 1.0 --outFilterScoreMin 0 --outFilterScoreMinOverLread 0.66 --outFilterMatchNmin 0 --outFilterMatchNminOverLread 0.66 --outSAMmultNmax -1 --outSAMtlen 1 --outBAMsortingThreadN 10 --outBAMsortingBinsN 50 --limitBAMsortRAM 10000200000" ::: $list
./scripts/STAR --genomeLoad Remove --genomeDir ./bin/genome_ind 


parallel -j 3 "python3.7 ./scripts/htseq-count --stranded=no STAR_Align/{}Aligned.sortedByCoord.out.bam "$gff"" ::: $list

#parallel -j 3 "cufflinks -G "$gff" STAR_Align/{}Aligned.sortedByCoord.out.bam" ::: $list

#STAR

## Check for installed modules 

 
##python3.6 fpkm.py -d "$1"/fpkm/

##python3.6 tpm_to_C.py -d "$1"/tpm/ -f "$2"

# cd "$1"; 

# files=$(find Conc/ -type f | sort)
# files=$(readlink -f $files)

# count=0
# touch Conc.csv

# for J in $files 
# do 
#     count=$((count+1))
#     iden=$(echo "$J" | rev | cut -d '/' -f 2 | rev | cut -d '_' -f 1)
#     cat "$J" > temp 
#     sed "1d" temp > tempfile ; mv tempfile temp
#     sed "1 i\Gene,$iden" temp > tempfile ; mv tempfile temp 
#     if (($count == 1))
#     then 
#         cat temp > Conc.csv
#     else 
#         cut -d ',' -f 2 temp  > count_ind
#         paste -d ',' Conc.csv count_ind > temp2 && mv temp2 Conc.csv 
#         rm -f temp2 
#     fi 
# done 

# rm count_ind temp 

# exit 0 