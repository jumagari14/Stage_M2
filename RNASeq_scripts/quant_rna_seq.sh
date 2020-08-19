#!/bin/bash

#############################
# les directives Slurm vont ici:

# Your job name (displayed by the queue)
#SBATCH -J Fmol

# walltime (hh:mm::ss)
#SBATCH -t 10:00:00

# Specify the number of nodes(nodes=) and the number of cores per nodes(tasks-pernode=) to be used
#SBATCH -N 1
#SBATCH --tasks-per-node=8

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
module load bioinfo/cufflinks-2.2.1: /usr/local/bioinfo/src/Cufflink/cufflinks-2.2.1.Linux_x86_64/cufflinks
module load bioinfo/Trimmomatic-0.38 : java -jar $TRIM_HOME/trimmomatic.jar
module load bioinfo/HTSeq-0.9.1: 
module load bioinfo/STAR-2.6.0c
module load bioinfo/FastQC_v0.11.7

check=$(python3.7 -c "import HTSeq" | echo $?)
if (($check==1))
    then
    pip3.7 install --user HTSeq 
fi

export LD_LIBRARY_PATH=/gpfs/softs/contrib/apps/gcc/7.3.0/lib64/:/gpfs/softs/contrib/apps/gcc/7.3.0/lib
export LD_LIBRARY_PATH=/gpfs/softs/contrib/apps/python/3.7.2/lib:$LD_LIBRARY_PATH

cd /gpfs/home/juagarcia/ 
list=$(ls data/*fq | xargs -n 1 basename | sed 's/\(.*\)_.*/\1/' | sort -u)

mkdir -p -m 755 bin/trimm_data 
mkdir -p -m 755 bin/trimm_data/quality

gff=$(find data/ -not -path '*/\R*' -name "*gff*")
fasta=$(find data/ -name "*.fa")
~/scripts/gffread "$gff" -o "$gff".gtf
tail -n +4 "$gff".gtf > ~/data/tmp ; mv ~/data/tmp "$gff".gtf

if [ ! -d "bin/genome_ind" ] 
 then  ## Genome indexes are generated if necessary 
    mkdir -m 755 -p bin/genome_ind 
    ./scripts/STAR \
    --runThreadN 64 \
    --runMode genomeGenerate \
    --genomeDir ./bin/genome_ind \
    --genomeFastaFiles "$fasta" \
    --genomeSAindexNbases 14 \
    --sjdbOverhang 149 \
    --sjdbGTFfile "$gff".gtf \
    --genomeChrBinNbits 18 
fi  

mkdir -p -m 755 bin/STAR_Align  
mkdir -p -m 755 bin/counts


./scripts/STAR --genomeLoad LoadAndExit --genomeDir ./bin/genome_ind
for I in $list
do 
    ./scripts/FastQC/fastqc ./data/"$I"_R1.fq ./data/"$I"_R2.fq --outdir=./bin/trimm_data/quality
    java -jar ./scripts/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 64 -phred33 ./data/"$I"_R1.fq ./data/"$I"_R2.fq ./bin/trimm_data/"$I"_R1.par.fq "$I"_R1.unp.fq ./bin/trimm_data/"$I"_R2.par.fq "$I"_R2.unp.fq ILLUMINACLIP:./scripts/Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
    rm -rf "$I"_R1.unp.fq "$I"_R2.unp.fq

    ./scripts/STAR --runThreadN 64 --genomeDir ./bin/genome_ind --outFileNamePrefix ./bin/STAR_Align/"$I" --runMode alignReads --genomeLoad LoadAndKeep --readFilesIn ./bin/trimm_data/"$I"_R1.par.fq  ./bin/trimm_data/"$I"_R2.par.fq  --outSAMtype BAM SortedByCoordinate --twopassMode None --quantMode - --outSAMstrandField intronMotif --outSAMattrIHstart 1 --outSAMattributes NH HI AS nM NM MD MC jM jI ch --outSAMprimaryFlag OneBestScore --outSAMmapqUnique 60 --outSAMunmapped Within --outFilterType Normal --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 10 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.3 --outFilterMismatchNoverReadLmax 1.0 --outFilterScoreMin 0 --outFilterScoreMinOverLread 0.66 --outFilterMatchNmin 0 --outFilterMatchNminOverLread 0.66 --outSAMmultNmax -1 --outSAMtlen 1 --outBAMsortingThreadN 10 --outBAMsortingBinsN 50 --limitBAMsortRAM 10000200000
    python3.7 ./scripts/htseq-count -f bam -t mRNA --stranded=no -i geneID ./bin/STAR_Align/"$I"Aligned.sortedByCoord.out.bam "$gff".gtf > ./bin/counts/HTSeq_"$I".txt
done 
./scripts/STAR --genomeLoad Remove --genomeDir ./bin/genome_ind

python3.6 fpkm.py -d "$1"/fpkm/

python3.6 tpm_to_C.py -d "$1"/tpm/ -f "$2"

cd "$1"

files=$(find tpm/ -type f | sort)
files=$(readlink -f $files)
count=0
touch TPM_"$2".csv

for J in $files 
do 
    count=$((count+1))
    iden=$(echo "$J" | rev | cut -d '/' -f 2 | rev | cut -d '_' -f 1)
    cat "$J" > temp 
    sed "1d" temp > tempfile ; mv tempfile temp
    sed "1 i\Gene,$iden" temp > tempfile ; mv tempfile temp 
    if (($count == 1))
    then 
        cat temp > TPM_"$2".csv
    else 
        cut -d ',' -f 2 temp  > count_ind
        paste -d ',' TPM_"$2".csv count_ind > temp2 && mv temp2 TPM_"$2".csv 
        rm -f temp2 
    fi 
done 

rm count_ind temp

files=$(find Conc/ -type f | sort)
files=$(readlink -f $files)

count=0
touch Conc_"$2".csv

for J in $files 
do 
    count=$((count+1))
    iden=$(echo "$J" | rev | cut -d '/' -f 2 | rev | cut -d '_' -f 1)
    cat "$J" > temp 
    sed "1d" temp > tempfile ; mv tempfile temp
    sed "1 i\Gene,$iden" temp > tempfile ; mv tempfile temp 
    if (($count == 1))
    then 
        cat temp > Conc_"$2".csv
    else 
        cut -d ',' -f 2 temp  > count_ind
        paste -d ',' Conc_"$2".csv count_ind > temp2 && mv temp2 Conc_"$2".csv 
        rm -f temp2 
    fi 
done 

rm count_ind temp 


files=$(find readCount/ -type f | sort)
files=$(readlink -f $files)

count=0
touch readCount.csv

for J in $files 
do 
    count=$((count+1))
    iden=$(echo "$J" | rev | cut -d '/' -f 1 | rev | cut -d '_' -f 1 | cut -d '-' -f 2)
    cat "$J" > temp 
    sed $"1 i\Read\t$iden" temp > tempfile ; mv tempfile temp 
    if (($count == 1))
    then 
        paste -d ',' readCount.csv temp > file1 && mv file1 readCount.csv
    else 
        cut -f 2 temp  > count_ind
        paste -d ',' readCount.csv count_ind > temp2 && mv temp2 readCount.csv 
        rm -f temp2 
    fi 
done 

rm count_ind temp 

exit 0 