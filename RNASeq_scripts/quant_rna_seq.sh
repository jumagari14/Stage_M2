#!/bin/bash

#############################
# les directives Slurm vont ici:

# Your job name (displayed by the queue)
#SBATCH -J Kiwi_red5

# walltime (hh:mm::ss)
#SBATCH -t 96:00:00

# Specify the number of nodes(nodes=) and the number of cores per nodes(tasks-pernode=) to be used
#SBATCH -N 1
#SBATCH --tasks-per-node=4


# change working directory
#SBATCH --chdir=.


#SBATCH  -p workq
#SBATCH --mem-per-cpu=40000
# fin des directives PBS
#############################


# useful informations to print
echo "#############################" 
echo "User:" $USER
echo "Date:" `date`
echo "Host:" `hostname`
echo "Directory:" `pwd`
echo "SLURM_JOBID:" $SLURM_JOBID
echo $MaxMemPerCPU
echo "SLURM_SUBMIT_DIR:" $SLURM_SUBMIT_DIR
echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST
echo "#############################" 

#############################

## Load necessary modules
module load system/parallel-20180122
module load system/Python-3.7.4
module load bioinfo/cufflinks-2.2.1
module load bioinfo/Trimmomatic-0.38 
module load bioinfo/HTSeq-0.9.1 
module load bioinfo/STAR-2.5.1b
module load bioinfo/FastQC_v0.11.7

cd ~/ 
list=$(ls work/data/*fastq* | xargs -n 1 basename | sed 's/\(.*\)_.*/\1/' | sort -u )
list1=$(echo $list | cut -f 14 -d " ")
list2=$(echo $list | cut -f 15-27 -d " ")
mkdir -p -m 755 bin/trimm_data 
mkdir -p -m 755 bin/trimm_data/quality

gff=$(find work/data/ -name "*gff[3]*")
fasta=$(find work/data/ -name "*.fa")
fake_fasta=$(find work/data/ -name "*.fasta")

if [ ! -d "work/bin/genome_ind" ] 
 then  ## Genome indexes are generated if necessary 
    mkdir -m 755 -p work/bin/genome_ind 
    STAR \
    --runThreadN 4 \
    --runMode genomeGenerate \
    --genomeDir /work/jgarcia/bin/genome_ind \
    --genomeFastaFiles "$fasta" "$fake_fasta" \
    --sjdbGTFtagExonParentGene ID \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbOverhang 299 \
    --genomeSAindexNbases 13 \
    --sjdbGTFfile "$gff"
fi  

mkdir -p -m 755 work/bin/STAR_Align  
mkdir -p -m 755 work/bin/counts 
mkdir -p -m 755 work/bin/fpkm 

for I in $list1
do 
    if [ ! -s "work/bin/trimm_data/quality/"$I"_R1_fastqc.html" ] && [ ! -s "work/bin/trimm_data/quality/"$I"_R2_fastqc.html" ]; then 
        fastqc work/data/"$I"_R1.fastq.gz work/data/"$I"_R2.fastq.gz --outdir=work/bin/trimm_data/quality
    fi 
    if [ ! -s "work/bin/trimm_data/"$I"_R1.par.fastq.gz" ] && [ ! -s "work/bin/trimm_data/"$I"_R2.par.fastq.gz" ]; then  
        java -jar $TRIM_HOME/trimmomatic.jar PE -threads 4 -phred33 work/data/"$I"_R1.fastq.gz work/data/"$I"_R2.fastq.gz work/bin/trimm_data/"$I"_R1.par.fastq.gz "$I"_R1.unp.fastq.gz work/bin/trimm_data/"$I"_R2.par.fastq.gz "$I"_R2.unp.fastq.gz ILLUMINACLIP:$TRIM_HOME/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
        rm -rf "$I"_R1.unp.fastq.gz "$I"_R2.unp.fastq.gz
    fi
    if [ ! -s "work/bin/STAR_Align/"$I"Aligned.sortedByCoord.out.bam" ]; then
        STAR --runThreadN 4 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --genomeDir work/bin/genome_ind --outFileNamePrefix work/bin/STAR_Align/"$I" --runMode alignReads --readFilesIn work/bin/trimm_data/"$I"_R1.par.fastq.gz  work/bin/trimm_data/"$I"_R2.par.fastq.gz  --outSAMstrandField intronMotif --outSAMattributes All 
    fi
    if [ ! -s "work/bin/counts/HTSeq_"$I".txt" ]; then
        htseq-count -f bam -t exon --stranded=no -i Parent work/bin/STAR_Align/"$I"Aligned.sortedByCoord.out.bam "$gff" > work/bin/counts/HTSeq_"$I".txt
    fi
    if [ ! -s "work/bin/fpkm/"$I"/genes.fpkm_tracking" ]; then 
        cufflinks -G "$gff" -o work/bin/fpkm/"$I"/ work/bin/STAR_Align/"$I"Aligned.sortedByCoord.out.bam
    fi 
done 
STAR --genomeLoad Remove --genomeDir work/bin/genome_ind

cd $1
cp work/bin/fpkm $1

python3.6 count_to_tpm.py -d "$2"/fpkm/

python3.6 tpm_to_C.py -d "$2"/tpm/ -f "$3"

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
    if (($count == 1)); then 
        cat temp > Conc.csv
    else 
        cut -d ',' -f 2 temp  > count_ind
        paste -d ',' Conc.csv count_ind > temp2 && mv temp2 Conc.csv 
        rm -f temp2 
    fi 
done 

rm count_ind temp 

exit 0 