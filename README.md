This repository contains all the scripts and data files used during the Masters final project of Juan Manuel Garcia. 

## Protein turnover GUI

The protein turnover model is found in the *model* folder. Files necessary to run the interface (*ui.r*, *server.r*, *global.r*, ...) are found in this folder. Users can run the interface by clicking on *run.bat* if they are working under a Windowd environment, or run the following command on a terminal under a UNIX environment after downloading the files: 
```bash
cd Stage_M2
Rscript ./runShinyApp.R
```
It is necessary to have *R* installed and the packages required to use the interface are automatically installed. 

In *data_kiwi*, data files for the kiwifruit specie and a script, *main_clust.r*, to run the model without any interface are found. Default values for weight (double sigmoid) and mRNA fitting (polynomial of third degree on log-transformed values) are set to run the model. This script is to be run on a cluster to make use of more ressources for the parallel computing implemented in the model. To run it, users have to run the *run_main.sh* file on a terminal: 
```bash
./run_main.sh -w WORKINGDIR -we WEIGHTFILE -co CONCFILE -ot OUTPUT
```
, where *WORKINGDIR* is the directory where *main_clust.r* is located (it is usually ths same directory where *run_main.sh* is). *WEIGHTFILE* is the weight data file, a *csv* file with 2 columns (time instances and weight values in gFW). *CONCFILE* is a *xlsx* file with transcript and protein concentration tables in different tabs. Name of tabs is set to be *Transcripts* and *Proteins* by default, resepctively. However, users can change the names in the script if needed. *OUTPUT* is the output filename where all the results are stored in a *csv* file.   

## RNA-Seq files 

In this folder all the required files to run the quantitative RNA-Seq analysis are found. 2 Python scripts, *fpkm_to_tpm.py* and *tpm_to_C.py* are found. These 2 scripts are integrated in the main bash file that is to be run, *quant_rna_seq.sh* and calculate tpm and concentration of transcripts from the RNA-Seq results. In order to calculate concentration values, a file containing information about concentartion of RNA spikes is necessary. This file is included in *dataKentaro*. *count_to_tpm.py* is specially bound for changes, since the reading of spike data depends exclusively on the file. *quant_rna_seq.sh* is written to be run in Genotoul cluster, under a SLURM working scheduler. Users might change the paths where raw transcriptomic files are stored. The line to be run on a terminal to run this script shoud be:  
```bash
./quant_rna_seq.sh DIRECTORY RESULTSFOLDER SPECIE 
```
, where *DIRECTORY* is the path to the directory where the python files are found,  *RESULTSFOLDER* denotes the new folder where *tpm* and concentration files will be saved. Finally, *SPECIE* denotes the name of the excel tab where spike information is found. This argument is specific to the data used for the kiwifruit and is subject to future changes.  

In addition, a simple bash file, *extract_fasta_samples.sh* is added. This scripts extracts a number of lines (*LINES*) given by the user from all the *fasta* or *fastq* files located in a specific directory (*DIRECTORY*).   
```bash
./extract_fasta_samples.sh LINES DIRECTORY 
```
## R files for statistic analyses

Scripts used in *R* for stastistical analyses for both transcriptomic and proteomic dara are found in *Stat_R_files*. *script.r* and *functions.r* contain all the steps followed to get valid mRNA concentration table as well as all the graphs. In addition, 2 Python files are included in this folder: *mapman_parser.py* and *gff_gene_trans_parser.py*. These files parse files created from MapMan and *gff* files, respectively. Both are built using *argparse* module in **Python3** and more information about their use can be found by running `python3 ./mapman_parser.py -h` and `python3 ./gff_gene_trans_parser.py -h`. 

In the *Proteo* subfolder, a script in R to do quantitative proteomic analysis (*Proteo_Juanma.R*) is found. This script is to be modified since filenames should be changed for the user's convenience. A *RData* file is generated where the protein concentration table is saved. This table is later used on *stats_proteo.r* to get the graphs as in *script.r* for the mRNA.   