#!/bin/sh
#SBATCH -mmem-per-cpu=5000
#SBATCH -c 20
#SBATCH -t 96:00:00
#SBATCH -J Turnover_model
module load compiler/gcc-7.2.0 libraries/gdal-2.4.0_gcc-7.2.0 libraries/proj-5.2.0_gcc-7.2.0 libraries/geos-3.4.2_gcc-7.2.0 ; module load system/R-3.6.2_gcc-7.2.0

# module load R/3.6.3
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -w|--workingDir)
    WORK="$2"
    shift # past argument
    shift # past value
    ;;
    -we|--fileweight)
    WEIGHT_F="$2"
    shift # past argument
    shift # past value
    ;;
    -co|--fileConcentration)
    MAINFILE="$2"
    shift # past argument
    shift # past value
    ;;
    -ot|--outputFile)
    OUT="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters
 echo ${OUT}

Rscript $WORK/main_clust.r -o $WORK -f $MAINFILE -w $WEIGHT_F -a $OUT -n $SLURM_CPUS_PER_TASK
