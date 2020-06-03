#!/bin/sh

module load compiler/gcc-7.2.0 libraries/gdal-2.4.0_gcc-7.2.0 libraries/proj-5.2.0_gcc-7.2.0 libraries/geos-3.4.2_gcc-7.2.0 ; module load system/R-3.6.2_gcc-7.2.0

Rscript $1/main.r -o $1