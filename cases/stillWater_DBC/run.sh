#!/bin/bash

#========================================================
#source code path & case path
src_path='/home/tomas/Documents/__sovler/tinySPH_double_mr/src'
case_path='/home/tomas/Documents/__sovler/tinySPH_double_mr/cases/stillWater_DBC/'

#simulation file
simulation='SPH_simulationSymplectic.h'

#========================================================
#Create output directories if necessary
mkdir OUTPUT OUTPUT/resultsAll OUTPUT/resultsFluidOnly OUTPUT/resultsInterpolated

#Create temporary run folder
rm -r case_compiled
mkdir case_compiled

#Compile and rund case
cd ${src_path}
sed -i '2 s|.*|#include "'${case_path}${simulation}'"|' main.cu
make
cp main ${case_path}'case_compiled'
cd ${case_path}'case_compiled'
./main


#Finish run, clear temp. files
#rm -r case_compiled
#========================================================
