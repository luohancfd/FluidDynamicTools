#!/bin/bash -e
# 
# syntax: ./run_ireac_ds1v.sh (no argument is needed)
# comments: set ireac=2 in DS1_openmp.f90
#
# --------------------------------------------------------------

# set machine (1 local run; 2 for job submission at carter)
NMAC=1    

# set run parameters
NPROC=4     #number of processors

# --------------------------------------------------------------

# set run directories
RUN_DIR=`pwd`
CASE=`basename ${RUN_DIR}`

# --------------------------------------------------------------

if [ ${NMAC} == 2 ]; then
  # load modules
  module load gcc
fi
echo "Setting script: done"

# --------------------------------------------------------------

# compile executable
gfortran -O3 -o ds1v_${CASE}.exe DS1_openmp_reac.f90
#ifort -O3 -fopenmp -o ds1v_${CASE}.exe DS1_openmp_reac.f90
echo "Compiling executable: done"

# --------------------------------------------------------------

# do for all temperatures
i=0
#TEMPS="3000 5000 7500 10000 12500 15000 17500 20000" #"1000 2000 3000 4000 5000 7500 10000 15000 20000"
TEMPS="8000 10000 15000 20000"
VTEMP="1000 3000 5000 7500 10000 12500 15000 17500 20000"
for TEMP in ${TEMPS}; do
i=$((i+1))
ii=$(printf "%03d" $i)

# create directory
cd ${RUN_DIR}
mkdir -p ${RUN_DIR}/"T="${TEMP}

j=0
for VT in ${VTEMP}; do
j=$((j+1))
jj=$(printf "%03d" $j)
mkdir -p ${RUN_DIR}/"T="${TEMP}/"Tv="${VT}
#echo $(pwd)
cp ${RUN_DIR}/ds1v_${CASE}.exe ${RUN_DIR}/"T="${TEMP}/"Tv="${VT}/ds1v_neq_${CASE}_${TEMP}_${VT}.exe
cp ${RUN_DIR}/*.mod ${RUN_DIR}/DS1VD.in ${RUN_DIR}/"T="${TEMP}/"Tv="${VT}

# create input file
cat << EOF > ${RUN_DIR}/"T="${TEMP}/"Tv="${VT}/input.txt
3       !new run
0       !seed
0       !continue with current cells
${TEMP} !translational temperature
${VT} !vibrational temperature
EOF

# local run
cd ${RUN_DIR}/"T="${TEMP}/"Tv="${VT}
ulimit -s unlimited
export OMP_NUM_THREADS=${NPROC}
if [[ "$TEMP" == "$VT" ]]; then
  echo "skip T="${VT}
else
  echo "run T="${TEMP}" VT="${VT}
  nohup ./ds1v_neq_${CASE}_${TEMP}_${VT}.exe <input.txt > ds1v.log 2>&1  &
fi

done
done
# --------------------------------------------------------------

# delete exec and module files
cd ${RUN_DIR}
#rm  *.mod *.exe
#echo "Deleting exec and module files: done"

# --------------------------------------------------------------

exit

