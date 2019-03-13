#!/bin/bash -e
# 
# syntax: ./run_ireac_ds1v.sh (no argument is needed)
# comments: set ireac=2 in DS1_openmp.f90
#
# --------------------------------------------------------------

# set machine (1 local run; 2 for job submission at carter)
NMAC=1    

# set run parameters
NPROC=2     #number of processors

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
gfortran -O3 -fopenmp -o ds1v_${CASE}.exe DS1_openmp_reac.f90
echo "Compiling executable: done"

# --------------------------------------------------------------

# do for all temperatures
i=-1
TEMPS="600 1000 1500 2000 4000 6000 8000 10000 15000 1750 20000"
#TEMPS="15000 20000"
for TEMP in ${TEMPS}; do
i=$((i+1))
ii=$(printf "%03d" $i)

# create directory
cd ${RUN_DIR}
if [ -d "${RUN_DIR}/${CASE}_${ii}" ]; then
  rm *.mod *.exe
  echo "Copy ${ii} already exists: stop script"
  exit
else
  mkdir ${RUN_DIR}/${CASE}_${ii}
  cp ds1v_${CASE}.exe ${RUN_DIR}/${CASE}_${ii}/ds1v_${CASE}_${ii}.exe
  cp *.mod DS1VD.in ${RUN_DIR}/${CASE}_${ii}
  echo "Creating copy ${ii}: done"
fi

# create input file
cat <<EOF> ${RUN_DIR}/${CASE}_${ii}/input.txt
3       !new run
0       !seed
0       !continue with current cells
${TEMP} !temperature
500     !vibrational temperature
EOF

if [ ${NMAC} == 2 ]; then

# create submission job
cat <<EOF> ${RUN_DIR}/${CASE}_${ii}/${CASE}_${ii}.sub
#!/bin/bash
#PBS -l nodes=1:ppn=${NPROC}
#PBS -l walltime=24:00:00
#PBS -q alexeenk
#PBS -V
#PBS -N ${CASE}_${ii}
#PBS -j oe

module load gcc
cd ${RUN_DIR}/${CASE}_${ii}
export OMP_NUM_THREADS=${NPROC}
./ds1v_${CASE}_${ii}.exe <input.txt
EOF

# submit job
chmod 755 ${RUN_DIR}/${CASE}_${ii}/${CASE}_${ii}.sub
qsub ${RUN_DIR}/${CASE}_${ii}/${CASE}_${ii}.sub
echo "Submitting job: done"

else

# local run
cd ${RUN_DIR}/${CASE}_${ii}
ulimit -s unlimited
export OMP_NUM_THREADS=${NPROC}
nohup ./ds1v_${CASE}_${ii}.exe <input.txt >>${CASE}_${ii}.log &

fi

done

# --------------------------------------------------------------

# delete exec and module files
cd ${RUN_DIR}
rm  *.mod *.exe
echo "Deleting exec and module files: done"

# --------------------------------------------------------------

exit
