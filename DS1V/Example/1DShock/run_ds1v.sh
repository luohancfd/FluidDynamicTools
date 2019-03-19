#!/bin/bash -e
# 
# syntax: ./run_ds1v.sh input3-10.txt
#
# --------------------------------------------------------------

# set machine (1 local run; 2 for job submission at carter)
NMAC=1    

# set run parameters
NPROC=20    #number of processors
INPUT=$1

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
cp DS1_openmp_reac.f90 DS1_shock_unsteady.f90
sed -i 's|^NSCELLS=[0-9]\+|NSCELLS=1|' DS1_shock_unsteady.f90
sed -i 's|^IREAC=[0-9]\+|IREAC=0|' DS1_shock_unsteady.f90
sed -i 's|^IRELAX=[0-9]\+|IRELAX=1|' DS1_shock_unsteady.f90
sed -i 's|^IMFS=[0-9]\+|IMFS=0|' DS1_shock_unsteady.f90
gfortran -O3 -fopenmp -o ds1v_${CASE}.exe DS1_shock_unsteady.f90
echo "Compiling executable: done"
# --------------------------------------------------------------

if [ ${NMAC} == 2 ]; then

# create submission job
cat <<EOF> ${RUN_DIR}/${CASE}.sub
#!/bin/bash
#PBS -l nodes=1:ppn=${NPROC}
#PBS -l walltime=180:00:00
#PBS -q alexeenk
#PBS -V
#PBS -N ${CASE}
#PBS -j oe

module load gcc
cd ${RUN_DIR}
export OMP_NUM_THREADS=${NPROC}
./ds1v_${CASE}.exe <${INPUT}
EOF
echo "Creating job: done"

# submit job
chmod 755 ${RUN_DIR}/${CASE}.sub
qsub ${RUN_DIR}/${CASE}.sub
echo "Submitting job: done"

else

# local run
cd ${RUN_DIR}
ulimit -s unlimited
export OMP_NUM_THREADS=${NPROC}
nohup ./ds1v_${CASE}.exe < ${INPUT} >> ${CASE}.log &

fi

# --------------------------------------------------------------

exit

