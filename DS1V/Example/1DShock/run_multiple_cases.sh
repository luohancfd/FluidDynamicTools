#!/bin/bash -e
# 
# syntax: ./run_multiple_cases.sh (no argument is required)
#
# --------------------------------------------------------------

# set machine (1 local run; 2 for job submission at carter)
NMAC=1

# set run parameters
NPROC=10          #number of processors
ICOPY=101        #initial run/copy
FCOPY=${ICOPY}   #final

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
#gfortran -O3 -fopenmp -o ds1v_${CASE}.exe DS1_openmp_202.f90
echo "Compiling executable: done"

# --------------------------------------------------------------

# create first copy
ii=$(printf "%03d" $ICOPY)
mkdir -p ${RUN_DIR}/${CASE}_${ii}
cp *.mod *.f90 ${RUN_DIR}/${CASE}_${ii}/
cp ds1v_${CASE}.exe ${RUN_DIR}/${CASE}_${ii}/ds1v_${CASE}_${ii}.exe
mv *.DAT *.TXT *.BIN ${RUN_DIR}/${CASE}_${ii}/
echo "Creating copy ${ii}: done"

# create remaining copies
cd ${RUN_DIR}
SCOPY=ICOPY+1
for (( i=${SCOPY}; i <= ${FCOPY}; i++ )); do
  iii=$(printf "%03d" $i)
  if [ -d "${RUN_DIR}/${CASE}_${iii}" ]; then
    cp *.mod  ${RUN_DIR}/${CASE}_${ii}/
    cp ${RUN_DIR}/ds1v_${CASE}.exe ${RUN_DIR}/${CASE}_${iii}/ds1v_${CASE}_${iii}.exe
    echo "Updating executable for copy ${iii}: done"
  else
    cp -r ${RUN_DIR}/${CASE}_${ii} ${RUN_DIR}/${CASE}_${iii}
    mv ${RUN_DIR}/${CASE}_${iii}/ds1v_${CASE}_${ii}.exe ${RUN_DIR}/${CASE}_${iii}/ds1v_${CASE}_${iii}.exe
    echo "Creating copy ${iii}: done"
  fi
done

# delete exec and module files
cd ${RUN_DIR}
#rm *.mod *.exe
echo "Deleting exec and module files: done"

# --------------------------------------------------------------

# do for all copies
for (( i=${ICOPY}; i <= ${FCOPY}; i++ )); do
ii=$(printf "%03d" $i)

# create input file with different seeds
cat <<EOF> ${RUN_DIR}/${CASE}_${ii}/input2.txt
2    !restart
$i   !seed
0    !continue with current cells
0    !no molecule removal
EOF

if [ ${NMAC} == 2 ]; then

# creat submission job
cat <<EOF> ${RUN_DIR}/${CASE}_${ii}/${CASE}_${ii}.sub
#!/bin/bash
#PBS -l nodes=1:ppn=${NPROC}
#PBS -l walltime=360:00:00
#PBS -q alexeenk
#PBS -V
#PBS -N ${CASE}_${ii}
#PBS -j oe

module load gcc
cd ${RUN_DIR}/${CASE}_${ii}
export OMP_NUM_THREADS=${NPROC}
./ds1v_${CASE}_${ii}.exe <input2.txt
EOF

# submit job
chmod 755 ${RUN_DIR}/${CASE}_${ii}/${CASE}_${ii}.sub
qsub ${RUN_DIR}/${CASE}_${ii}/${CASE}_${ii}.sub
echo "Submitting job ${ii}: done"

else

# local run
cd ${RUN_DIR}/${CASE}_${ii}
ulimit -s unlimited
export OMP_NUM_THREADS=${NPROC}
nohup ./ds1v_${CASE}_${ii}.exe < input2.txt > ${CASE}_${ii}.log 2>&1 &

fi

done

# --------------------------------------------------------------

exit

