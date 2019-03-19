#!/bin/bash -e
# 
# syntax: ./run_restart_shock.sh (no argument is required)
#
# --------------------------------------------------------------

# set machine (1 local run; 2 for job submission at carter)
NMAC=1

# set run parameters
NPROC=10          #number of processors
OLDICOPY=$1      #unsteady shock run which reach the center of domain
ICOPY=$(($OLDICOPY+1))       #initial run/copy
FCOPY=${ICOPY}   #final

# --------------------------------------------------------------

# set run directories
RUN_DIR=`pwd`
CASE=`basename ${RUN_DIR}`
NEW_DIR=${RUN_DIR}/${CASE}_${ICOPY}
OLD_DIR=${RUN_DIR}/${CASE}_${OLDICOPY}

if [ ! -d $OLD_DIR ]; then
  echo "$OLD_DIR doesn't exist"
  exit
fi

if [ -d $NEW_DIR ]; then
  echo "$NEW_DIR shouldn't exist"
  exit
else
  mkdir -p $NEW_DIR
fi

# --------------------------------------------------------------

# create dir
rsync -av --include="*.f90" --include="*.exe" --include="RESTART.DAT" \
          --include="PARAMETERS.DAT" --include="*.mod" \
          --exclude="*" $OLD_DIR/ $NEW_DIR/
cp ${RUN_DIR}/DS1VD.in ${NEW_DIR}/
OLD_EXE=$(find ${NEW_DIR} -name "*.exe" | head -1)
OLD_EXE=$(basename ${OLD_EXE})
NEW_EXE=$(echo ${OLD_EXE} | sed "s|${OLDICOPY}|${ICOPY}|")
mv ${NEW_DIR}/${OLD_EXE} ${NEW_DIR}/${NEW_EXE}
echo "${CASE}_${ICOPY}: unsteady run to stablelize the shock" >> ${RUN_DIR}/run.log
echo "${CASE}_${ICOPY}: ${ICOPY}" >> ${RUN_DIR}/seed.log

# create input file with different seeds
cat <<EOF> ${NEW_DIR}/input2shock.txt
2      !restart
$ICOPY     !seed
0      !continue with current cells
1      !iretrem
1      !full removal
-0.005 !xrem
EOF

if [ ${NMAC} == 2 ]; then

# creat submission job
cat <<EOF> ${RUN_DIR}/${CASE}_${ICOPY}/${CASE}_${ICOPY}.sub
#!/bin/bash
#PBS -l nodes=1:ppn=${NPROC}
#PBS -l walltime=360:00:00
#PBS -q alexeenk
#PBS -V
#PBS -N ${CASE}_${ICOPY}
#PBS -j oe

module load gcc
cd ${RUN_DIR}/${CASE}_${ICOPY}
export OMP_NUM_THREADS=${NPROC}
./ds1v_${CASE}.exe <input2shock.txt
EOF

# submit job
chmod 755 ${RUN_DIR}/${CASE}_${ICOPY}/${CASE}_${ICOPY}.sub
qsub ${RUN_DIR}/${CASE}_${ICOPY}/${CASE}_${ICOPY}.sub
echo "Submitting job ${ICOPY}: done"

else

# local run
cd ${NEW_DIR}
ulimit -s unlimited
export OMP_NUM_THREADS=${NPROC}
nohup ./${NEW_EXE} <input2shock.txt >>${CASE}_${ICOPY}.log &

fi


# --------------------------------------------------------------

exit

