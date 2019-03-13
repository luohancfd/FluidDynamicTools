#!/bin/bash -e
#
# syntax: ./run_ireac_ds1v.sh (no argument is needed)
# comments: set ireac=2 in DS1_openmp.f90
# Modify IRM in DS1VD.in before run
# Modify LS+MS in code before run
# --------------------------------------------------------------
# set machine (1 local run; 2 for job submission at carter)
NMAC=1

# set run parameters
NPROC=2    #number of processors for one case
MNPROC=20   #total number of processors
# --------------------------------------------------------------
# set run directories and running range
RUN_DIR=`pwd`
DS1V_DIR=`pwd`
CASE=`basename ${RUN_DIR}`
TEMPS="3000 4000 5000 7500 10000 12500 15000 17500 20000"
#TEMPS="4000 5000 7500 10000 12500 1500 17500"
#TEMPS="17500"
# --------------------------------------------------------------
#module load gcc
if [ ${NMAC} == 2 ]; then
  # create header of submission job
  cat <<EOF> ${RUN_DIR}/${CASE}.sub
#!/bin/bash
#PBS -l nodes=1:ppn=${MNPROC}
#PBS -l walltime=12:00:00
#PBS -q alexeenk
#PBS -V
#PBS -N ${CASE}

module load gcc
cd ${RUN_DIR}
export OMP_NUM_THREADS=${NPROC}
pids=""
EOF
else
  export OMP_NUM_THREADS=${NPROC}
fi
# --------------------------------------------------------------
cd ${RUN_DIR}
# compile executable
if [ ! -f DS1_openmp_reac.f90 ]; then
  echo "Copy code file from" ${DS1V_DIR}
  cp ${DS1V_DIR}/DS1_openmp_reac.f90 ${RUN_DIR}/
fi
if [ ! -f ds1v_${CASE}.exe ]; then
  echo "Compiling executable: start"
  gfortran -O3 -fopenmp -o ds1v_${CASE}.exe DS1_openmp_reac.f90
  echo "Compiling executable: done"
fi
if [ ! -f ./DS1VD.in ]; then
  echo "DS1VD.in not found"
  return
fi
# --------------------------------------------------------------
# do for all temperatures
i=0
for TEMP in ${TEMPS}; do
  i=$((i+1))
  ii=$(printf "%03d" $i)

  # create directory
  cd ${RUN_DIR}
  CASE_DIR=${RUN_DIR}/"T="${TEMP}
  mkdir ${CASE_DIR}
  # ln -sf ${DS1V_DIR}/GAS_INFO/n2.bin ${CASE_DIR}/n2.bin
  # ln -sf ${DS1V_DIR}/GAS_INFO/no.bin ${CASE_DIR}/no.bin
  # ln -sf ${DS1V_DIR}/GAS_INFO/oo.bin ${CASE_DIR}/oo.bin
  ln -sf ${RUN_DIR}/ds1v_${CASE}.exe ${CASE_DIR}/ds1v_${CASE}_T_${TEMP}.exe
  cp ${RUN_DIR}/*.mod ${RUN_DIR}/DS1VD.in ${CASE_DIR}/

  # create input file
  cat << EOF > ${CASE_DIR}/input.txt
3       !new run
0       !seed
0       !continue with current cells
${TEMP} !translational temperature
${TEMP} !vibrational temperature
EOF

  if [ ${NMAC} == 2 ]; then
  # append to end of sub script for NMAC = 2
    cat <<EOF>> ${RUN_DIR}/${CASE}.sub
cd ${CASE_DIR}
nohup ./ds1v_${CASE}_T_${TEMP}.exe < input.txt > ${CASE}_${ii}.log 2>&1 &
pids="\$pids \$!"
EOF
  else
    # local run
    cd ${CASE_DIR}
    ulimit -s unlimited
    echo "Start running case: "$(pwd)
    nohup $(pwd)/ds1v_${CASE}_T_${TEMP}.exe < input.txt > ${CASE}_${ii}.log 2>&1 &
  fi
done

# submit job
if [ ${NMAC} == 2 ]; then
  cat <<EOF>> ${RUN_DIR}/${CASE}.sub
wait \$pids
EOF
  chmod 755 ${RUN_DIR}/${CASE}.sub
 # qsub ${RUN_DIR}/${CASE}.sub
  echo "Submitting job manually"
fi

# --------------------------------------------------------------

exit

