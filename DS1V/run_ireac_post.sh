#!/bin/bash -e
# 
# syntax: ./run_ireac_post.sh (no argument is needed)
#
# --------------------------------------------------------------

# set run directories
RUN_DIR=`pwd`
CASE=`basename ${RUN_DIR}`

# --------------------------------------------------------------

# delete old file
if [ -f "${RUN_DIR}/DS1REAC_${CASE}.DAT" ]; then
  rm ${RUN_DIR}/DS1REAC_${CASE}.DAT
  echo "Deleting old file: done"
fi

# --------------------------------------------------------------

# do for all temperatures
ICOPY=1
FCOPY=2
for (( i=${ICOPY}; i <= ${FCOPY}; i++ )); do
  ii=$(printf "%03d" $i)
  tail -n 1 ${RUN_DIR}/${CASE}_${ii}/DS1REAC.DAT >> ${RUN_DIR}/DS1REAC_${CASE}.DAT
done
echo "Creating new file: done"
echo " "

# --------------------------------------------------------------

more DS1REAC_${CASE}.DAT

exit

