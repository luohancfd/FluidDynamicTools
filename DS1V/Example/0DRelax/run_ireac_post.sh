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
if [ -f "${RUN_DIR}/RELAX_${CASE}.DAT" ]; then
  rm ${RUN_DIR}/RELAX_${CASE}.DAT
  echo "Deleting old file: done"
fi

# --------------------------------------------------------------

# do for all temperatures
ICOPY=0
FCOPY=10
for (( i=${ICOPY}; i <= ${FCOPY}; i++ )); do
  ii=$(printf "%03d" $i)
  cp ${RUN_DIR}/${CASE}_${ii}/RELAX.DAT ${RUN_DIR}/RELAX_${ii}.DAT
  echo "Creating new file: ${ii}"
done
echo " "

# --------------------------------------------------------------

#more DS1REAC_${CASE}.DAT

exit

