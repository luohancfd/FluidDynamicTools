#!/bin/bash

COMP=$2
FILE=$1
WORK_DIR=`pwd`
EXENAME="ds1v_"$(basename $WORK_DIR)".exe"

if [ "$COMP" = "ifort" ]; then
  ifort -O3 -init=zero -qopenmp -o $EXENAME $WORK_DIR/$FILE
else
  gfortran -O3 -fopenmp -o $EXENAME $WORK_DIR/$FILE
fi