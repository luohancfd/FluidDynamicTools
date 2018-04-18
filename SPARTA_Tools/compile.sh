#!/bin/bash -e
module load intel
ROOT_DIR=$HOME/Soft/DSMC
mkdir -p $ROOT_DIR
GREP_OLD=$GREP_OPTIONS
export GREP_OPTIONS='--color=auto'
## Get SPARTA
cd $ROOT_DIR
#wget -N http://www.sandia.gov/~sjplimp/tars/sparta.tar.gz
#tar -zxf sparta.tar.gz
SRC_DIR=$(find . -maxdepth 1 -type d -name "sparta*" )
SRC_DIR="$ROOT_DIR${SRC_DIR:1}"
ln -sf $SRC_DIR $HOME/Soft/DSMC/sparta 

## Get PNG and compile
cd $SRC_DIR
mkdir -p lib
LIB_DIR=$SRC_DIR/lib
cd lib
mkdir -p src
cd src
#wget -N http://prdownloads.sourceforge.net/libpng/libpng-1.6.34.tar.gz?download -O libpng.tar.gz
#tar -zxf libpng.tar.gz
PNG_SRC=$(find . -maxdepth 1 -type d -name "libpng*")
PNG_SRC="$(pwd)${PNG_SRC:1}"
cd $PNG_SRC
#./configure --prefix=$LIB_DIR
make
make install 

## create a new makefile for impi
#cp $SRC_DIR/src/MAKE/Makefile.mpi $SRC_DIR/src/MAKE/Makefile.impi
#cd $SRC_DIR/src/MAKE
#sed -i 's/mpic++/mpiicpc/g' Makefile.impi
#sed -i '/^SPARTA_INC / s/$/ -DSPARTA_PNG/' Makefile.impi
#sed -i "s|^JPG_INC.*|JPG_INC = -I$LIB_DIR\/include|" Makefile.impi
#sed -i "s|^JPG_LIB.*|JPG_LIB = -L$LIB_DIR\/lib -lpng|" Makefile.impi

## Compile
module load intel
cd $SRC_DIR/src
make impi
export GREP_OPTIONS=$GREP_OLD

## link 
ln -sf $SRC_DIR/src/spa_impi $SRC_DIR/spa

### modulefile
#VER=${SRC_DIR#*-}

#if [[ $HOST_NAME == *"ecn"* ]]; then
    #mkdir -p $HOME/.local/modulefiles/sparta
    #modfile=$HOME/.local/modulefiles/sparta/$VER.lua
#else
    #mkdir -p $HOME/privatemodules/sparta
    #modfile=$HOME/privatemodules/sparta/$VER.lua
#fi 
  

#cat << EOF > $modfile
#whatis("Version: ${VER}")
#whatis("Description: Sets the environment for running Sparta")
#always_load("intel")
#prepend_path ('LD_LIBRARY_PATH', "${LIB_DIR}/lib") 
#prepend_path ('PATH', "$HOME/Soft/DSMC/sparta") 
#EOF



