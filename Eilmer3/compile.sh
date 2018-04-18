#!/bin/bash -e
compile_popt () {
   OLD_DIR=$(pwd)
   cd ${ROOT_DIR}
   wget -N http://rpm5.org/files/popt/popt-1.16.tar.gz
   mkdir -p ${ROOT_DIR}/Addons
   mkdir -p ${POPT_SRC}
   tar -zxf popt-1.16.tar.gz --strip-components 1 -C ${POPT_SRC}
   cd ${POPT_SRC}
   ./configure  --prefix=${POPT_ROOT} --disable-static
   make install
   cd $OLD_DIR
}

install_cea () {
   OLD_DIR=$(pwd)
   cd ${ROOT_DIR}
   wget -N http://www.grc.nasa.gov/WWW/CEAWeb/CEA+Fortran.tar.Z 
   mkdir -p ${CFCFD3_SRC}/extern/cea2
   tar -zxf CEA+Fortran.tar.Z --strip-components 1 -C ${CFCFD3_SRC}/extern/cea2
   cd $OLD_DIR
}

update_eilmer3 () {
   OLD_DIR=$(pwd)
   mkdir -p $CFCFD3_SRC
   #update files
   if [[ $HOST_NAME == *"ecn"* ]]; then
      cd $GIT_DIR
      git checkout $GIT_BRANCH 
      if [ "$1" == "install" ]; then
         rsync -av --exclude='.git*' --exclude='.hg*' $GIT_DIR/ $CFCFD3_SRC/
      else
         PATCH_FILE=$(git format-patch -n HEAD^)
         PATCH_FILE=$GIT_DIR"/"$PATCH_FILE
         cd $CFCFD3_SRC
         patch -p1 --verbose --dry-run < $PATCH_FILE
         rm $PATCH_FILE
      fi
  else
      if [ "$1" == "install" ]; then
         cd $CFCFD3_SRC/../
         git clone -b dev --depth 1 git@github.com:luohancfd/CFCFD-NG.git CFCFD3
      else
         cd $CFCFD3_SRC
         git fetch --all
         git reset --hard origin/dev
      fi
  fi
  cd $OLD_DIR
}

build_eilmer3 () {
  OLD_DIR=$(pwd)
  cd ${EILMER3_BUILD}
  sed -i "s|LLIB +=.*|LLIB += -lpopt -ldl -L\/opt\/local\/lib -L${POPT_ROOT}\/lib|g" "makefile"
  sed -i "s/gfortran -m32 -std=legacy\(.*cea2.f\)/gfortran -std=legacy\1/g" ${CFCFD3_SRC}/lib/radiation/build/makefile
  sed -i "s/gfortran -m32 -std=legacy\(.*cea2.f\)/gfortran -std=legacy\1/g" ${CFCFD3_SRC}/lib/gas/build/makefile
  make TARGET=for_intel_mpi INSTALL_DIR=${EILMER3_DIR} install

  cd ${POSHAX_BUILD}
  make TARGET=for_intel INSTALL_DIR=${EILMER3_DIR} install
  cd $OLD_DIR
}

install_mod () {
  OLD_DIR=$(pwd)
  if [[ $HOST_NAME == *"ecn"* ]]; then
     modfile=$HOME/.local/modulefiles/cfcfd3/3.lua
  else
     modfile=$HOME/privatemodules/cfcfd3/3.lua
  fi 
  mkdir -p $(dirname $modfile)
cat << EOF > $modfile
whatis("Version: 3")
whatis("Description: Sets the environment for running CFCFD3")
always_load("intel")
always_load("$PYTHON_MOD")
prepend_path ('PATH', "${EILMER3_DIR}") 
prepend_path ('LD_LIBRARY_PATH', "${EILMER3_DIR}") 
prepend_path ('PYTHONPATH',"${EILMER3_DIR}") 
prepend_path ('LD_LIBRARY_PATH', "${POPT_ROOT}") 
setenv("LUA_PATH","${EILMER3_DIR}/?.lua")
setenv("LUA_CPATH","${EILMER3_DIR}/?.so")
setenv("E3BIN","${EILMER3_DIR}")
EOF
  cd $OLD_DIR
}


   
GLIBCV=$(ldd --version)
OVERWRITE=1
OVERWRITEBIN=1
module load intel
HOST_NAME=$(hostname)
ROOT_DIR=${HOME}/Soft/CFD/Eilmer3
if [[ $HOST_NAME == *"ecn"* ]]; then
  PYTHON_MOD="python/py27"
  module load $PYTHON_MOD
else
  PYTHON_MOD="python"
  module load $PYTHON_MOD
fi
GIT_DIR=XXXX
GIT_BRANCH=dev
CGNS_SRC=${ROOT_DIR}/Addons/cgns
POPT_SRC=${ROOT_DIR}/Addons/popt
POPT_ROOT=${ROOT_DIR}/lib
CFCFD3_SRC=${ROOT_DIR}/CFCFD3
EILMER3_BUILD=${CFCFD3_SRC}/app/eilmer3/build
EILMER3_DIR=${CFCFD3_SRC}/e3bin
POSHAX_BUILD=${CFCFD3_SRC}/app/poshax3/build

if [ -z "$CPATH" ]; then  
   export CPATH=${POPT_ROOT}/include
else
   export CPATH=${POPT_ROOT}/include:$CPATH
fi

if [ -z "$LD_LIBRARY_PATH" ]; then  
   export LD_LIBRARY_PATH=${POPT_ROOT}/lib
else
   export LD_LIBRARY_PATH=${POPT_ROOT}/lib:$LD_LIBRARY_PATH
fi

if [ -z "$LIBRARY_PATH" ]; then  
   export LIBRARY_PATH=/usr/lib64
else
   export LIBRARY_PATH=${LIBRARY_PATH}:/usr/lib64
fi

#compile_popt 
#update_eilmer3 update
#install_cea
build_eilmer3 
#install_mod 

