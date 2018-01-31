#!/bin/bash
ROOT_DIR=${PWD}
POPT_SRC=${ROOT_DIR}/Addons/popt
CGNS_SRC=${ROOT_DIR}/Addons/cgns
POPT_ROOT=${ROOT_DIR}/lib
CFCFD3_SRC=${ROOT_DIR}/cfcfd3
EILMER3_BUILD=${CFCFD3_SRC}/app/eilmer3/build
EILMER3_DIR=${ROOT_DIR}/e3bin
POSHAX_BUILD=${CFCFD3_SRC}/app/poshax3/build
MODULE_FILE=${HOME}/.local/modules
DOWNLOAD=1

#source /opt/intel/16.0.3/bin/compilervars.sh intel64
module load intel
cd ${ROOT_DIR}
if [ $DOWNLOAD = 1 ]
then
   hg clone https://source.eait.uq.edu.au/hg/cfcfd3 cfcfd3
else
   scp rice:~/Soft/CFCFD/cfcfd3.tar.gz ${ROOT_DIR}/
   mkdir -p ${CFCFD3_SRC}
   tar -zxf cfcfd3.tar.gz --strip-components 1 -C ${CFCFD3_SRC}
fi
#
wget http://rpm5.org/files/popt/popt-1.16.tar.gz
curl -O https://www.grc.nasa.gov/WWW/CEAWeb/CEA+Fortran.tar.Z
#wget wget https://github.com/CGNS/CGNS/archive/v3.1.4.tar.gz -O cgns.tar.gz

##### extracting lib
mkdir -p ${ROOT_DIR}/Addons
mkdir -p ${POPT_SRC}
# mkdir -p ${CGNS_SRC}
tar -zxf popt-1.16.tar.gz --strip-components 1 -C ${POPT_SRC}
# tar -zxf cgns.tar.gz --strip-components 1 -C ${CGNS_SRC}

#### compile popt and cgns
mkdir -p ${POPT_ROOT}
cd ${POPT_SRC}
./configure  --prefix=${POPT_ROOT} --disable-static
make install
# cd ${CGNS_SRC}
# cmake ../cgns

export CPATH=${POPT_ROOT}/include:$CPATH
export LD_LIBRARY_PATH=${POPT_ROOT}/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=${LIBRARY_PATH}:/usr/lib64
##### extracting CFCFD3 code
echo "Extracting src file"
cd ${CFCFD3_SRC}
hg pull -u https://source.eait.uq.edu.au/hg/cfcfd3
## get cea
cd ${ROOT_DIR}
mkdir -p ${CFCFD3_SRC}/extern/cea2
tar -zxf CEA+Fortran.tar.Z --strip-components 1 -C ${CFCFD3_SRC}/extern/cea2
#### compile eilmer3
mkdir -p ${EILMER3_DIR}
cd ${EILMER3_BUILD}
## change some directory
sed -i "s|LLIB +=.*|LLIB += -lpopt -ldl -L\/opt\/local\/lib -L${POPT_ROOT}\/lib|g" "makefile"
sed -i "s/gfortran -m32 -std=legacy\(.*cea2.f\)/gfortran -std=legacy\1/g" ${CFCFD3_SRC}/lib/radiation/build/makefile
sed -i "s/gfortran -m32 -std=legacy\(.*cea2.f\)/gfortran -std=legacy\1/g" ${CFCFD3_SRC}/lib/gas/build/makefile
make TARGET=for_intel_mpi INSTALL_DIR=${EILMER3_DIR} install

#### compile poshax 3
cd ${POSHAX_BUILD}
make TARGET=for_intel INSTALL_DIR=${EILMER3_DIR} install


### add to module
cd ${MODULE_FILE}
mkdir -p cfcfd3
cd cfcfd3
cat << EOF > 3.lua
whatis("Version: 3")
whatis("Description: Sets the environment for running CFCFD3")
load("intel")
prepend_path ('PATH', "${EILMER3_DIR}") 
append_path  ('LD_LIBRARY_PATH', "${EILMER3_DIR}") 
prepend_path ('PYTHONPATH',"${EILMER3_DIR}") 
prepend_path("LUA_PATH","${EILMER3_DIR}/?.lua",";")
prepend_path("LUA_CPATH","${EILMER3_DIR}/?.so",";")
setenv("E3BIN","${EILMER3_DIR}")
EOF


