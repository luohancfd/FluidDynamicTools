#!/bin/bash
#if [ -z "$STY" ]; then exec screen -dm -S InstallUS3D /bin/bash "$0"; fi
GLIBCV=$(ldd --version)
module load intel
MODULE_DIR=$HOME/.local/modulefiles
LICENSE_FILE=XXX
LSERVER1=XXXX
HOST1=XXX
LSERVER2=XXX
SRC_ID=XXX
ROOT_DIR=$HOME/Soft/CFD/US3D
US3D_src=$ROOT_DIR/Src/$SRC_ID
INSTALL_DIR=$ROOT_DIR/Install
ParMETIS_DIR=${ROOT_DIR}/Addons/ParMETIS


SET_UP_SERVER=0
DOWNLOAD=0

### download US3D code to $US3D_src, put parmetis to $ParMETIS_src
if [ $DOWNLOAD = '1' ]
then
   mkdir -p $ROOT_DIR/Src
   rsync -av --progress cts:~/Soft/CFD/US3D/Src/$SRC_ID $ROOT_DIR/Src/
fi

US3D_ver=$(echo ${US3D_src} | sed -E 's/.*RC([0-9]+\.[0-9]+).*/\1/')
mkdir -p $INSTALL_DIR

### get ParMETIS
ParMETIS_src=$(find ${ROOT_DIR}'/Addons' -name "parmetis*.tar.gz")
if [ -z "$ParMETIS_src"];
then 
  mkdir -p $ParMETIS_DIR
  cd $ROOT_DIR/Addons
  wget -N http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
  ParMETIS_src=$(find ${ROOT_DIR}'/Addons' -name "parmetis*.tar.gz")
fi
   
ParMETIS_tmp=${ROOT_DIR}'/Addons/tmp'
mkdir -p $ParMETIS_tmp

### extact and compile ParMETIS
if [ -f ${ParMETIS_DIR}"/lib/libmetis.a" ] && [ -f ${ParMETIS_DIR}"/lib/libparmetis.a" ]
then
   echo "ParMETIS directory: "${ParMETIS_DIR}
else
   tar -zxf ${ParMETIS_src} --strip-components 1 -C ${ParMETIS_tmp}
   sed -i 's/#define IDXTYPEWIDTH [0-9]\+/#define IDXTYPEWIDTH 32/g' ${ParMETIS_tmp}/metis/include/metis.h
   cd ${ParMETIS_tmp}
   make config prefix=${ParMETIS_DIR}
   make
   make install
   cp ${ParMETIS_tmp}/build/Linux-x86_64/libmetis/libmetis.a ${ParMETIS_DIR}/lib/
   rm -rf ${ParMETIS_tmp}
   cd ${ROOT_DIR}
fi
export METIS_HOME=${ParMETIS_DIR}

### compile US3D
mkdir -p ${INSTALL_DIR}
echo "Start compiling US3D"
cd $US3D_src
./install.pl --cfam='intel' --prefix=${INSTALL_DIR} --debug

### set up license server
if [ $SET_UP_SERVER = '1' ]
then
   mkdir -p $ROOT_DIR/License
   cp $US3D_src/files/license/* $ROOT_DIR/License
   cd $ROOT_DIR/License
fi   

### add license switcher
cat << EOF > $INSTALL_DIR/us3d/bin/us3d-switch.exe
#!/bin/sh -l
if [ -z "\$1" ]; then
  lserver="$LSERVER1
elif [ "\$1" == "local" ]; then
  lserver="localhost"
elif [ "\$1" == "us3d" ]; then
  lserver="$LSERVER2
else
  lserver="\$1"
fi
if [ -z "\$1" ]; then
  lport=29750
elif [ -z "\$2" ]; then
  lport=29750
else
  lport=\$2
fi

export US3D_LICENSE=\$lserver:\$lport
export SG_LICENSE_FILE=\$lport@\$lserver
EOF
chmod +x $INSTALL_DIR/us3d/bin/us3d-switch.exe

### add to module
cd $MODULE_DIR
mkdir -p us3d
cd us3d
cat << EOF > ${US3D_ver}.lua
whatis("Version: ${US3D_ver}")
whatis("Description: Sets the environment for running US3D")
always_load("intel")
local  hostname  = capture("hostname")
local  us3dhome  =   "${INSTALL_DIR}/us3d"
local  us3dbin   =   "${INSTALL_DIR}/us3d/bin"

setenv('METIS_HOME',  '${ParMETIS_DIR}')
setenv('US3D_HOME'  ,  us3dhome         )
set_alias('us3d-switch', 'source ${INSTALL_DIR}/us3d/bin/us3d-switch.exe')
prepend_path ('PATH', us3dbin) 

-- If using a license file, set US3D_LICENSE to the path to the
-- license file.  If using a license server, then US3D_LICENSE
-- can either be the path to the license file, or the host and
-- port of the license server.  In this case replace the string
-- <hostname> with the actual hostname of the license server. 
-- setenv US3D_LICENSE <hostname>:29750
-- setenv US3D_LICENSE us3dhome/license/us3d.lic

-- You can add any required modules here
-- module load (your compiler) -or-
-- module prereq  (your compiler)
if string.find(hostname,"$HOST1") then
  setenv("US3D_LICENSE","$LSERVER1:29750")
  setenv("SG_LICENSE_FILE","29750@$LSERVER1")
else
  setenv("US3D_LICENSE","$LSERVER2:29750")
  setenv("SG_LICENSE_FILE","29750@$LSERVER2")

EOF


