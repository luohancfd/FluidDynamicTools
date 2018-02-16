#!/bin/bash

# clean Nan in Vrat

clean (){
echo "Processing "$2
cd $1
PWD=$(pwd)
files=$(find . -name 'IMF_Vrate.dat')
for i in $files; do
  foldername=$(dirname $i)
  foldername=${foldername:2}
  echo $foldername
  sed -i "s|^\s*#|#|" $i
  sed -i "s|^\s*VARIABLES|VARIABLES|g" $i
  sed -i "s|^\s*ZONE|ZONE|g" $i
  sed -i "s|\s*NaN|    0.000000E+00|g" $i

  if [ "$2" == "O2+N2" ]; then
    sed -i "s|reac           1 |O2+N2$foldername|" $i
    sed -i "s|reac           2 |N2+O2$foldername|" $i
  fi

  if [ "$2" == "NO+N2" ]; then
    sed -i "s|reac           1 |NO+N2$foldername|" $i
  fi
done

files=$(find . -name 'IMF_EV.DAT')
for i in $files; do
  foldername=$(dirname $i)
  foldername=${foldername:2}
  echo $foldername
  sed -i "s|^\s*#|#|" $i
  sed -i "s|^\s*VARIABLES|VARIABLES|g" $i
  sed -i "s|^\s*ZONE|ZONE|g" $i
  sed -i "s|PASSIVEVARLISTWRITE|PASSIVEVARLIST|" $i
done
}

clean "/mnt/e/Research/NonequilibriumGas/Macheret-Fridmann/DSMC_MF/Diatom/Eq/O2N2Double" "O2+N2"
clean "/mnt/e/Research/NonequilibriumGas/Macheret-Fridmann/DSMC_MF/Diatom/Eq/NON2" "NO+N2"

clean "/mnt/e/Research/NonequilibriumGas/Macheret-Fridmann/DSMC_MF/Diatom/NEq/O2N2Double" "O2+N2"
clean "/mnt/e/Research/NonequilibriumGas/Macheret-Fridmann/DSMC_MF/Diatom/NEq/NON2" "NO+N2"

clean "/mnt/e/Research/NonequilibriumGas/Macheret-Fridmann/DSMC_MF/Eq2/O3_wysong" "O3_wysong"

#cd /mnt/e/Research/NonequilibriumGas/Macheret-Fridmann/Vrate
#matlab.exe -nodisplay -nodesktop -nosplash -r VrateNON2
#wait
#cp DSMCVrateNON2.plt /mnt/e/Research/NonequilibriumGas/Macheret-Fridmann/DSMC_MF/Diatom/Eq/NON2/DSMCVrateNON2.plt
