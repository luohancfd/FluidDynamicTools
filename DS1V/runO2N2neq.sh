#!/bin/bash
PWD=`pwd`
ROOTFILE=$PWD'/runO2N2neq.sub'


folder='/scratch/rice/l/luo160/0DRate/Diatom/NEq/O2N2Double'
cd $folder
TEMPS="8000 10000 15000 20000"
VTEMP="1000 3000 4000 5000 6500 8000 10000 12500 15000 17500 20000"
EXEfile=$(find . -maxdepth 2 -name "*.exe")
EXEbase=${EXEfile:2}
EXEbase=${EXEbase%.exe}
EXEfile=$folder'/'${EXEfile:2}
cat << EOF > $ROOTFILE
#PBS -l nodes=1:ppn=20
#PBS -l walltime=4:00:00
#PBS -q standby
#PBS -V
#PBS -N O2N2Neq

module load gcc
EOF

for TEMP in $TEMPS; do
for VT in $VTEMP; do
wd=$folder"/T="$TEMP"/Tv="$VT
mkdir -p $wd
CEXE=$wd/$EXEbase"_"$TEMP"_"$VT".exe"
CEXEbase=$EXEbase"_"$TEMP"_"$VT".exe"
cp $EXEfile $CEXE
cp $folder/*.mod $folder/DS1VD.in $wd/

cat << EOF > ${wd}/input.txt
3       !new run
0       !seed
0       !continue with current cells
${TEMP} !translational temperature
${VT} !vibrational temperature
EOF
cat << EOF >> $ROOTFILE
cd $wd
nohup ./$CEXEbase < input.txt > ds1v.log 2>&1  &
EOF
done
done

