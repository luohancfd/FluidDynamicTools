#!/bin/sh -e
# syntax: ./plot_hist.sh FILE I_BLOCK J_BLOCK
F_HEAD=$1
T_START=$2
F_PATTERN=${F_HEAD}.hist

# findout block to process
for i in $(ls); do
#  echo ${i%.b*}
#  echo ${i##*.}
  if [ "${i%.b*}" == ${F_PATTERN} ]; then
    if [[  ${i##*.} == b* ]]; then
      bk=${i##*.b}
      awk -v block=${bk} '!/^#/ {printf "%s\t%d\t%d\t%8.4e\t%8.4e\n",block,$2,$3,$5,$6}' $i >> 1.txt
  #   $5 is pos_x,$6 is pos_y
    fi
  fi
done
awk '!seen[$0]++' 1.txt >hist_point.txt


# create awk file
cat <<EOF> extract-history.awk
BEGIN {print "# t(ms)\t Mach\t p(kPa)\t mass_N2\t mass_N\t T(K)\t Tv(K)";}
\$2 == idx && \$3 == jdx {
t = \$1 ;
rho = \$9;
u = \$10 ;
v = \$11 ;
a = \$14 ;
p = \$13 ;
rho = \$9; 
mass_N2 = \$23;
mass_N = \$24;
T = \$27;
TV = \$29;
print t *1000.0 , sqrt ( u * u + v * v )/ a , p /1000.0,mass_N2*100.0,mass_N*100.0,T,TV, rho;
}
EOF

rm -f *.pdf *.ps

while read p; do
  IFS=$'\t' read -r -a array  <<< "$p"
  echo "blk=${array[0]} i=${array[1]} j=${array[2]}"
  DAT_FILE=b${array[0]}.i${array[1]}.j${array[2]}.dat
  DAT_HEAD=${DAT_FILE%.dat}
  TITLE="B=${array[0]} x=${array[3]} y=${array[4]}"
  awk -v idx="${array[1]}" -v jdx="${array[2]}" -f extract-history.awk \
   < ${F_PATTERN}".b"${array[0]} \
   > ${DAT_FILE} 
  
gnuplot << EOF
set terminal postscript enhanced color font "Helvetica" 20 
set output '${DAT_HEAD}.ps'
set style line 1 lw 5
set title 'Density history ${TITLE}'
set xlabel 'time (ms)'
set ylabel 'Rho (kg/m^3)'
set xrange [${T_START}:]
plot '${DAT_FILE}' using 1:8 lw 4 title "Rho (kg/m^3)" with lines
set title 'Pressure history ${TITLE}'
set xlabel 'time (ms)'
set ylabel 'P (kPa)'
set xrange [${T_START}:]
plot '${DAT_FILE}' using 1:3 lw 4 title "P (kPa)" with lines
set title 'Temperature history ${TITLE}'
set xlabel 'time (ms)'
set ylabel 'Temperature (K)'
set xrange [${T_START}:]
plot '${DAT_FILE}' using 1:6 lw 4 title "T" with lines,\
     '${DAT_FILE}' using 1:7 lw 4 title "Tv" with lines
set title 'Mass Fraction ${TITLE}'
set xlabel 'time (ms)'
set ylabel 'MF (%)'
set xrange [${T_START}:]
#set yrange [-10:110]
plot '${DAT_FILE}' using 1:4 lw 4 title "N_2" with lines
#,\
#     '${DAT_FILE}' using 1:5 lw 4 title "N" with lines
EOF
  ps2pdf ${DAT_HEAD}.ps ${DAT_HEAD}.pdf
done < hist_point.txt

gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=$1.pdf b*.pdf

while read p; do
  IFS=$'\t' read -r -a array  <<< "$p"
#  echo "blk=${array[0]} i=${array[1]} j=${array[2]}"
  DAT_FILE=b${array[0]}.i${array[1]}.j${array[2]}.dat
  DAT_HEAD=${DAT_FILE%.dat}
  TITLE="B=${arrayp[0]} I=${array[1]} J=${array[2]}"
  rm ${DAT_HEAD}.pdf
done < hist_point.txt
rm 1.txt
nohup evince $1.pdf >/dev/null 2>&1 &
