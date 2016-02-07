#!/bin/bash
#usage: ./a.sh label
if [ -z "$1" ]; then
    echo "No argument supplied"
	exit 1
else
    label=$1
fi
outfile1=bbeat_$label.eps
outfile2=dquad_$label.eps
outfile3=harms_$label.eps
awk '{ print $3," ",$8," ",$11 }' FODO0_$label.twiss > b0
awk '{ print $8," ",$11 }' FODO1_$label.twiss > b1
awk '{ print $8," ",$11 }' FODO2_$label.twiss > b2

paste b0 b1 b2 | uniq > ball
rm b0 b1 b2

cat <<HERE1 | gnuplot
set term post eps enhanced color font 'Helvetica,28' \
dashed dashlength 2.0 linewidth 1.0 butt
set key center bottom samplen 2 spacing .85 font ',28' maxrows 2
set output '$outfile1'

#set label "N=$N" at screen 0.3,0.9
set xlabel "s (m)"
set ylabel "Beta beating (%)" offset 1,0
set xrange [0:324]
set yrange [-1.6:1.6]
set ytics 0.5
set xzeroaxis
p 'ball' u (\$1):(100*(\$4-\$2)/\$2) w l lw 3 lt 1 lc 1 t "x-before",\
      '' u (\$1):(100*(\$5-\$3)/\$3) w l lw 3 lt 1 lc 3 t "z-before",\
      '' u (\$1):(100*(\$6-\$2)/\$2) w l lw 3 lt 3 lc 1 t "x-after",\
      '' u (\$1):(100*(\$7-\$3)/\$3) w l lw 3 lt 3 lc 3 t "z-after"
HERE1

cp dqba_$label.dat dquad.dat
cat <<HERE2 | gnuplot
set term post eps enhanced color font 'Helvetica,24'
set output '$outfile2'
set style fill solid 0.5 noborder
set xlabel 'Quadrupole index'
set ylabel '{/Symbol D}K_1L (10^{-3} m^{-1})' offset 1,0
set xrange [0:37]
set yrange [-1.0:1.0]
set xzeroaxis
set xtics 5 font ",20"

#set label "15 trim quads used" at screen 0.3,0.9
set boxwidth 0.3
set grid
p '<cat -n dquad.dat' u (\$1-0.15):(\$2*1000) w boxes lc 1 t 'before',\
                  '' u (\$1+0.15):(\$3*1000) w boxes lc 2 t 'after'
HERE2

cat harm_$label.dat > harm.dat
echo >> harm.dat
echo >> harm.dat
cat dnxz_$label.dat >> harm.dat
cat <<HERE3 | gnuplot
set term post eps enhanced color font 'Helvetica,28'
set output '$outfile3'
set multiplot
set origin 0,0
set size 1,1
set xlabel 'Harmonic index'
set ylabel 'q_{nx} / q_{nz}' offset 2,0
set xrange [0:12]
set yrange [0:0.012]
set xtics 1 nomirror out 1,1,10 font ",28" 
set ytics 0,0.002,0.01

set boxwidth 0.2
#set grid
set key font ",28" Left reverse samplen 2 at 4.8,0.0115
p '<cat -n harm.dat' u (\$1-0.30):2 w boxes lc 1 lt 1 t 'x-before' fs solid noborder,\
                  '' u (\$1+0.10):3 w boxes lc 3 lt 1 t 'z-before' fs solid noborder,\
                  '' u (\$1-0.10):4 w boxes lc 1 lt 1 t 'x-after'  fs pattern 4 border,\
                  '' u (\$1+0.30):5 w boxes lc 3 lt 1 t 'z-after'  fs pattern 5 border 

set origin 0.53, 0.5
set size 0.45,0.45
unset grid
unset ylabel
set title '{/Symbol D}{/Symbol n}_x / {/Symbol D}{/Symbol n}_z' font ",28" offset 0,-0.9
set ytics font ",26" offset 0.5
unset xtics
unset xlabel
unset key
set xzeroaxis
set xrange [12:14]
set yrange [-0.002:0.002]
set ytics 0.001
replot
unset multiplot
HERE3

rm ball harm.dat dquad.dat
echo "$outfile1, $outfile2, $outfile3 plotted!"
