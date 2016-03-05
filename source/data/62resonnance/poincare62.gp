set terminal postscript eps size 4.5,3.37 enhanced color font 'Helvetica,20' linewidth 2
set output 'poincare62.eps'

set xlabel "x [mm]"
set ylabel "p_x  [mm]"
set xrange [-5:5]
set yrange [-5:5]
#set size 1,1
#set origin 0,0
#set grid
#set key opaque
#set key left top

# This plots the big plot
plot "raw_K0.dat" using (($1)*1000):(($3)*1000) ps 0.1  notitle
