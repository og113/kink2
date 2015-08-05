#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
if (outFile ne 'gui') set term png size 1600,800; set output outFile;


unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "E vs N"
set xlabel "N"
set ylabel "E"
plot "results/mainResults2.dat" using ($11/10.5):($10/18.9) with points

pause -1
