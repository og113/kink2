#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
if (outFile ne 'gui') set term png size 1600,800; set output outFile;

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "E vs T"
set xlabel "T"
set ylabel "E"
plot "results/mainResults.dat" using 5:($8/18.9) with points

pause -1
