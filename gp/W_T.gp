#gnuplot program to plot W against T from mainAction.dat

#if you want to save directly to a file, use the following two lines of code
if (outFile ne 'gui') set term png size 1600,800; set output 'outFile';

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "W(T,N)"
set xlabel "T"
set ylabel "W"
set grid
plot "results/main_pot_3.dat" using 5:11 with points

pause -1
