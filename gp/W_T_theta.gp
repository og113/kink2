#gnuplot program to plot W against T and theta from mainAction.dat

#if you want to save directly to a file, use the following two lines of code
if (outFile ne 'gui') set term png size 1600,800; set output outFile;

unset log
unset label
unset key
set autoscale
#set xrange [0.5:1.2]
#set yrange [0.0:0.1]
#set zrange [-2.0:20.0]
set xtic auto
set ytic auto
set title "W(T,theta)"
set xlabel "T"
set ylabel "theta"
set zlabel "W"
set grid
splot "results/mainResults2.dat" using 7:9:13 with points

pause -1
