#gnuplot program to plot W against E from mainAction.dat

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "W(E,N)"
set xlabel "E"
set ylabel "W"
set grid
plot "results/mainResults_constW_1.dat" using 10:13 with points

pause -1
