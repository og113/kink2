#gnuplot program to plot theta against T from mainAction.dat

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
set title "T vs theta"
set xlabel "T"
set ylabel "theta"
set grid
plot "results/mainResults_constW_2.dat" using 7:9 with points

pause -1
