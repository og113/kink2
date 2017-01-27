#gnuplot program to plot energy versus number of particles, or something similar

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
set title "N vs theta"
set xlabel "theta"
set ylabel "N"
plot "results/mainResults_constW_1.dat" using 7:9 with points

pause -1
