#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
#set terminal postscript eps color enhanced size 10,5
#set output 'figure.eps';

unset log
unset label
unset key
set autoscale
set xtic auto
set ytic auto
set title "E vs T"
set xlabel "T"
set ylabel "E"
plot "results/main_pot_3.dat" using 5:($8/18.9) with points

pause -1
