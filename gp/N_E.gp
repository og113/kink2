#gnuplot program to plot energy versus number of particles, or something similar

#if you want to save directly to a file, use the following two lines of code
if (exists("outFile")) \
if (outFile ne 'gui') \
set term png size 1600,800; \
set output outFile; \


unset log
unset label
unset key
#set key bottom
set autoscale
#set xrange [0.5:1.2]
#set yrange [0.5:1.2]
set xtic auto
set ytic auto
set title "N vs E"
set xlabel "E"
set ylabel "N"
#f(x) = m*x+c	           			# define the function to be fit
#m = 1.0; c = 0.01;            		# initial guess for m and c
#fit f(x) "results/mainResults_constT_1.dat" using ($10/18.9):($11/10.5) via m, c
#fitTitle = sprintf("%.2fx+%.2f", m, c)
plot \
	"results/mainResults_constT_1.dat" using ($10/18.9):($11/10.5) with points#, \
	#f(x) title fitTitle

pause -1
