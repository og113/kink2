# gnuplot program to test various bits of code

a = 10
if (exists("a")) print "a is defined"
if (!exists("b")) print "b is not defined"
