#!/bin/bash

sflag=false

# checking if outFile required and getting filename if so
options=':s'
OPTIND=1
while getopts $options option
do
	case $option in
		s  ) sflag=true;;
		\? ) echo "Unknown option argument -$OPTARG" >&2; exit 1;;
		:  )
	esac
done
shift $((OPTIND-1))

if [ -z "$1" ]
then
	echo "must supply input file"
else
	echo "changing mathematica output to c++ input in $1"
	# generic changes
	#sed -i 's/\([1-9]\+\)\*/\1.0\*/g' $1
	#sed -i 's/\([1-9]\+\)\.\*/\1.0\*/g' $1
	sed -i 's/\*\([1-9]\+\)/\*\1.0/g' $1
	sed -i 's/\*\([1-9]\+\)\./\*\1.0/g' $1
	sed -i 's/\/\([1-9]\+\)/\/\1.0/g' $1
	sed -i 's/\/\([1-9]\+\)\./\/\1.0/g' $1
	sed -i 's/\([1-9]\+\)\//\1.0\//g' $1
	sed -i 's/\([1-9]\+\)\.\//\1.0\//g' $1
	sed -i 's/Pi/PI/g' $1
	sed -i 's/Im(/imag(/g' $1
	sed -i 's/Re(/real(/g' $1
	sed -i 's/Power(/pow(/g' $1
	sed -i 's/pow(E,/exp(/g' $1
	sed -i 's/Sqrt(/sqrt(/g' $1
	sed -i 's/Cos(/cos(/g' $1
	sed -i 's/Sin(/sin(/g' $1
	sed -i 's/Tan(/tan(/g' $1
	sed -i 's/Cosh(/cosh(/g' $1
	sed -i 's/Sinh(/sinh(/g' $1
	sed -i 's/Tanh(/tanh(/g' $1
	sed -i 's/\\\[Mu\]/mu/g' $1
	sed -i 's/\\\[Nu\]/nu/g' $1
	sed -i 's/\\\[Rho\]/rho/g' $1
	sed -i 's/\\\[Sigma\]/sigma/g' $1
	sed -i 's/\\\[Delta\]/delta/g' $1
	sed -i 's/\\\[Kappa\]/kappa/g' $1
	sed -i 's/\\\[Tau\]/tau/g' $1
	sed -i 's/\\\[Beta\]/beta/g' $1
	sed -i ':a;N;$!ba;s/\(\s*\)\([+-]\)\(\s*\)\n/ \\\n\2 /g' $1
	sed -i 's/^\([+-]\)\(\s*\)/\1 /g' $1
	sed -i 's/^\(\s*\)/ /g' $1
	sed -i 's/\/$/\/\\/g' $1
	sed -i 's/\*$/\*\\/g' $1
	# specific changes, for worldline n-r calculations
	if $sflag
	then echo "making specific changes, for n-r calculations"
			sed -i 's/(-1 + \([c-v]\))/(m\1)/g' $1
			sed -i 's/(1 + \([c-v]\))/(p\1)/g' $1
			sed -i 's/2\*(\([c-v]\+\))/2*\1/g' $1
			sed -i 's/a(\([c-v]\+\))/real(f(\1))/g' $1
			sed -i 's/b(\([c-v]\+\))/imag(f(\1))/g' $1
			sed -i 's/Complex(0,1)/ii/g' $1
	fi
fi
