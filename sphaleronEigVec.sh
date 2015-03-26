#!/bin/bash

if [ -z "$1" ]
  then
    echo "Must supply argument for L"
else
	echo "finding negative eigenvector of sphaleron with L = $1"
	echo ""
	./sphaleron -r1 $1 #1st argument is L
	cp data/stable/D1_L_$1.dat ../mpi/data/D1.dat
	cp data/stable/D2_L_$1.dat ../mpi/data/D2.dat
	cd ../mpi
	matlab -nodisplay -r "[D1,D2,diff,maximum] = compareDDS;[V,D] = eigs(D2,2,-20);V0 = V(:,1);printVector(V0,'../kink2/data/stable/sphaleronEigVec_L_5.dat');quit"
	echo "negative eigenvector printed to data/stable/eigVec_pot_3_L_5.dat"
	echo ""
	cd ../kink2
fi
