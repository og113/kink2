#!/bin/bash

echo "finding modes about soliton"
echo ""
cd ../mpi
matlab -nodisplay -r "M = loadSquareMatrix('../kink2/data/staticShootingDDS.dat');M = -M;[V,D] = eig(M);V0 = V(:,1);printVector(V0,'../kink2/data/stable/staticModes.dat');quit"
echo "soltion modes printed to data/stable/staticModes.dat"
echo ""
cd ../kink2
