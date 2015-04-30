#!/bin/bash

echo "finding modes about soliton"
echo ""
cp data/staticNewtonDDS.dat ../mpi/data/D1.dat
cd ../mpi
matlab -nodisplay -r "[D1,D2,diff,maximum] = compareDDS;[V,D] = eigs(D1,2,-20);V0 = V(:,1);printVector(V0,'../kink2/data/stable/staticModes.dat');quit"
echo "soltion modes printed to data/stable/staticModes.dat"
echo ""
cd ../kink2
