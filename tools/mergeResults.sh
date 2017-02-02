#!/bin/bash

results="results/nr/nr1.csv results/nr/nr1error.csv"

echo "removing lines with <<<, >>> and === from:"

for f in $results;
do
	echo $f
	sed -i '/<<</d' $f
	sed -i '/>>>/d' $f
	sed -i '/===/d' $f
done
