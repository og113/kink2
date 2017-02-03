#!/bin/bash
# update the sqlite database with new results
## to do: should get it to take databaseFile and table as inputs (with the given defaults if not input)

if [ -z "$1" ]
then
	echo "must supply input file to updateResults";
	exit 1;
elif [ "$#" -eq 1 ]
then
	echo "must supply table in database";
	exit 1;
	databaseFile="db/nr/kink2.sqlite";
elif [ "$#" -eq 2 ]
then
	resultsFile=$1;	
	table=$2;
	databaseFile="db/nr/kink2.sqlite";
elif [ "$#" -eq 3 ]
then
	resultsFile=$1;	
	table=$2;
	databaseFile=$3;
else
	echo "too many arguments for updateResults";
	exit 1;
fi

# filenames generated from input
tempResultsFile=$(echo "$resultsFile" | sed 's/results\/nr/temp/');
databaseResultsFile=$(echo "$resultsFile" | sed 's/results/db/');
deltaResultsFile=$(echo "$databaseResultsFile" | sed 's/\([^\/]\+\)$/delta_\1/');

echo "input      : $resultsFile";
echo "updating   : $databaseResultsFile";
echo "           : $deltaResultsFile";
echo "           : $databaseFile - $table";

# copying results into temp file
cat $resultsFile > $tempResultsFile

# getting end part of filename
resultsFileID=${resultsFile##*/}

# adding $resultsFileID to end of temp file
sed -i "s/$/,$resultsFileID/" $tempResultsFile

# putting difference between temp and database files in delta
touch $databaseResultsFile
comm -23 <(sort $tempResultsFile) <(sort $databaseResultsFile) > $deltaResultsFile

deltaLines=$(wc -l < $deltaResultsFile)
echo "new results: $deltaLines";

if [ "$deltaLines" -gt 0 ]
then

# some "here documents" follow. note that there must be no trailing space after the last "HERE"
# this imports results to the sqlite database
sqlite3 $databaseFile <<HERE
.separator ","
.import "$deltaResultsFile" "$table"
HERE

if [ "$?" -eq 0 ]
then
	cat $deltaResultsFile >> $databaseResultsFile
fi

fi

# views to export
views=""
if [ "$table" == "nr1" ];
then
	views="W_Tb_Theta_E_Num_LoR_5_pot_3"
fi
	
echo "updating views:"

# exporting views
for v in $views;
do
echo $v
sqlite3 $databaseFile <<HERE
.mode csv
.separator " "
.header off
.output "db/nr/plotsHand/${v}.dat"
SELECT * FROM "${v}";
HERE
done

