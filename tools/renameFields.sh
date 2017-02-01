#!/bin/bash

# program to change filenames of field files

if [ -z "$1" ]
  then
    echo "must supply input folder"
else
	test=$(ls $1/*_theta_*.data | wc -l)
	if (("$test" > "0")); then
		echo "step 1, renaming theta -> Theta"
		for f in $1/*_theta_*.data;
			do fo=$(echo "$f"|sed "s/_theta_/_Theta_/")
			echo $f "renamed as" $fo
			mv $f $fo;
		done
	else
		echo "no files to rename theta -> Theta"
	fi
	
	test=$(ls $1/*_reg_*.data | wc -l)
	if (("$test" > "0")); then
		echo "step 1, renaming reg -> Reg"
		for f in $1/*_reg_*.data;
			do fo=$(echo "$f"|sed "s/_reg_/_Reg_/")
			echo $f "renamed as" $fo
			mv $f $fo;
		done
	else
		echo "no files to rename reg -> Reg"
	fi

	test=$(ls $1/*LoR_?_Tb_*.data | wc -l)
	if (("$test" > "0")); then
		echo "step 1, renaming to introduce _DE_0"
		for f in $1/*LoR_?_Tb_*.data;
			do fo=$(echo "$f"|sed "s/_Tb/_DE_0_Tb/")
			echo $f "renamed as" $fo
			mv $f $fo;
		done
	else
		echo "no files to introduce _DE_0"
	fi
fi

