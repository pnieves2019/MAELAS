#!/bin/bash

atomsk cell.lmp pos

mv POSCAR POSCAR_rlx

maelas -g -mode 1 -i POSCAR_rlx -s 0.005 -n 9

for i in {1..2}

do

for j in {1..9}

do

atomsk POSCAR_${i}_${j} -alignx lmp -ow


myFile="POSCAR_${i}_${j}.lmp"


n=1
while read line; do
	if [ $n -gt 15 ];
	then

	echo  $line "  0.0 0.0 0.0 0.0" >> POSCAR_${i}_${j}_new.lmp
	else
	echo  $line  >> POSCAR_${i}_${j}_new.lmp
	fi

	n=$((n+1))
done < $myFile

mv POSCAR_${i}_${j}_new.lmp POSCAR_${i}_${j}.lmp


done

done
