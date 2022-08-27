#!/bin/bash

atomsk INPOS.lmp pos


maelas -mode 3 -g -n 7 -s 0.01

for i in {1..1} # cubic crystal

do

for j in {1..7}

do

cp POS_${i}_${j}.vasp POSCAR_${i}_${j}

atomsk POSCAR_${i}_${j} -alignx lmp -ow

rm POSCAR_${i}_${j}


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



