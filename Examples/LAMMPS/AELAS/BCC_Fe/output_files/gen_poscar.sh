#!/bin/bash

atomsk INPOS.lmp pos

mv POSCAR INPOS

AELAS -g

for i in {1..3}

do

for j in {1..9}

do

cp POS_0${i}_00${j}.vasp POSCAR_0${i}_00${j}

atomsk POSCAR_0${i}_00${j} -alignx lmp -ow

rm POSCAR_0${i}_00${j}


myFile="POSCAR_0${i}_00${j}.lmp"


n=1
while read line; do
	if [ $n -gt 15 ];
	then

	echo  $line "  0.0 0.0 0.0 0.0" >> POSCAR_0${i}_00${j}_new.lmp
	else
	echo  $line  >> POSCAR_0${i}_00${j}_new.lmp
	fi
	
	n=$((n+1))
done < $myFile

mv POSCAR_0${i}_00${j}_new.lmp POSCAR_0${i}_00${j}.lmp

done

done




for i in {1..3}

do

for j in {10..13}

do


		cp POS_0${i}_0${j}.vasp POSCAR_0${i}_0${j}

		atomsk POSCAR_0${i}_0${j} -alignx lmp -ow

		rm POSCAR_0${i}_0${j}


		myFile="POSCAR_0${i}_0${j}.lmp"


		n=1
		while read line; do
	if [ $n -gt 15 ];
	then

       echo  $line "  0.0 0.0 0.0 0.0" >> POSCAR_0${i}_0${j}_new.lmp
	        else	
	echo  $line  >> POSCAR_0${i}_0${j}_new.lmp
       fi

    n=$((n+1))
	done < $myFile

		mv POSCAR_0${i}_0${j}_new.lmp POSCAR_0${i}_00${j}.lmp

		done

		done
