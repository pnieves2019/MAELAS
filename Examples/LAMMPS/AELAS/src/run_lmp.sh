#!/bin/bash

inputlmp=in.Fe    # input file name for LAMMPS

exelmp=lmp_serial  # name of the executable file for LAMMPS


for i in {1..3}

do

for j in {1..13}
 
do
	

cp ${inputlmp} in.lmp_0


sed -i "s/AA/$i/"  ./in.lmp_0
sed -i "s/BB/$j/"  ./in.lmp_0


${exelmp} < in.lmp_0 > results_${i}_${j}


ene=`grep "energy" energy_${i}_${j}.dat | awk '{print $2}'`


if [ $j -gt 9 ];
then
cp OSZICAR OSZICAR_0${i}_0${j}

sed -i "s/energy_lmp/$ene/"  ./OSZICAR_0${i}_0${j}

else
	
cp OSZICAR OSZICAR_0${i}_00${j}

sed -i "s/energy_lmp/$ene/"  ./OSZICAR_0${i}_00${j}
fi


done

done


AELAS -d 
