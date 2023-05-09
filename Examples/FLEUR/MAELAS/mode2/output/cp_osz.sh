#!/bin/bash

SCRDIR=/scratch/ # main folder where we performed FLEUR calculations. It should be consistent with the jsub file

offset=69000.0 # offset in total energy to avoid a very large number and allow MAELAS to read it correctly

smax=0.018  #maximum applied strain. WARNING! It should be equal to the value set in the generation step in file run.sh

nmax=7 #number of deformed structures. WARNING! It should be equal to the value set in the generation step in file run.sh


SG=`grep "Space group number =" output_generation.dat | awk '{printf"%i",$5}'`


if [[ $SG -le 230 ]] && [[ $SG -ge 207 ]]; then
	nb=2 #number of anisotropic magnetoleastic constants
elif [[ $SG -le 194 ]] && [[ $SG -ge 177 ]]; then
	nb=4
elif [[ $SG -le 167 ]] && [[ $SG -ge 149 ]]; then
	nb=6
elif [[ $SG -le 142 ]] && [[ $SG -ge 89 ]]; then
	nb=5
elif [[ $SG -le 74 ]] && [[ $SG -ge 16 ]]; then
	nb=9
else
	echo "The space group number ${SG} is not supported by current version of MAELAS"
	exit 1
fi

for (( j=1 ; j<=${nb} ; j++ ))
do

for (( i=1 ; i<=${nmax} ; i++ ))
do 

for k in 1 2
do 

osz=OSZICAR_${j}_${i}_${k}

cp OSZICAR ./${osz}

grep -H "total energy" ${SCRDIR}/P_${j}_${i}_${k}/out | tail -1 > ene.dat  # extraction of total energy from FLEUR output file

E0=`grep "total energy=" ene.dat | awk '{printf"%.10f\n",$5}'`

a=27.2114079527 # conversion factor between htr and eV
E0ev=`echo "$E0*$a+$offset" | bc -l | awk '{printf"%.10f\n",$1}'`
sed -i "s/energy/${E0ev}/"  ./${osz}
rm ene.dat


done

done

done

ml Python/3.8.2-GCCcore-9.3.0

maelas -mode 2 -d -s ${smax} -n ${nmax} -b > output_derivation.dat  # derivation of anisotropic magnetoelastic constants and anisotropic magnetostrictive coefficients with -mode 2 of MAELAS, the results are printed in file output_derivation.dat

