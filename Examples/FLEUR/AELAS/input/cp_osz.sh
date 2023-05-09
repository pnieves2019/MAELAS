#!/bin/bash

ml purge
ml intel

SCRDIR=/scratch/   # main folder where we performed FLEUR calculations. It should be consistent with the jsub file

offset=69000.0 # offset in total energy to avoid a very large number and allow AELAS to read it correctly

nmax=7 #number of deformed structures (it depends on your AELAS compilation default = 13)

SG=`head -1 DEFMOD | awk '{printf"%i",$1}'`

if [[ $SG -le 230 ]] && [[ $SG -ge 195 ]]; then
	nc=3 #number of independent elastic constants
elif [[ $SG -le 194 ]] && [[ $SG -ge 168 ]]; then
	nc=5
elif [[ $SG -le 167 ]] && [[ $SG -ge 149 ]]; then
	nb=6
elif [[ $SG -le 148 ]] && [[ $SG -ge 143 ]]; then
	nb=7
elif [[ $SG -le 142 ]] && [[ $SG -ge 89 ]]; then
	nb=6
elif [[ $SG -le 88 ]] && [[ $SG -ge 75 ]]; then
	nb=7
elif [[ $SG -le 74 ]] && [[ $SG -ge 16 ]]; then
	nb=9
elif [[ $SG -le 15 ]] && [[ $SG -ge 3 ]]; then
	nb=13
elif [[ $SG -le 2 ]] && [[ $SG -ge 1 ]]; then
	nb=21
else
	echo "The space group number ${SG} is not in the range 230-1"
	exit 1
fi

for (( j=1 ; j<=${nc} ; j++ ))
do

for (( i=1 ; i<=${nmax} ; i++ ))
do 

if [[ $i -ge 10 ]] ; then
	osz=OSZICAR_0${j}_0${i}
else
	osz=OSZICAR_0${j}_00${i}
fi

cp OSZICAR ${osz}

grep -H "total energy" ${SCRDIR}/P_${j}_${i}/out | tail -1 > ene.dat   # extraction of total energy from FLEUR output file

E0=`grep "total energy=" ene.dat | awk '{printf"%.10f\n",$5}'`

a=27.2114079527 # conversion factor between htr and eV
E0ev=`echo "$E0*$a+$offset" | bc -l | awk '{printf"%.10f\n",$1}'`
sed -i "s/energy/${E0ev}/"  ./${osz}
rm ene.dat



done

done

ml purge
ml intel


AELAS -d   # derivation of elastic constants with AELAS, the calculated elastic tensor (Cij) is written in file ELADAT

