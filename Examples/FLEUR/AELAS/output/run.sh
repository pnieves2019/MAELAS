#!/bin/bash

ml purge
ml intel

nmax=7 #number of deformed structures (it depends on your AELAS compilation default = 13)

mv POSCAR INPOS

AELAS -g   # generate deformed unit cells with AELAS

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
	cp POS_0${j}_0${i}.vasp POSCAR
	posc=POSCAR
else
	cp POS_0${j}_00${i}.vasp POSCAR
	posc=POSCAR
fi


mkdir P_${j}_${i}

cp ./${posc} ./P_${j}_${i}/
cp ./jsub ./P_${j}_${i}/
cd P_${j}_${i}

sed -i "s/AA/$j/"  ./jsub
sed -i "s/BB/$i/"  ./jsub


ml purge
ml intel
ml libxml2
ml ifort
ml Python/3.8.2-GCCcore-9.3.0

vasp2cif POSCAR   # conversion from VASP to cif

cif2cell -f POSCAR.cif -p fleur -o inp_file  # conversion from cif to FLEUR

sed -i "/oldfleur=f/ c\&input cartesian=f oldfleur=f film=f /" inp_file

echo -e "\n&comp gmax=12.0 kmax=5.0 /" > gmax0.dat  # we change plane-wave cutoff

cat inp_file gmax0.dat >> inp_file0

mv inp_file0 inp_file

inpgen -f inp_file -kpt testset#gamma@grid=20,20,20   # we define new k-point mesh


sed -i 's/itmax="15"/itmax="60"/' inp.xml   # we increase maximum number of iteration in self-consistent field (scf) cycle

sed -i 's/kPointListSelection listName="default"/kPointListSelection listName="testset"/' inp.xml   # we set our new k-point mesh

qsub jsub   # run FLEUR 


cd ..

rm $posc

done

done



