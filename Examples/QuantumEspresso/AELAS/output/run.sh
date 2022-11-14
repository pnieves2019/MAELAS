#!/bin/bash

nmax=13 #number of deformed structures (it depends on your AELAS compilation, default nmax=13)

atomsk structure.pw -frac pos   #conversion between QuantumEspresso and VASP format

sed -i "/Direct/ c\Cartesian" POSCAR

mv POSCAR INPOS

AELAS -g

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
	osz=OSZICAR_0${j}_0${i}
else
	cp POS_0${j}_00${i}.vasp POSCAR
	posc=POSCAR
	osz=OSZICAR_0${j}_00${i}
fi


atomsk ${posc} pw  #conversion between QuantumEspresso and VASP format
sed -i "/ATOMIC_POSITIONS/ c\ATOMIC_POSITIONS angstrom" ${posc}.pw
sed -i "/2 2 2  0 0 0/ c\10 10 10  0 0 0" ${posc}.pw    # k-points
sed -i "/pseudo_dir/ c\ pseudo_dir = '/home/$USER/qe/pseudo/'" ${posc}.pw   # pseudopotential directory
sed -i "/conv_thr/ c\conv_thr =  1.0d-6\n electron_maxstep=300" ${posc}.pw   # convergence tolerance in scf and max number scf steps 
sed -i "/Fe  55.845/ c\Fe  55.845  Fe.pbesol-spn-kjpaw_psl.1.0.0.UPF" ${posc}.pw #pseudopotential
sed -i "/&SYSTEM/ c\&SYSTEM \n nspin = 2\n starting_magnetization(1)=0.7\n occupations='smearing'\n smearing='gauss'\n degauss=0.02" ${posc}.pw  # spin-polirized calculation, starting magnetization, smearing
sed -i "/ecutwfc/ c\ ecutwfc = 60.0\n ecutrho = 480.0" ${posc}.pw   # cutoff



mkdir P_${j}_${i}

cp ${posc}.pw ./P_${j}_${i}/
cp OSZICAR ./P_${j}_${i}/${osz}
cd P_${j}_${i}

mpirun pw.x < ${posc}.pw > scf_${j}_${i}.out

E0=`grep "!    total energy              =" scf_${j}_${i}.out | awk '{printf"%.10f\n",$5}'`
a=13.605698066 # conversion factor between Ry and eV
E0ev=`echo "$E0*$a" | bc -l | awk '{printf"%.6f\n",$1}'`
sed -i "s/energy/${E0ev}/"  ./${osz}
cp ${osz} ../

cd ..

rm $posc

done

done


AELAS -d

