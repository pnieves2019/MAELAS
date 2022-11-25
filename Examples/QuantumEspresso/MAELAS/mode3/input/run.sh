#!/bin/bash


smax=0.01  #maximum applied strain

nmax=7 #number of deformed structures

atomsk structure.pw -frac pos   #conversion between QuantumEspresso and VASP format

sed -i "/Direct/ c\Cartesian" POSCAR

maelas -mode 3 -g -s ${smax} -n ${nmax} > output_generation.dat

SG=`grep "Space group number =" output_generation.dat | awk '{printf"%i",$5}'`

if [[ $SG -le 230 ]] && [[ $SG -ge 207 ]]; then
	nb=1 #number of isotropic magnetoleastic constants
elif [[ $SG -le 194 ]] && [[ $SG -ge 177 ]]; then
	nb=2
elif [[ $SG -le 167 ]] && [[ $SG -ge 149 ]]; then
	nb=2
elif [[ $SG -le 142 ]] && [[ $SG -ge 89 ]]; then
	nb=2
elif [[ $SG -le 74 ]] && [[ $SG -ge 16 ]]; then
	nb=3
else
	echo "The space group number ${SG} is not supported by current version of MAELAS"
	exit 1
fi

for (( j=1 ; j<=${nb} ; j++ ))
do

for (( i=1 ; i<=${nmax} ; i++ ))
do 

atomsk POSCAR_${j}_${i} pw  #conversion between QuantumEspresso and VASP format
sed -i "/ATOMIC_POSITIONS/ c\ATOMIC_POSITIONS angstrom" POSCAR_${j}_${i}.pw
sed -i "/2 2 2  0 0 0/ c\10 10 10  0 0 0" POSCAR_${j}_${i}.pw    # k-points
sed -i "/pseudo_dir/ c\ pseudo_dir = '/home/$USER/qe/pseudo/'" POSCAR_${j}_${i}.pw   # pseudopotential directory
sed -i "/conv_thr/ c\conv_thr =  1.0d-6\n electron_maxstep=300" POSCAR_${j}_${i}.pw   # convergence tolerance in scf and max number scf steps 
sed -i "/Fe  55.845/ c\Fe  55.845  Fe.pbesol-spn-kjpaw_psl.1.0.0.UPF" POSCAR_${j}_${i}.pw #pseudopotential
sed -i "/&SYSTEM/ c\&SYSTEM \n nspin = 2\n starting_magnetization(1)=0.7\n occupations='smearing'\n smearing='gauss'\n degauss=0.02" POSCAR_${j}_${i}.pw  # spin-polirized calculation, starting magnetization, smearing
sed -i "/ecutwfc/ c\ ecutwfc = 60.0\n ecutrho = 480.0" POSCAR_${j}_${i}.pw   # cutoff



mkdir P_${j}_${i}

cp POSCAR_${j}_${i}.pw ./P_${j}_${i}/
cp OSZICAR ./P_${j}_${i}/OSZICAR_${j}_${i}_1
cd P_${j}_${i}

mpirun pw.x < POSCAR_${j}_${i}.pw > scf_${j}_${i}.out

E0=`grep "!    total energy              =" scf_${j}_${i}.out | awk '{printf"%.10f\n",$5}'`
a=13.605698066 # conversion factor between Ry and eV
E0ev=`echo "$E0*$a" | bc -l | awk '{printf"%.6f\n",$1}'`
sed -i "s/energy/${E0ev}/"  ./OSZICAR_${j}_${i}_1
cp OSZICAR_${j}_${i}_1 ../

cd ..


done

done


maelas -mode 3 -d -s ${smax} -n ${nmax} -b > output_derivation.dat

