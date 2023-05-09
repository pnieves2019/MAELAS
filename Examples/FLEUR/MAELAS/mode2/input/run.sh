#!/bin/bash

smax=0.018  #maximum applied strain

nmax=7 #number of deformed structures

ml Python/3.8.2-GCCcore-9.3.0

maelas -mode 2 -g -s ${smax} -n ${nmax} > output_generation.dat  # generation of deformed unit cells with -mode 2 of MAELAS

SG=`grep "Space group number =" output_generation.dat | awk '{printf"%i",$5}'`

#declare -A mx
#declare -A my
#declare -A mz
declare -A theta
declare -A phi


# We set spin directions (-mode 2 of MAELAS) according to the crsytal symmetry

if [[ $SG -le 230 ]] && [[ $SG -ge 207 ]]; then
	nb=2 #number of anisotropic magnetoleastic constants
	#mx[1,1]=1
	#my[1,1]=0
	#mz[1,1]=0
	theta[1,1]=1.57079632679
	phi[1,1]=0.0
	#mx[1,2]=1
	#my[1,2]=1
	#mz[1,2]=0
	theta[1,2]=1.57079632679
	phi[1,2]=0.785398163397
	#mx[2,1]=1
	#my[2,1]=1
	#mz[2,1]=0
	theta[2,1]=1.57079632679
	phi[2,1]=0.785398163397
	#mx[2,2]=-1
	#my[2,2]=1
	#mz[2,2]=0
	theta[2,2]=1.57079632679
	phi[2,2]=2.35619449019
elif [[ $SG -le 194 ]] && [[ $SG -ge 177 ]]; then
	nb=4
	#mx[1,1]=0
	#my[1,1]=0
	#mz[1,1]=1
	theta[1,1]=0.0
	phi[1,1]=0.0
	#mx[1,2]=1
	#my[1,2]=1
	#mz[1,2]=0
	theta[1,2]=1.57079632679
	phi[1,2]=0.785398163397
	#mx[2,1]=0
	#my[2,1]=0
	#mz[2,1]=1
	theta[2,1]=0.0
	phi[2,1]=0.0
	#mx[2,2]=1
	#my[2,2]=1
	#mz[2,2]=0
	theta[2,2]=1.57079632679
	phi[2,2]=0.785398163397
	#mx[3,1]=1
	#my[3,1]=0
	#mz[3,1]=0
	theta[3,1]=1.57079632679
	phi[3,1]=0.0
	#mx[3,2]=0
	#my[3,2]=1
	#mz[3,2]=0
	theta[3,2]=1.57079632679
	phi[3,2]=1.57079632679
	#mx[4,1]=1
	#my[4,1]=0
	#mz[4,1]=1
	theta[4,1]=0.785398163397
	phi[4,1]=0.0
	#mx[4,2]=-1
	#my[4,2]=0
	#mz[4,2]=1
	theta[4,2]=0.785398163397
	phi[4,2]=3.1415926535898
elif [[ $SG -le 167 ]] && [[ $SG -ge 149 ]]; then
	nb=6
	#mx[1,1]=0
	#my[1,1]=0
	#mz[1,1]=1
	theta[1,1]=0.0
	phi[1,1]=0.0
	#mx[1,2]=1
	#my[1,2]=1
	#mz[1,2]=0
	theta[1,2]=1.57079632679
	phi[1,2]=0.785398163397
	#mx[2,1]=0
	#my[2,1]=0
	#mz[2,1]=1
	theta[2,1]=0.0
	phi[2,1]=0.0
	#mx[2,2]=1
	#my[2,2]=1
	#mz[2,2]=0
	theta[2,2]=1.57079632679
	phi[2,2]=0.785398163397
	#mx[3,1]=1
	#my[3,1]=0
	#mz[3,1]=0
	theta[3,1]=1.57079632679
	phi[3,1]=0.0
	#mx[3,2]=0
	#my[3,2]=1
	#mz[3,2]=0
	theta[3,2]=1.57079632679
	phi[3,2]=1.57079632679
	#mx[4,1]=1
	#my[4,1]=0
	#mz[4,1]=1
	theta[4,1]=0.785398163397
	phi[4,1]=0.0
	#mx[4,2]=-1
	#my[4,2]=0
	#mz[4,2]=1
	theta[4,2]=0.785398163397
	phi[4,2]=3.1415926535898
	#mx[5,1]=1
	#my[5,1]=1
	#mz[5,1]=0
	theta[5,1]=1.57079632679
	phi[5,1]=0.785398163397
	#mx[5,2]=-1
	#my[5,2]=1
	#mz[5,2]=0
	theta[5,2]=1.57079632679
	phi[5,2]=2.35619449019
	#mx[6,1]=1
	#my[6,1]=0
	#mz[6,1]=1
	theta[6,1]=0.785398163397
	phi[6,1]=0.0
	#mx[6,2]=-1
	#my[6,2]=0
	#mz[6,2]=1
	theta[6,2]=0.785398163397
	phi[6,2]=3.1415926535898
elif [[ $SG -le 142 ]] && [[ $SG -ge 89 ]]; then
	nb=5
	#mx[1,1]=0
	#my[1,1]=0
	#mz[1,1]=1
	theta[1,1]=0.0
	phi[1,1]=0.0
	#mx[1,2]=1
	#my[1,2]=1
	#mz[1,2]=0
	theta[1,2]=1.57079632679
	phi[1,2]=0.785398163397
	#mx[2,1]=0
	#my[2,1]=0
	#mz[2,1]=1
	theta[2,1]=0.0
	phi[2,1]=0.0
	#mx[2,2]=1
	#my[2,2]=1
	#mz[2,2]=0
	theta[2,2]=1.57079632679
	phi[2,2]=0.785398163397
	#mx[3,1]=1
	#my[3,1]=0
	#mz[3,1]=0
	theta[3,1]=1.57079632679
	phi[3,1]=0.0
	#mx[3,2]=0
	#my[3,2]=1
	#mz[3,2]=0
	theta[3,2]=1.57079632679
	phi[3,2]=1.57079632679
	#mx[4,1]=1
	#my[4,1]=0
	#mz[4,1]=1
	theta[4,1]=0.785398163397
	phi[4,1]=0.0
	#mx[4,2]=-1
	#my[4,2]=0
	#mz[4,2]=1
	theta[4,2]=0.785398163397
	phi[4,2]=3.1415926535898
	#mx[5,1]=1
	#my[5,1]=1
	#mz[5,1]=0
	theta[5,1]=1.57079632679
	phi[5,1]=0.785398163397
	#mx[5,2]=-1
	#my[5,2]=1
	#mz[5,2]=0
	theta[5,2]=1.57079632679
	phi[5,2]=2.35619449019
elif [[ $SG -le 74 ]] && [[ $SG -ge 16 ]]; then
	nb=9
	#mx[1,1]=1
	#my[1,1]=0
	#mz[1,1]=0
	theta[1,1]=1.5707963267949
	phi[1,1]=0.0
	#mx[1,2]=0
	#my[1,2]=0
	#mz[1,2]=1
	theta[1,2]=0.0
	phi[1,2]=0.0
	#mx[2,1]=0
	#my[2,1]=1
	#mz[2,1]=0
	theta[2,1]=1.5707963267949
	phi[2,1]=1.5707963267949
	#mx[2,2]=0
	#my[2,2]=0
	#mz[2,2]=1
	theta[2,2]=0.0
	phi[2,2]=0.0
	#mx[3,1]=1
	#my[3,1]=0
	#mz[3,1]=0
	theta[3,1]=1.5707963267949
	phi[3,1]=0.0
	#mx[3,2]=0
	#my[3,2]=0
	#mz[3,2]=1
	theta[3,2]=0.0
	phi[3,2]=0.0
	#mx[4,1]=0
	#my[4,1]=1
	#mz[4,1]=0
	theta[4,1]=1.5707963267949
	phi[4,1]=1.5707963267949
	#mx[4,2]=0
	#my[4,2]=0
	#mz[4,2]=1
	theta[4,2]=0.0
	phi[4,2]=0.0
	#mx[5,1]=1
	#my[5,1]=0
	#mz[5,1]=0
	theta[5,1]=1.5707963267949
	phi[5,1]=0.0
	#mx[5,2]=0
	#my[5,2]=0
	#mz[5,2]=1
	theta[5,2]=0.0
	phi[5,2]=0.0
	#mx[6,1]=0
	#my[6,1]=1
	#mz[6,1]=0
	theta[6,1]=1.5707963267949
	phi[6,1]=1.5707963267949
	#mx[6,2]=0
	#my[6,2]=0
	#mz[6,2]=1
	theta[6,2]=0.0
	phi[6,2]=0.0
	#mx[7,1]=1
	#my[7,1]=1
	#mz[7,1]=0
	theta[7,1]=1.5707963267949
	phi[7,1]=0.78539816339745
	#mx[7,2]=0
	#my[7,2]=0
	#mz[7,2]=1
	theta[7,2]=0.0
	phi[7,2]=0.0
	#mx[8,1]=1
	#my[8,1]=0
	#mz[8,1]=1
	theta[8,1]=0.78539816339745
	phi[8,1]=0.0
	#mx[8,2]=0
	#my[8,2]=0
	#mz[8,2]=1
	theta[8,2]=0.0
	phi[8,2]=0.0
	#mx[9,1]=0
	#my[9,1]=1
	#mz[9,1]=1
        theta[9,1]=0.78539816339745
	phi[9,1]=1.5707963267949
	#mx[9,2]=0
	#my[9,2]=0
	#mz[9,2]=1
	theta[9,2]=0.0
	phi[9,2]=0.0
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


mkdir P_${j}_${i}_${k}

cp POSCAR_${j}_${i} ./P_${j}_${i}_${k}/POSCAR
cp jsub ./P_${j}_${i}_${k}/
cd P_${j}_${i}_${k}


sed -i "s/AA/$j/"  ./jsub
sed -i "s/BB/$i/"  ./jsub
sed -i "s/CC/$k/"  ./jsub

ml purge
ml intel
ml libxml2
ml ifort
ml Python/3.8.2-GCCcore-9.3.0

vasp2cif POSCAR  # conversion from VASP to cif

cif2cell -f POSCAR.cif -p fleur -o inp_file  # conversion from cif to FLEUR

sed -i "/oldfleur=f/ c\&input cartesian=f oldfleur=f film=f /" inp_file

echo -e "\n&comp gmax=12.0 kmax=5.0 /" > gmax0.dat  # we change plane-wave cutoff

echo -e "\n&soc 0.10 0.37 /" > soc0.dat  # we include spin-orbit coupling

cat inp_file gmax0.dat soc0.dat >> inp_file0  

mv inp_file0 inp_file

inpgen -f inp_file -kpt testset#gamma@grid=20,20,20  # we define new k-point mesh


sed -i 's/itmax="15"/itmax="60"/' inp.xml # we increase maximum number of iteration in self-consistent field (scf) cycle

sed -i 's/kPointListSelection listName="default"/kPointListSelection listName="testset"/' inp.xml  # we set our new k-point mesh

sed -i '/soc l_soc=/ c\<soc l_soc="T" theta="thetasoc" phi="phisoc" spav="F"/>' inp.xml   

sed -i "s/thetasoc/${theta[${j},${k}]}/"  ./inp.xml  # We set spin direction (theta angle)

sed -i "s/phisoc/${phi[${j},${k}]}/"  ./inp.xml  # We set spin direction (phi angle)

qsub jsub   # run fleur

cd ..

done

done

done


