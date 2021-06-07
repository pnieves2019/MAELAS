#!/bin/bash


inputlmp=in.Fe    # input file name for LAMMPS

exelmp=lmp_serial  # name of the executable file for LAMMPS


i=1

for j in {1..9}
 
do
	

cp ${inputlmp} in.lmp_${i}_1_0

cp ${inputlmp} in.lmp_${i}_2_0

sed -i "s/AA/$i/"  ./in.lmp_${i}_1_0
sed -i "s/BB/$j/"  ./in.lmp_${i}_1_0
sed -i "s/CC/1/"  ./in.lmp_${i}_1_0

mx=0.0
my=0.0
mz=1.0

sed -i "s/spinx/$mx/"  ./in.lmp_${i}_1_0
sed -i "s/spiny/$my/"  ./in.lmp_${i}_1_0
sed -i "s/spinz/$mz/"  ./in.lmp_${i}_1_0


mx=1.0
my=0.0
mz=0.0


sed -i "s/spinx/$mx/"  ./in.lmp_${i}_2_0
sed -i "s/spiny/$my/"  ./in.lmp_${i}_2_0
sed -i "s/spinz/$mz/"  ./in.lmp_${i}_2_0

sed -i "s/AA/$i/"  ./in.lmp_${i}_2_0
sed -i "s/BB/$j/"  ./in.lmp_${i}_2_0
sed -i "s/CC/2/"  ./in.lmp_${i}_2_0


${exelmp} < in.lmp_${i}_1_0 > results_${i}_${j}_1
${exelmp} < in.lmp_${i}_2_0 > results_${i}_${j}_2


ene_1=`grep "energy" energy_${i}_${j}_1.dat | awk '{print $2}'`
ene_2=`grep "energy" energy_${i}_${j}_2.dat | awk '{print $2}'`

cp OSZICAR OSZICAR_${i}_${j}_1
cp OSZICAR OSZICAR_${i}_${j}_2

sed -i "s/energy_lmp/$ene_1/"  ./OSZICAR_${i}_${j}_1
sed -i "s/energy_lmp/$ene_2/"  ./OSZICAR_${i}_${j}_2


done


#################################### 111

i=2

	for j in {1..9}

	do

	

		cp ${inputlmp} in.lmp_${i}_1_0

		cp ${inputlmp} in.lmp_${i}_2_0

		sed -i "s/AA/$i/"  ./in.lmp_${i}_1_0
		sed -i "s/BB/$j/"  ./in.lmp_${i}_1_0
		sed -i "s/CC/1/"  ./in.lmp_${i}_1_0


		sed -i "s/AA/$i/"  ./in.lmp_${i}_2_0
		sed -i "s/BB/$j/"  ./in.lmp_${i}_2_0
		sed -i "s/CC/2/"  ./in.lmp_${i}_2_0

	
		if [ $j == 5 ]; 
		then
			mx=1.0
			my=1.0
			mz=1.0

			sed -i "s/spinx/$mx/"  ./in.lmp_${i}_1_0
			sed -i "s/spiny/$my/"  ./in.lmp_${i}_1_0
			sed -i "s/spinz/$mz/"  ./in.lmp_${i}_1_0


			mx=1.0
			my=0.0
			mz=-1.0
			
			sed -i "s/spinx/$mx/"  ./in.lmp_${i}_2_0
			sed -i "s/spiny/$my/"  ./in.lmp_${i}_2_0
			sed -i "s/spinz/$mz/"  ./in.lmp_${i}_2_0
		else

		ax=`grep "xlo" POSCAR_${i}_${j}.lmp | awk '{print $2}'`
		ay=0.0
		az=0.0
		bx=`grep "xy" POSCAR_${i}_${j}.lmp | awk '{print $1}'`
		by=`grep "ylo" POSCAR_${i}_${j}.lmp | awk '{print $2}'`
		bz=0.0
		cx=`grep "xy" POSCAR_${i}_${j}.lmp | awk '{print $2}'`
		cy=`grep "xy" POSCAR_${i}_${j}.lmp | awk '{print $3}'`
		cz=`grep "zlo" POSCAR_${i}_${j}.lmp | awk '{print $2}'`


		mx=`echo "scale=8; $ax+$bx+$cx" | bc` 
		my=`echo "scale=8; $ay+$by+$cy" | bc` 
		mz=`echo "scale=8; $az+$bz+$cz" | bc`

		sed -i "s/spinx/$mx/"  ./in.lmp_${i}_1_0
		sed -i "s/spiny/$my/"  ./in.lmp_${i}_1_0
		sed -i "s/spinz/$mz/"  ./in.lmp_${i}_1_0



		echo $j
		echo "new spin direction for trigonal deformation: mx, my, mz"
		echo $mx $my $mz
		
		ccx=`echo "scale=8; -1.0*$cx" | bc`
		ccy=`echo "scale=8; -1.0*$cy" | bc`
		ccz=`echo "scale=8; -1.0*$cz" | bc`

		mx=`echo "scale=8; $ax+$ccx" | bc` 
		my=`echo "scale=8; $ay+$ccy" | bc` 
		mz=`echo "scale=8; $az+$ccz" | bc` 

		echo $j
		echo "mx, my, mz"
		echo $mx $my $mz

		sed -i "s/spinx/$mx/"  ./in.lmp_${i}_2_0
		sed -i "s/spiny/$my/"  ./in.lmp_${i}_2_0
		sed -i "s/spinz/$mz/"  ./in.lmp_${i}_2_0
		
		fi


		${exelmp} < in.lmp_${i}_1_0 > results_${i}_${j}_1
		${exelmp} < in.lmp_${i}_2_0 > results_${i}_${j}_2


		ene_1=`grep "energy" energy_${i}_${j}_1.dat | awk '{print $2}'`
		ene_2=`grep "energy" energy_${i}_${j}_2.dat | awk '{print $2}'`

		cp OSZICAR OSZICAR_${i}_${j}_1
		cp OSZICAR OSZICAR_${i}_${j}_2

		sed -i "s/energy_lmp/$ene_1/"  ./OSZICAR_${i}_${j}_1
		sed -i "s/energy_lmp/$ene_2/"  ./OSZICAR_${i}_${j}_2

	done




maelas -d -mode 1 -i POSCAR_rlx -n 9 -b -e ELADAT
