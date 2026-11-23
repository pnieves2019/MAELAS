#!/bin/bash

 maelas -mode 2 -g -n 5 -nc > results_generate.out
 
 lmp_file="in.Fe"
 
 s1xb1=1.0
 s1yb1=0.0
 s1zb1=0.0
 
 s2xb1=1.0
 s2yb1=1.0
 s2zb1=0.0
 
 s1xb2=1.0
 s1yb2=1.0
 s1zb2=0.0
 
 s2xb2=-1.0
 s2yb2=1.0
 s2zb2=0.0
 
 for i in {1..2}
 do
 for j in {1..5}
 do
 for k in {1..2}
 do
 
 cp ${lmp_file} ${lmp_file}_${i}_${j}_${k}
 cp OSZICAR OSZICAR_${i}_${j}_${k}
 
 if [ $i == 1 ]
 then
 line1=$(sed -n 3p POSCAR_${i}_${j})
 line2=$(sed -n 4p POSCAR_${i}_${j})
 line3=$(sed -n 5p POSCAR_${i}_${j})
 else
 line1=$(sed -n 3p POSCAR_2_3)
 line2=$(sed -n 4p POSCAR_2_3)
 line3=$(sed -n 5p POSCAR_2_3)
 fi
 
 sed -i "s/lmpline1/${line1}/" ${lmp_file}_${i}_${j}_${k}
 sed -i "s/lmpline2/${line2}/" ${lmp_file}_${i}_${j}_${k}
 sed -i "s/lmpline3/${line3}/" ${lmp_file}_${i}_${j}_${k}

 if [ $i == 1 ] && [ $k == 1 ]
 then
 sed -i "s/sssx/${s1xb1}/" ${lmp_file}_${i}_${j}_${k}
 sed -i "s/sssy/${s1yb1}/" ${lmp_file}_${i}_${j}_${k}
 sed -i "s/sssz/${s1zb1}/" ${lmp_file}_${i}_${j}_${k}
 fi
 
 if [ $i == 1 ] && [ $k==2 ]
 then
 sed -i "s/sssx/${s2xb1}/" ${lmp_file}_${i}_${j}_${k}
 sed -i "s/sssy/${s2yb1}/" ${lmp_file}_${i}_${j}_${k}
 sed -i "s/sssz/${s2zb1}/" ${lmp_file}_${i}_${j}_${k}
 fi
 
 if [ $i == 2 ] && [ $k == 1 ]
 then
 sed -i "s/sssx/${s1xb2}/" ${lmp_file}_${i}_${j}_${k}
 sed -i "s/sssy/${s1yb2}/" ${lmp_file}_${i}_${j}_${k}
 sed -i "s/sssz/${s1zb2}/" ${lmp_file}_${i}_${j}_${k}
 if [ $j == 1 ]
 then
yes n | atomsk POSCAR_${i}_${j} lmp
 xytilt1b=$(sed -n 9p POSCAR_${i}_${j}.lmp | awk '{print $1}')
 xytilt1=$(echo "2.0*${xytilt1b}" | bc)
 sed -i "s|create_box	1 whole|create_box	1 whole \n change_box all triclinic \n change_box all xy final ${xytilt1} units box|" ${lmp_file}_${i}_${j}_${k}
 fi
 if [ $j == 2 ]
 then
yes n | atomsk POSCAR_${i}_${j} lmp
 xytilt2b=$(sed -n 9p POSCAR_${i}_${j}.lmp | awk '{print $1}')
 xytilt2=$(echo "2.0*${xytilt2b}" | bc)
 sed -i "s|create_box	1 whole|create_box	1 whole \n change_box all triclinic \n change_box all xy final ${xytilt2} units box|" ${lmp_file}_${i}_${j}_${k}
 fi
 if [ $j == 4 ]
 then
yes n | atomsk POSCAR_${i}_${j} lmp
 xytilt4b=$(sed -n 9p POSCAR_${i}_${j}.lmp | awk '{print $1}')
 xytilt4=$(echo "2.0*${xytilt4b}" | bc)
 sed -i "s|create_box	1 whole|create_box	1 whole \n change_box all triclinic \n change_box all xy final ${xytilt4} units box|" ${lmp_file}_${i}_${j}_${k}
 fi
 if [ $j == 5 ]
 then
yes n | atomsk POSCAR_${i}_${j} lmp
 xytilt5b=$(sed -n 9p POSCAR_${i}_${j}.lmp | awk '{print $1}')
 xytilt5=$(echo "2.0*${xytilt5b}" | bc)
 sed -i "s|create_box	1 whole|create_box	1 whole \n change_box all triclinic \n change_box all xy final ${xytilt5} units box|" ${lmp_file}_${i}_${j}_${k}
 fi
 fi
 
 if [ $i == 2 ] && [ $k == 2 ]
 then
 sed -i "s/sssx/${s2xb2}/" ${lmp_file}_${i}_${j}_${k}
 sed -i "s/sssy/${s2yb2}/" ${lmp_file}_${i}_${j}_${k}
 sed -i "s/sssz/${s2zb2}/" ${lmp_file}_${i}_${j}_${k}
 if [ $j == 1 ]
 then
 sed -i "s|create_box	1 whole|create_box	1 whole \n change_box all triclinic \n change_box all xy final ${xytilt1} units box|" ${lmp_file}_${i}_${j}_${k}
 fi
 if [ $j == 2 ]
 then
 sed -i "s|create_box	1 whole|create_box	1 whole \n change_box all triclinic \n change_box all xy final ${xytilt2} units box|" ${lmp_file}_${i}_${j}_${k}
 fi
 if [ $j == 4 ]
 then
 sed -i "s|create_box	1 whole|create_box	1 whole \n change_box all triclinic \n change_box all xy final ${xytilt4} units box|" ${lmp_file}_${i}_${j}_${k}
 fi
 if [ $j == 5 ]
 then
 sed -i "s|create_box	1 whole|create_box	1 whole \n change_box all triclinic \n change_box all xy final ${xytilt5} units box|" ${lmp_file}_${i}_${j}_${k}
 fi
 fi
 
 lmp_serial < ${lmp_file}_${i}_${j}_${k}
 
 ene=$(sed -n 1p energy_0.dat)
 
 
 
 sed -i "s/energy_lmp/${ene}/"  OSZICAR_${i}_${j}_${k}
 
 
 
 done
 done
 done
 
 maelas -mode 2 -d -n 5 -nc > results_derive.out
 
