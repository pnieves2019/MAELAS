#!/bin/bash

# load Python module in HPC facilities

ml Python/3.8.2-GCC-8.3.0-2.32-base

# folder where the POTCAR file is placed

path_potcar=/home/$USER/vasp/potentials/POTPAW_PBE_5.4/Fe_pv

# folder where we would like to run VASP calculations

path_cal=/scratch/work/project/open-17-14/$USER/maelas/example_Fe_qprod_k60

rm -f -r ${path_cal}/*

rm -f -r results/*

############################
# step 1: cell relaxation
############################

echo running step1 ...

mkdir results

mkdir results/Step1_cell_relaxation

cp POSCAR_Fe_bcc ./results/Step1_cell_relaxation
 
cp ${path_potcar}/POTCAR ./results/Step1_cell_relaxation

cd results/Step1_cell_relaxation

maelas -r -i POSCAR_Fe_bcc -k 60 -f ${path_cal}/step1 -a OPEN-19-14 -q qexp -t 1 > output.dat

qsub vasp_jsub_rlx > pid_job

###--- Check job status ---------------------------

read line < pid_job

echo pid_number=$line

qstat -f ${line} -x > job_summary

jobstate=`grep "job_state = " job_summary | awk '{print $3}'`

echo job_state = $jobstate

while [ $jobstate != "F" ]
 do
   rm job_summary
   qstat -f ${line} -x > job_summary
   jobstate=`grep "job_state = " job_summary | awk '{print $3}'`
   sleep 5s
 done

mv job_summary job_summary_step1


############################
# step 2: test MAE
############################

echo running step2 ... 

mkdir ../Step2_test_MAE

cd ../Step2_test_MAE

cp ${path_cal}/step1/CONTCAR ./POSCAR_Fe_bcc_rlx 

cp ${path_potcar}/POTCAR ./

maelas -m -i POSCAR_Fe_bcc_rlx -k 60 -s1 1 0 0 -s2 1 0 1 -f ${path_cal}/step2 -a OPEN-19-14 -q qprod > output.dat 

sed -i 's/qsub vasp_mae_jsub/qsub vasp_mae_jsub > ..\/\pid_job/g' vasp_mae

./vasp_mae

###--- Check job status -------------------------

read line < pid_job

echo pid_number=$line

qstat -f $line -x > job_summary

jobstate=`grep "job_state = " job_summary | awk '{print $3}'`

echo job_state = $jobstate

while [ "${jobstate}" != "F" ]
 do
   rm job_summary
   qstat -f $line -x > job_summary
   jobstate=`grep "job_state = " job_summary | awk '{print $3}'`
   sleep 5s
 done 

mv job_summary job_summary_step2

###--- copy OSZICAR files to check MAE -----

./vasp_mae_cp_oszicar

echo Test MAE:

echo Total energy of spin direction 1 is in file OSZICAR_0_0_1

echo Total energy of spin direction 2 is in file OSZICAR_0_0_2


############################
# step 3: generation VASP input files
############################

echo running step3 ...

mkdir ../Step3_generation_VASP_input_files

cd ../Step3_generation_VASP_input_files

cp ${path_cal}/step1/CONTCAR ./POSCAR_Fe_bcc_rlx

cp ${path_potcar}/POTCAR ./

maelas -g -i POSCAR_Fe_bcc_rlx -k 60 -s 0.01 -n 7 -f ${path_cal}/step3 -a OPEN-19-14 -q qprod > output.dat


############################
# step 4: run VASP jobs
############################

echo running step4 ...

sed -i 's/qsub vasp_jsub/qsub vasp_jsub > ..\/\pid_job_${i}_${j}/g' vasp_maelas

./vasp_maelas 

###--- Check job status --------------------

ndist=`grep "Number of distorted states for each magnetostriction mode =" output.dat | awk '{print $10}'`

nlamb=`grep "Number of anisotropic magnestostriction coefficients =" output.dat | awk '{print $7}'`

for ((ii=1;ii<=nlamb;ii++))
do
  for ((jj=1;jj<=ndist;jj++))
  do
    read a < pid_job_${ii}_${jj}

    jobid+=($a)

    echo pid_number_${ii}_${jj}=$a
  done
done

count=0
dim0=$(( ${#jobid[@]} - 1 ))

until [ $count -eq ${#jobid[@]} ]
do
    count=0
    for ((k=0;k<=dim0;k++))
    do       
      qstat -f ${jobid[$k]} -x > job_summary
      jobstate0="`grep "job_state = " job_summary | awk '{print $3}'`"
      
      if [ "${jobstate0}" = "F" ]
      then
        count=$(( $count + 1 ))       
      fi
      sleep 5s
      rm job_summary
    done
done


############################
# step 5a: Derivation of anisotropic magnetostriction coefficients
############################

echo running step5a ... 

mkdir ../Step5a_derivation_magnetostriction_coefficients

cp vasp_cp_oszicar ../Step5a_derivation_magnetostriction_coefficients

cp POSCAR* ../Step5a_derivation_magnetostriction_coefficients

cd ../Step5a_derivation_magnetostriction_coefficients

./vasp_cp_oszicar 

maelas -d -i POSCAR_Fe_bcc_rlx -n 7 > output.dat



############################
# step 5b: Derivation of magnetoelastic constants
############################

echo running step5b ...

mkdir ../Step5b_derivation_magnetoelastic_constants

cp OSZICAR* ../Step5b_derivation_magnetoelastic_constants

cp POSCAR* ../Step5b_derivation_magnetoelastic_constants

cd ../Step5b_derivation_magnetoelastic_constants

cp ../../ELADAT ./

maelas -d -i POSCAR_Fe_bcc_rlx -n 7 -b -e ELADAT > output.dat




echo program finished


