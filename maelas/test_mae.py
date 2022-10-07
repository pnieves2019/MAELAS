import os
import stat

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar

from maelas.data import SymmetryData

# _            _     __  __    _    _____ 
#| |_ ___  ___| |_  |  \/  |  / \  | ____|
#| __/ _ \/ __| __| | |\/| | / _ \ |   _|  
#| ||  __/\__ \ |_  | |  | |/ ___ \|  |___ 
# \__\___||___/\__| |_|  |_|_/   \_\_____/
#                                         
class TestMAE:
    delec_list = ['Sc', 'Y', 'Ti', 'Zr', 'Hf', 'V', 'Nb', 'Ta', 'Cr', 'Mo', 'W', 'Mn', 'Tc', 'Re', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'Hg', 'Au', 'Ir', 'Pt', 'Os']
    felec_list = ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu','Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'U', 'Ac', 'Th', 'Pa', 'Np', 'Pu', 'Am']
    def __init__(self,args):
        """   Generation of VASP files to test MAE   """
        self.lmax         =  2
        self.inc_std_list = []
        self.inc_ncl_list = []
        self.args         = args
        self.symData      = SymmetryData()

    def poscar(self):
        print('--------------------------------------------------------------------------------------------------------')
        print("Generation of VASP files to test MAE")
        print('--------------------------------------------------------------------------------------------------------')
    
        self.symData.structure = Structure.from_file(self.args.pos[0])
        sym1 = float(self.args.sympre[0])
        sym2 = float(self.args.symang[0])
        aa = SpacegroupAnalyzer(self.symData.structure,symprec=sym1, angle_tolerance=sym2)
        self.symData.space_group = aa.get_space_group_number()
        print("Space group number =", self.symData.space_group)
        self.symData.point_group = aa.get_space_group_symbol()
        print("Space group symbol =", str(self.symData.point_group))
        self.symData.number_of_species = len(self.symData.structure.species)
        print("Number of atoms = {}".format(len(self.symData.structure.species)))
    
        print("First spin direction =", self.args.spin1[0],self.args.spin1[1],self.args.spin1[2])
        print("Second spin direction =", self.args.spin2[0],self.args.spin2[1],self.args.spin2[2])
    
        # POSCAR file
        pos_name = "POSCAR_0_0"
        self.symData.structure0 = Poscar(self.symData.structure)
        self.symData.structure0.write_file(filename = pos_name,significant_figures=16)
    
    def kpoints(self):
        """   KPOINT file   """
        path_kp = 'KPOINTS'
        kp_file = open(path_kp,'w')
        kp_file.write('k-points\n')
        kp_file.write('0\n')
        kp_file.write('Auto\n')
        kp_file.write(str(self.args.kp[0]))
        kp_file.close()
    
    def incar(self):
        """   INCAR file   """
        for i in range(self.symData.number_of_species):
            for j in range(len(self.delec_list)):
                if str(self.symData.structure.species[i]) == str(self.delec_list[j]):
                    #print('Material contains a d-element =', str(structure2.species[i]))
                    self.lmax = 4
                    break
            for j in range(len(self.felec_list)):
                if str(self.symData.structure.species[i]) == str(self.felec_list[j]):
                    #print('Material contains a f-element =', str(structure2.species[i]))
                    self.lmax = 6
                    break
    
        # INCAR file for collinear calculation
        self.inc_std_list = ['ISTART = 0\n', 'LORBIT = 11\n', 'ISYM = -1\n', 'ENCUT = 520\n','PREC  = Accurate\n', 'EDIFF  = 1.e-09\n', 'NELM   = 100\n', 'NELMIN = 4\n','# LDAU = .TRUE.\n', '# LDAUL =\n', '# LDAUU =\n', '# LDAUJ = \n', '# LDAUTYPE = 2\n', 'ADDGRID = TRUE\n', 'ISMEAR = -5\n', '# SIGMA  = 0.10\n', 'ISPIN  = 2\n', 'LMAXMIX = ', self.lmax, ' ! for d-elements increase LMAXMIX to 4, f-elements: LMAXMIX = 6\n', 'GGA_COMPAT = FALSE\n', 'LREAL = FALSE\n', 'LCHARG = TRUE\n', 'LWAVE = TRUE\n']
        path_inc_std = 'INCAR_std'
        inc_std = open(path_inc_std,'w')
        for j in range(len(self.inc_std_list)):
            inc_std.write(str(self.inc_std_list[j]))
        mom_std = 'MAGMOM = ' + str(self.symData.number_of_species) + '*5'
        inc_std.write(mom_std)
        inc_std.close()
    
        # INCAR file for non-collinear calculation spin 1
        self.inc_ncl_list = ['LORBIT = 11\n', 'ISYM = -1\n', 'ENCUT = 520\n', 'PREC  = Accurate\n', 'EDIFF  = 1.e-09\n', 'NELM   = 100\n', 'NELMIN = 4\n','# LDAU = .TRUE.\n', '# LDAUL =\n', '# LDAUU =\n', '# LDAUJ = \n', '# LDAUTYPE = 2\n' ,'ADDGRID = TRUE\n', 'ISMEAR = -5\n', '# SIGMA  = 0.10\n', 'ISPIN  = 2\n', 'LMAXMIX = ', self.lmax, ' ! for d-elements increase LMAXMIX to 4, f-elements: LMAXMIX = 6\n', 'GGA_COMPAT = FALSE\n', 'LREAL = FALSE\n', 'ICHARG = 11\n', 'LCHARG = TRUE\n', 'LWAVE = TRUE\n', 'LORBMOM = TRUE\n', 'LSORBIT = TRUE\n', 'NBANDS = nbands ! 2 * number of bands of collinear run\n' ,'MAGMOM = ']
    
        for i in range(self.symData.number_of_species):
            self.inc_ncl_list += ['0 0 4 ']
        self.inc_ncl_list += ['\n']
    
        path_inc_ncl_1_1 = 'INCAR_0_1'
        inc_ncl_1_1 = open(path_inc_ncl_1_1,'w')
        self.inc_ncl_list_1_1 = self.inc_ncl_list[:]
        self.inc_ncl_list_1_1 += ['SAXIS = ',self.args.spin1[0],' ', self.args.spin1[1],' ', self.args.spin1[2]]
        self.inc_ncl_list_1_1 += ['\n']
    
        for entry in self.inc_ncl_list_1_1:
            inc_ncl_1_1.write(str(entry))
    
        inc_ncl_1_1.close()
    
        # INCAR file for non-collinear calculation spin 2
        path_inc_ncl_1_1 = 'INCAR_0_2'
        inc_ncl_1_1 = open(path_inc_ncl_1_1,'w')
        self.inc_ncl_list_1_1 = self.inc_ncl_list[:]
        self.inc_ncl_list_1_1 += ['SAXIS = ',self.args.spin2[0],' ', self.args.spin2[1],' ', self.args.spin2[2]]
        self.inc_ncl_list_1_1 += ['\n']
    
        for entry in self.inc_ncl_list_1_1:
            inc_ncl_1_1.write(str(entry))
    
        inc_ncl_1_1.close()
    
    def scripts(self): 
        """   Generation of bash scripts for running vasp easily   """
       # bash script: vasp_mae
        nmag = 0
        self.args.ndist[0] = 0
    
        path_vasp_mag = 'vasp_mae'
        vasp_mag = open(path_vasp_mag,'w')
    
        vasp_mag.write('#! /bin/bash\n')
        vasp_mag.write('\n')
        vasp_mag.write('for i in {0..')
        vasp_mag.write(str(nmag))
        vasp_mag.write('}\n')
        vasp_mag.write('do\n')
        vasp_mag.write('for j in {0..')
        vasp_mag.write(str(self.args.ndist[0]))
        vasp_mag.write('}\n')
        vasp_mag.write('do\n')
        vasp_mag.write('mkdir P_${i}_${j}\n')
        vasp_mag.write('mkdir P_${i}_${j}/ncl_1\n')
        vasp_mag.write('mkdir P_${i}_${j}/ncl_2\n')
        vasp_mag.write('cp POSCAR_${i}_${j} ./P_${i}_${j}/POSCAR\n')
        vasp_mag.write('cp KPOINTS ./P_${i}_${j}\n')
        vasp_mag.write('cp POTCAR ./P_${i}_${j}\n')
        vasp_mag.write('cp INCAR_std ./P_${i}_${j}/INCAR\n')
        vasp_mag.write('cp vasp_mae_jsub ./P_${i}_${j}/\n')
        vasp_mag.write('sed -i "s/AA/$i/"  ./P_${i}_${j}/vasp_mae_jsub\n')
        vasp_mag.write('sed -i "s/BB/$j/"  ./P_${i}_${j}/vasp_mae_jsub\n')
        vasp_mag.write('cp POSCAR_${i}_${j} ./P_${i}_${j}/ncl_1/POSCAR\n')
        vasp_mag.write('cp KPOINTS ./P_${i}_${j}/ncl_1\n')
        vasp_mag.write('cp POTCAR ./P_${i}_${j}/ncl_1\n')
        vasp_mag.write('cp INCAR_${i}_1 ./P_${i}_${j}/ncl_1/INCAR\n')
        vasp_mag.write('cp POSCAR_${i}_${j} ./P_${i}_${j}/ncl_2/POSCAR\n')
        vasp_mag.write('cp KPOINTS ./P_${i}_${j}/ncl_2\n')
        vasp_mag.write('cp POTCAR ./P_${i}_${j}/ncl_2\n')
        vasp_mag.write('cp INCAR_${i}_2 ./P_${i}_${j}/ncl_2/INCAR\n')
        vasp_mag.write('cp vasp_mae_0 ./P_${i}_${j}/\n')
        vasp_mag.write('cd ./P_${i}_${j}/\n')
        vasp_mag.write('qsub vasp_mae_jsub\n')
        vasp_mag.write('cd ..\n')
        vasp_mag.write('done\n')
        vasp_mag.write('done\n')
        vasp_mag.write('\n')
    
        vasp_mag.close()
    
       # bash script: vasp_mae_jsub
    
        path_vasp_jsub = 'vasp_mae_jsub'
        vasp_jsub = open(path_vasp_jsub,'w')
    
        vasp_jsub.write('#!/bin/bash\n')
        vasp_jsub.write('#PBS -A ')
        vasp_jsub.write(str(self.args.p_id[0]))
        vasp_jsub.write('\n')
        vasp_jsub.write('#PBS -q ')
        vasp_jsub.write(str(self.args.queue[0]))
        vasp_jsub.write('\n')
        vasp_jsub.write('#PBS -l select=1:ncpus=')
        vasp_jsub.write(str(self.args.core[0]))
        vasp_jsub.write(':mpiprocs=')
        vasp_jsub.write(str(self.args.core[0]))
        vasp_jsub.write(':ompthreads=1\n')
        vasp_jsub.write('#PBS -l walltime=')
        vasp_jsub.write(str(self.args.time[0]))
        vasp_jsub.write(':00:00\n')
        vasp_jsub.write('#PBS -N job_AA_BB\n')
        vasp_jsub.write('#PBS -j oe\n')
        vasp_jsub.write('#PBS -S /bin/bash\n')
        vasp_jsub.write('\n')
        vasp_jsub.write('cd ${PBS_O_WORKDIR}\n')
        vasp_jsub.write('SCRDIR=')
        vasp_jsub.write(str(self.args.vasp_fold[0]))
        vasp_jsub.write('/P_AA_BB\n')
        vasp_jsub.write('\n')
        vasp_jsub.write('mkdir -p $SCRDIR\n')
        vasp_jsub.write('cd $SCRDIR || exit\n')
        vasp_jsub.write('cp -f -r $PBS_O_WORKDIR/* .\n')
        vasp_jsub.write('ml purge\n')
        vasp_jsub.write('ml ')
        vasp_jsub.write(str(self.args.load_module[0]))
        vasp_jsub.write('\n')
        vasp_jsub.write('./vasp_mae_0 >> log\n')
        vasp_jsub.write('exit\n')
    
        vasp_jsub.close()
    
       # bash script: vasp_mae_0
    
        path_vasp_0 = 'vasp_mae_0'
        vasp_0 = open(path_vasp_0,'w')
    
        vasp_0.write('#!/bin/bash\n')
        vasp_0.write(str(self.args.mpi[0]))
        vasp_0.write(' -np ')
        vasp_0.write(str(self.args.core[0]))
        vasp_0.write(' vasp_std > vasp.out\n')
        vasp_0.write('fold1=ncl_1\n')
        vasp_0.write('fold2=ncl_2\n')
        vasp_0.write('cp WAVECAR ./${fold1}/\n')
        vasp_0.write('cp CHGCAR  ./${fold1}/\n')
        vasp_0.write("nbands=`grep \"NBANDS\" OUTCAR | awk '{printf\"%d\",$15}'`\n")
        vasp_0.write("nbands=`echo \"2*$nbands\" | bc -l | awk '{printf\"%d\",$1}'`\n")
        vasp_0.write('cd ./${fold1}\n')
        vasp_0.write('sed -i "s/nbands/$nbands/" INCAR\n')
        vasp_0.write(str(self.args.mpi[0]))
        vasp_0.write(' -np ')
        vasp_0.write(str(self.args.core[0]))
        vasp_0.write(' vasp_ncl > vasp.out\n')
        vasp_0.write('cd ..\n')
        vasp_0.write('cp WAVECAR ./${fold2}/\n')
        vasp_0.write('cp CHGCAR  ./${fold2}/\n')
        vasp_0.write("nbands=`grep \"NBANDS\" OUTCAR | awk '{printf\"%d\",$15}'`\n")
        vasp_0.write("nbands=`echo \"2*$nbands\" | bc -l | awk '{printf\"%d\",$1}'`\n")
        vasp_0.write('cd ./${fold2}\n')
        vasp_0.write('sed -i "s/nbands/$nbands/" INCAR\n')
        vasp_0.write(str(self.args.mpi[0]))
        vasp_0.write(' -np ')
        vasp_0.write(str(self.args.core[0]))
        vasp_0.write(' vasp_ncl > vasp.out\n')
    
        vasp_0.close()
    
       # bash script: vasp_mae_cp_oszicar
    
        path_vasp_osz = 'vasp_mae_cp_oszicar'
        vasp_osz = open(path_vasp_osz,'w')
    
        vasp_osz.write('#!/bin/bash\n')
        vasp_osz.write('path_files=')
        vasp_osz.write(str(self.args.vasp_fold[0]))
        vasp_osz.write('\n')
        vasp_osz.write('\n')
        vasp_osz.write('for i in {0..')
        vasp_osz.write(str(nmag))
        vasp_osz.write('}\n')
        vasp_osz.write('do\n')
        vasp_osz.write('for j in {0..')
        vasp_osz.write(str(self.args.ndist[0]))
        vasp_osz.write('}\n')
        vasp_osz.write('do\n')
        vasp_osz.write('cp ${path_files}/P_${i}_$j/ncl_1/OSZICAR ./OSZICAR_${i}_${j}_1\n')
        vasp_osz.write('cp ${path_files}/P_${i}_$j/ncl_2/OSZICAR ./OSZICAR_${i}_${j}_2\n')
        vasp_osz.write('done\n')
        vasp_osz.write('done\n')
    
        vasp_osz.close()
    
        st = os.stat(path_vasp_osz)
        os.chmod(path_vasp_osz, st.st_mode | stat.S_IEXEC)
    
        st = os.stat(path_vasp_0)
        os.chmod(path_vasp_0, st.st_mode | stat.S_IEXEC)
    
        st = os.stat(path_vasp_jsub)
        os.chmod(path_vasp_jsub, st.st_mode | stat.S_IEXEC)
    
        st = os.stat(path_vasp_mag)
        os.chmod(path_vasp_mag, st.st_mode | stat.S_IEXEC)
