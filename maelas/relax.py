from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from   maelas.data import SymmetryData

import os
import stat

class Relaxation:
    delec_list = ['Sc', 'Y', 'Ti', 'Zr', 'Hf', 'V', 'Nb', 'Ta', 'Cr', 'Mo', 'W', 'Mn', 'Tc', 'Re', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'Hg', 'Au', 'Ir', 'Pt', 'Os']
    felec_list = ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu','Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'U', 'Ac', 'Th', 'Pa', 'Np', 'Pu', 'Am']

    def __init__(self,args):
        self.lmax         = 2
        self.inc_rlx_list = []
        self.args         = args
        self.symData      = SymmetryData()

    def poscar(self):
        """   generating poscar for relaxation calculations   """
        print('--------------------------------------------------------------------------------------------------------')
        print("Generation of VASP files for the cell relaxation:")
        print('--------------------------------------------------------------------------------------------------------')

        self.symData.structure = Structure.from_file(self.args.pos[0])
        sym1 = float(self.args.sympre[0])
        sym2 = float(self.args.symang[0])
        aa = SpacegroupAnalyzer(self.symData.structure,symprec=sym1, angle_tolerance=sym2)
        self.symData.space_group = aa.get_space_group_number()
        print("Space group number =", self.symData.space_group)
        spg = aa.get_space_group_symbol()
        print("Space group symbol =", str(spg))
        self.symData.number_of_species = len(self.symData.structure.species)
        print("Number of atoms = {}".format(len(self.symData.structure.species)))

        pos_name = "POSCAR"
        structure00 = Poscar(self.symData.structure)
        structure00.write_file(filename = pos_name,significant_figures=16)

        return self.symData

    def incar(self):
        """   generating INCAR file for cell relaxation   """
        for i in range(self.symData.number_of_species):
            for j in range(len(self.delec_list)):
                if str(self.symData.structure.species[i]) == str(self.delec_list[j]):
                    #print('Material contains a d-element =', str(structure2.species[i]))
                    self.lmax = 4

        for i in range(self.symData.number_of_species):
            for j in range(len(self.felec_list)):
                if str(self.symData.structure.species[i]) == str(self.felec_list[j]):
                    #print('Material contains a f-element =', str(structure2.species[i]))
                    self.lmax = 6

        self.inc_rlx_list = ['ISTART = 0\n', 'NSW = 40\n', 'ENCUT = 520\n','IBRION = 1\n', 'ISIF = 3\n', 'EDIFFG = -0.001\n', '# LDAU = .TRUE.\n', '# LDAUL =\n', '# LDAUU =\n', '# LDAUJ = \n', '# LDAUTYPE = 2\n', 'LCHARG = FALSE\n', 'LWAVE = FALSE\n', 'PREC  = Normal\n', 'EDIFF  = 1.e-06\n', 'NELM   = 100\n', 'NELMIN = 4\n', 'ISMEAR = 1\n', 'SIGMA  = 0.10\n', 'ISPIN  = 2\n', 'LMAXMIX = ', self.lmax, ' ! for d-elements increase LMAXMIX to 4, f-elements: LMAXMIX = 6\n']
        path_inc_rlx = 'INCAR'
        inc_rlx = open(path_inc_rlx,'w')
        for entry in self.inc_rlx_list:
            inc_rlx.write(str(entry))
        mom_rlx = 'MAGMOM = ' + str(self.symData.number_of_species) + '*5'
        inc_rlx.write(mom_rlx)
        inc_rlx.close()

    def kpoints(self):
        """   KPOINT file   """
        path_kp = 'KPOINTS'
        kp_file = open(path_kp,'w')
        kp_file.write('k-points\n')
        kp_file.write('0\n')
        kp_file.write('Auto\n')
        kp_file.write(str(self.args.kp[0]))
        kp_file.close()

    def scripts(self):
        """   bash script to run vasp: vasp_jsub_rlx   """

        path_vasp_jsub = 'vasp_jsub_rlx'
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
        vasp_jsub.write('#PBS -N job_rlx\n')
        vasp_jsub.write('#PBS -j oe\n')
        vasp_jsub.write('#PBS -S /bin/bash\n')
        vasp_jsub.write('\n')
        vasp_jsub.write('cd ${PBS_O_WORKDIR}\n')
        vasp_jsub.write('SCRDIR=')
        vasp_jsub.write(str(self.args.vasp_fold[0]))
        vasp_jsub.write('\n')
        vasp_jsub.write('mkdir -p $SCRDIR\n')
        vasp_jsub.write('cd $SCRDIR || exit\n')
        vasp_jsub.write('cp -f -r $PBS_O_WORKDIR/* .\n')
        vasp_jsub.write('ml purge\n')
        vasp_jsub.write('ml ')
        vasp_jsub.write(str(self.args.load_module[0]))
        vasp_jsub.write('\n')
        vasp_jsub.write(str(self.args.mpi[0]))
        vasp_jsub.write(' -np ')
        vasp_jsub.write(str(self.args.core[0]))
        vasp_jsub.write(' vasp_std > vasp.out\n')

        vasp_jsub.write('exit\n')
        vasp_jsub.close()

        st = os.stat(path_vasp_jsub)
        os.chmod(path_vasp_jsub, st.st_mode | stat.S_IEXEC)
