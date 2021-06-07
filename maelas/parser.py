import argparse

class MAELAS_Options:
    """   A class handling the commandline options   """
    def __init__(self):
        """   Constructor: input is the unpacked argv form system   """
        self.parser = argparse.ArgumentParser(description='MAELAS code v2.0')
        self.parser.add_argument('-mode', dest='mode', type=int, nargs=1, 
                            default=['2'], 
                            help='-mode 1: Scheme for the direct calculation of the magnetostrictive coefficients. -mode 2: Scheme for the direct calculation of the magnetoelastic constants. Mode 2 is more accurate for non-cubic symmetries (default: 2) ')
        self.parser.add_argument('-i', dest='pos', type=str, nargs=1,
                            default=['POSCAR'], 
                            help='Name of the initial non-distorted POSCAR file (default: POSCAR)')
        self.parser.add_argument('-n', dest='ndist', type=int, nargs=1, 
                            default=['7'], 
                            help='Number of distorted states for the direct calculation for each magnetostriction coefficient (-mode 1) or magentoelastic contant (-mode 2)(default: 7)')
        self.parser.add_argument('-s', dest='strain', type=float, nargs=1,
                            default=['0.01'],
                            help='Maximum value of the parameter s for the deformation gradient Fij(s) to generate the distorted POSCAR files (default: 0.01)')
        self.parser.add_argument('-k', dest='kp', type=int, nargs=1,
                            default=['60'],
                            help='VASP automatic k-point mesh generation to create the KPOINTS file (default: 60)')
        self.parser.add_argument('-g', dest='gen', action='store_true',
                            default=False,
                            help='Generation of required VASP files for the direct calculation of magnetostriction coefficients (-mode 1) or magnetoelastic constants (-mode 2). Notation of the generated output files: POSCAR_A_B (distorted cell where A=magnetostriction coefficient (-mode 1) or magnetoelastic constant (-mode 2), B=distorted cell), INCAR_A_C (non-collinear calculation where C=spin orientation case), INCAR_std (collinear calculation). How to run the VASP calculations: For each generated POSCAR_A_B one should run first a collinear calculation using INCAR_std and use the generated WAVECAR and CHGCAR files to run non-collinear calculations for each INCAR_A_C using the same POSCAR_A_B. It also generates bash scripts to run VASP calculations easily (vasp_maelas, vasp_jsub, vasp_0) and to get calculated OSZICAR_A_B_C files (vasp_cp_oszicar)' )
        self.parser.add_argument('-d', dest='der', action='store_true',
                            default=False,
                            help='Derivation of magnetostriction coefficients (-mode 1) or magnetoelastic contants (-mode 2) from the energy written in the OSZICAR files. WARNING!: OSZICAR files must be in the same folder where you run MAELAS using the notation OSZICAR_A_B_C obtained for POSCAR_A_B and INCAR_A_C. Distorted POSCAR files (POSCAR_A_B) must be in this folder too (jointly with the initial non-distorted POSCAR which should be specified using tag -i). Specify the number of distorted states to be considered in the calculation of magnetostriction coefficients using tag -n, and the same maximum applied strain (tag -s) used in the generation of VASP input files. Energy values extracted from OSZICAR_A_B_C files are shown in files ene_A_C.dat and fit_ene_A_C.png. The energy difference between the two spin configurations for each magnetostriction mode are shown in Figs. dE_A.png')
        self.parser.add_argument('-r', dest='rel', action='store_true',
                            default=False, help='Generation of required VASP files for the cell relaxation')
        self.parser.add_argument('-m', dest='mae', action='store_true',
                            default=False,
                            help='Generation of required VASP files to test MAE')
        self.parser.add_argument('-s1', dest='spin1', type=float, nargs=3,
                            default=['1','0','0'],
                            help='First spin direction to calculate MAE: s1x s1y s1z')
        self.parser.add_argument('-s2', dest='spin2', type=float, nargs=3,
                            default=['0','0','1'],
                            help='Second spin direction to calculate MAE: s2x s2y s2z')
        self.parser.add_argument('-b', dest='delas', action='store_true',
                            default=False,
                            help='Indirect calculation of the magnetoelastic constants from the calculated magnetostriction coefficients and provided elastic tensor (-mode 1) or indirect calculation of the magnetostrictive coefficients from the calculated magnetoelastic constants and provided elastic tensor (-mode 2). For this option the tag -d must be included as well as tag -e with the name of the elastic tensor file')
        self.parser.add_argument('-e', dest='elas', type=str, nargs=1,
                            default=['ELADAT'],
                            help='File with the elastic tensor data in the same format and units (GPa) as it is written by ELAS code (file ELADAT). You can check this format in the Examples folder')
        self.parser.add_argument('-sp', dest='sympre', type=float, nargs=1,
                            default=['0.01'],
                            help='Tolerance for symmetry finding (default: 0.01)')
        self.parser.add_argument('-sa', dest='symang', type=float, nargs=1,
                            default=['5.0'],
                            help='Angle tolerance for symmetry finding (default: 5.0)')
        self.parser.add_argument('-sg', dest='sg0', type=int, nargs=1,
                            default=['0'],
                            help='Space group number 1-230. If it is equal to 0, then it will be determined by a symmetry analysis (default: 0)')
        self.parser.add_argument('-c', dest='core', type=int, nargs=1,
                            default=['24'],
                            help='Number of cores for the VASP calculation (default: 24)')
        self.parser.add_argument('-t', dest='time', type=int, nargs=1,
                            default=['48'],
                            help='Number of maximum CPU hours for the VASP calculation (default: 48)')
        self.parser.add_argument('-f', dest='vasp_fold', type=str, nargs=1,
                            default=['/scratch'],
                            help='Folder where you will run VASP calculations (default: /scratch)')
        self.parser.add_argument('-mp', dest='mpi', type=str, nargs=1,
                            default=['mpiexec.hydra'],
                            help='Command for mpi run of VASP (default: mpiexec.hydra)')
        self.parser.add_argument('-a', dest='p_id', type=str, nargs=1,
                            default=['OPEN-X-X'],
                            help='Project id for running jobs in HPC facilities (default: OPEN-X-X)')
        self.parser.add_argument('-l', dest='load_module', type=str, nargs=1,
                            default=['VASP/5.4.4-intel-2017c-mkl=cluster'],
                            help='Module of VASP that should be loaded (default: VASP/5.4.4-intel-2017c-mkl=cluster)')
        self.parser.add_argument('-q', dest='queue', type=str, nargs=1,
                            default=['qprod'],
                            help='Type of queue to be used for VASP calculations in HPC facilities (default: qprod)')
    def __call__(self, *args,**kwargs):
        """  Parsing commandline argv. Omitting argv[0].   """
        return self.parser.parse_args(*args,**kwargs)
