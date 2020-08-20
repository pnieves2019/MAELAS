"""
    __  ______    ________    ___   _____
   /  |/  /   |  / ____/ /   /   | / ___/
  / /|_/ / /| | / __/ / /   / /| | \__ \ 
 / /  / / ___ |/ /___/ /___/ ___ |___/ / 
/_/  /_/_/  |_/_____/_____/_/  |_/____/  
                                         

usage: maelas [-h] [-i POS] [-n NDIST] [-s STRAIN] [-k KP] [-g] [-d] [-r] [-m]
              [-s1 SPIN1 SPIN1 SPIN1] [-s2 SPIN2 SPIN2 SPIN2] [-b] [-e ELAS]
              [-sp SYMPRE] [-sa SYMANG] [-sg SG0] [-c CORE] [-t TIME]
              [-f VASP_FOLD] [-mp MPI] [-a P_ID] [-l LOAD_MODULE] [-q QUEUE]

MAELAS code v1.0

optional arguments:
  -h, --help            show this help message and exit
  -i POS                Name of the initial non-distorted POSCAR file
                        (default: POSCAR)
  -n NDIST              Number of distorted states for each magnetostriction
                        mode (default: 7)
  -s STRAIN             Maximum value of the parameter epsilon for the strain
                        tensor to generate the distorted POSCAR files
                        (default: 0.01)
  -k KP                 VASP automatic k-point mesh generation to create the
                        KPOINTS file (default: 60)
  -g                    Generation of required VASP files for the calculation
                        of magnetostriction coefficients. Notation of the
                        generated output files: POSCAR_A_B (volume-conserving
                        distorted cell where A=magnetostriction mode,
                        B=distorted cell), INCAR_A_C (non-collinear
                        calculation where A=magnetostriction mode, C=spin
                        orientation case), INCAR_std (collinear calculation).
                        How to run the VASP calculations: For each generated
                        POSCAR_A_B one should run first a collinear
                        calculation using INCAR_std and use the generated
                        WAVECAR and CHGCAR files to run non-collinear
                        calculations for each INCAR_A_C using the same
                        POSCAR_A_B. It also generates bash scripts to run VASP
                        calculations easily (vasp_maelas, vasp_jsub, vasp_0)
                        and to get calculated OSZICAR_A_B_C files
                        (vasp_cp_oszicar)
  -d                    Derivation of magnetostriction coefficients from the
                        energy written in the OSZICAR files. WARNING!: OSZICAR
                        files must be in the same folder where you run MAELAS
                        using the notation OSZICAR_A_B_C obtained for
                        POSCAR_A_B and INCAR_A_C. Distorted POSCAR files
                        (POSCAR_A_B) must be in this folder too (jointly with
                        the initial non-distorted POSCAR which should be
                        specified using tag -i). Specify the number of
                        distorted states to be considered in the calculation
                        of magnetostriction coefficients using tag -n. Energy
                        values extracted from OSZICAR_A_B_C files are shown in
                        files ene_A_C.dat and fit_ene_A_C.png. The energy
                        difference between the two spin configurations for
                        each magnetostriction mode are shown in Figs. dE_A.png
  -r                    Generation of required VASP files for the cell
                        relaxation
  -m                    Generation of required VASP files to test MAE
  -s1 SPIN1 SPIN1 SPIN1
                        First spin direction to calculate MAE: s1x s1y s1z
  -s2 SPIN2 SPIN2 SPIN2
                        Second spin direction to calculate MAE: s2x s2y s2z
  -b                    Calculation of the magnetoelastic constants from the
                        calculated magnetostriction coefficients and provided
                        elastic tensor. For this option the tag -d must be
                        included as well as tag -e with the elastic tensor
                        file
  -e ELAS               File with the elastic tensor data in the same format
                        and units (GPa) as it is written by ELAS code (file
                        ELADAT). You can check this format in the Examples
                        folder
  -sp SYMPRE            Tolerance for symmetry finding (default: 0.01)
  -sa SYMANG            Angle tolerance for symmetry finding (default: 5.0)
  -sg SG0               Space group number 1-230. If it is equal to 0, then it
                        will be determined by a symmetry analysis (default: 0)
  -c CORE               Number of cores for the VASP calculation (default: 24)
  -t TIME               Number of maximum CPU hours for the VASP calculation
                        (default: 48)
  -f VASP_FOLD          Folder where you will run VASP calculations (default:
                        /scratch)
  -mp MPI               Command for mpi run of VASP (default: mpiexec.hydra)
  -a P_ID               Project id for running jobs in HPC facilities
                        (default: OPEN-X-X)
  -l LOAD_MODULE        Module of VASP that should be loaded (default:
                        VASP/5.4.4-intel-2017c-mkl=cluster)
  -q QUEUE              Type of queue to be used for VASP calculations in HPC
                        facilities (default: qprod)

"""
