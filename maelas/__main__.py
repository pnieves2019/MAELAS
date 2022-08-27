#!/bin/bash

import math
import os
import stat

import maelas.parser   as parser
import maelas.generate as generate
import maelas.relax    as relax
import maelas.test_mae as test_mae
import maelas.cubic_I_mode1 as cubic_I_mode1
import maelas.cubic_I_mode2 as cubic_I_mode2
import maelas.hex_I_tet_I_mode1 as hex_I_tet_I_mode1
import maelas.hex_I_tet_I_mode2 as hex_I_tet_I_mode2
import maelas.trig_I_mode1 as trig_I_mode1
import maelas.trig_I_mode2 as trig_I_mode2
import maelas.ort_mode1 as ort_mode1
import maelas.ort_mode2 as ort_mode2
import maelas.mode3 as mode3
from   maelas.data import SymmetryData

from pymatgen import Lattice, Structure
from pymatgen.transformations.standard_transformations import ConventionalCellTransformation,DeformStructureTransformation
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar
from pyfiglet import Figlet



##############################################

f = Figlet(font='slant')
print(f.renderText('MAELAS v3.0'))

# ____                _                                  _   _ _
#|  _ \ __ _ _ __ ___(_)_ __   __ _    ___ _ __ ___   __| | | (_)_ __   ___
#| |_) / _` | '__/ __| | '_ \ / _` |  / __| '_ ` _ \ / _` | | | | '_ \ / _ \
#|  __/ (_| | |  \__ \ | | | | (_| | | (__| | | | | | (_| | | | | | | |  __/
#|_|   \__,_|_|  |___/_|_| |_|\__, |  \___|_| |_| |_|\__,_| |_|_|_| |_|\___|
#                             |___/
args  = parser.MAELAS_Options()()


print("MAELAS version 3.0")
print(" ")
print("Authors: P. Nieves, S. Arapan, S.H. Zhang, A.P. Kądzielawa, R.F. Zhang and D. Legut ")
print(" ")


if not args.der and not args.gen and not args.rel and not args.mae:
    print("Please include tag -r or -g or -d or -m")
    exit(-1)

if (args.der and args.gen) or (args.gen and args.rel) or (args.der and args.rel):
    print("Please include tag -r or -g or -d. Only one of these tags.")
    exit(-1)

if (args.delas and not args.der):
    print("Tag -d should be included if you use tag -b")
    exit(-1)


if args.gen:

    if int(args.mode[0])==1:
      print('--------------------------------------------------------------------------------------------------------')
      print("Generation of VASP files for the calculation of anisotropic magnetostriction coefficients (mode=1):")
      print('--------------------------------------------------------------------------------------------------------')
    elif int(args.mode[0])==2:
      print('--------------------------------------------------------------------------------------------------------')
      print("Generation of VASP files for the calculation of anisotropic magnetoelastic constants (mode=2):")
      print('--------------------------------------------------------------------------------------------------------')
    elif int(args.mode[0])==3:
      print('--------------------------------------------------------------------------------------------------------')
      print("Generation of VASP files for the calculation of isotropic magnetoelastic constants (mode=3):")
      print('--------------------------------------------------------------------------------------------------------')
      
    generator = generate.VASP(args)
    symData   = generator.poscar()
    sg         = symData.space_group
    pg         = symData.point_group
    generator.incar()
    generator.kpoints()
    generator.scripts()


if args.rel:
    generator = relax.Relaxation(args)
    symData   = generator.poscar()
    generator.incar()
    generator.kpoints()
    generator.scripts()
    exit(0)

if args.mae:
    generator = test_mae.TestMAE(args)
    symData   = generator.poscar()
    generator.incar()
    generator.kpoints()
    generator.scripts()
    exit(0)

if args.der:
    
    if int(args.mode[0])==1:
        print('--------------------------------------------------------------------------------------------------------')
        print("Derivation of anisotropic magnetostrictive coefficients from the energy written in the OSZICAR files (mode=1):")
        print('--------------------------------------------------------------------------------------------------------')
    elif int(args.mode[0])==2:
        print('--------------------------------------------------------------------------------------------------------')
        print("Derivation of anisotropic magnetoelastic constants from the energy written in the OSZICAR files (mode=2):")
        print('--------------------------------------------------------------------------------------------------------')
    elif int(args.mode[0])==3:
        print('--------------------------------------------------------------------------------------------------------')
        print("Derivation of isotropic magnetoelastic constants from the energy written in the OSZICAR files (mode=3):")
        print('--------------------------------------------------------------------------------------------------------')

    structure0 = Structure.from_file(args.pos[0])
    sym1 = float(args.sympre[0])
    sym2 = float(args.symang[0])
    aa = SpacegroupAnalyzer(structure0,symprec=sym1, angle_tolerance=sym2)

    if int(args.sg0[0]) == 0:

        sg = aa.get_space_group_number()
        print("Space group number =", sg)

        spg = aa.get_space_group_symbol()
        print("Space group symbol =", str(spg))

        pg = aa.get_point_group_symbol()

    elif 0 < int(args.sg0[0]) <= 230:

        sg = int(args.sg0[0])
        print("Space group number (set by user)=", sg)

        spg = 'set by user'

        pg = 'set by user'

    else:
        print("Space group number must be in the range 1-230")
        exit


    if sg <= 15:
        print("Current version does not calculate magnetostriction for monoclinic and triclinic systems (space group < 16)")
        exit()
    elif 168 <= sg <= 176:
        print("Current version does not calculate magnetostriction for hexagonal (II) systems (167 < space group < 177)")
        exit()
    elif 143 <= sg <= 148:
        print("Current version does not calculate magnetostriction for trigonal (II) systems (142 < space group < 149)")
        exit()
    elif 75 <= sg <= 88:
        print("Current version does not calculate magnetostriction for tetragonal (II) systems (74 < space group < 89)")
        exit()
    elif 195 <= sg <= 206:
        print("Current version does not calculate magnetostriction for cubic (II) systems (193 < space group < 207)")
        exit()


    

#######################################################
#
##### CUBIC (I) ###########  230 >= space group >= 207 ############ MODE = 1
#
########################################################

if 230 >= sg >= 207 and int(args.mode[0])==1:
    
    print("Mode=",int(args.mode[0]), ", direct calculation of the magnetostrictive coefficients")
    print("Cubic (I) system")
    print("Point group =", str(pg))
    print("Number of anisotropic magnetostriction coefficients =", 2)
    
    cubic = cubic_I_mode1.run(args)
    
    if args.gen == True:
        cubic.gen()
    if args.der == True:
        cubic.der()



#######################################################
#
##### CUBIC (I) ###########  230 >= space group >= 207 ############ MODE = 2
#
########################################################

if 230 >= sg >= 207 and int(args.mode[0])==2:

    print("Mode=",int(args.mode[0]), ", direct calculation of the magnetoelastic constants (b)")
    print("Cubic (I) system")
    print("Point group =", str(pg))
    print("Number of anisotropic magnetoelastic constants =", 2)

    cubic = cubic_I_mode2.run(args)
    
    if args.gen == True:
        cubic.gen()
    if args.der == True:
        cubic.der()





########################################################################

### HEXAGONAL (I) and TETRAGONAL (I) ##### SG 177 - 194   &  SG 89 - 142  ###### MODE = 1

########################################################################


elif (177 <= sg <= 194 and int(args.mode[0])==1) or (89 <= sg <= 142 and int(args.mode[0])==1):


    if 177 <= sg <= 194:
        print("Mode=",int(args.mode[0]), ", direct calculation of the anisotropic magnetostrictive coefficients")
        print("Hexagonal (I) system")
        print("Point group =", str(pg))
        print("Number of anisotropic magnestostriction coefficients =", 4)



    if 89 <= sg <= 142:
        print("Mode=",int(args.mode[0]), ", direct calculation of the anisotropic magnetostrictive coefficients")
        print("Tetragonal (I) system")
        print("Point group =", str(pg))
        print("Number of anisotropic magnestostriction coefficients =", 5)
        
    hex_tet = hex_I_tet_I_mode1.run(args)
    
    if args.gen == True:
        hex_tet.gen()
    if args.der == True:
        hex_tet.der()



########################################################################

### HEXAGONAL (I) and TETRAGONAL (I) ##### SG 177 - 194   &  SG 89 - 142  ###### MODE = 2

########################################################################

elif (177 <= sg <= 194 and int(args.mode[0])==2) or (89 <= sg <= 142 and int(args.mode[0])==2):


    if 177 <= sg <= 194:
        print("Mode=",int(args.mode[0]), ", direct calculation of the anisotropic magnetoelastic constants")
        print("Hexagonal (I) system")
        print("Point group =", str(pg))
        print("Number of anisotropic magnetoelastic constants =", 4)



    if 89 <= sg <= 142:
        print("Mode=",int(args.mode[0]), ", direct calculation of the anisotropic magnetoelastic constants")
        print("Tetragonal (I) system")
        print("Point group =", str(pg))
        print("Number of anisotropic magnetoelastic constants =", 5)

    hex_tet = hex_I_tet_I_mode2.run(args)
    
    if args.gen == True:
        hex_tet.gen()
    if args.der == True:
        hex_tet.der()

        


#################################################################

##### TRIGONAL (I) ##### SG 149 - 167    ######## MODE = 1

#################################################################

elif 149 <= sg <= 167 and int(args.mode[0])==1:
    print("Mode=",int(args.mode[0]), ", direct calculation of the anisotropic magnetostrictive coefficients")
    print("Trigonal system")
    print("Point group =", str(pg))
    print("Number of anisotropic magnestostriction coefficients =", 6)
    
    tri = trig_I_mode1.run(args)
    
    if args.gen == True:
        tri.gen()
    if args.der == True:
        tri.der()

          
#################################################################

##### TRIGONAL (I) ##### SG 149 - 167    ######## MODE = 2

#################################################################



elif 149 <= sg <= 167 and int(args.mode[0])==2:
    print("Mode=",int(args.mode[0]), ", direct calculation of the anisotropic magnetoelastic constants")
    print("Trigonal system")
    print("Point group =", str(pg))
    print("Number of anisotropic magnetoelastic constants =", 6)
    
    
    tri = trig_I_mode2.run(args)
    
    if args.gen == True:
        tri.gen()
    if args.der == True:
        tri.der()
        

#################################################################

##### ORTHORHOMBIC #### SG 16 - 74    ####### MODE = 1

#################################################################


elif 16 <= sg <= 74 and int(args.mode[0])==1:
    print("Mode=",int(args.mode[0]), ", direct calculation of the anisotropic magnetostrictive coefficients")
    print("Orthorhombic system")
    print("Point group =", str(pg))
    print("Number of anisotropic magnestostriction coefficients =", 9)
    
    ort = ort_mode1.run(args)
    
    if args.gen == True:
        ort.gen()
    if args.der == True:
        ort.der()
    
    
#################################################################

##### ORTHORHOMBIC #### SG 16 - 74     ##### MODE = 2

#################################################################


elif 16 <= sg <= 74 and int(args.mode[0])==2:
    print("Mode=",int(args.mode[0]), ", direct calculation of the anisotropic magnetoelastic constants")
    print("Orthorhombic system")
    print("Point group =", str(pg))
    print("Number of anisotropic magnetoelastic constants =", 9)
    
    ort = ort_mode2.run(args)
    
    if args.gen == True:
        ort.gen()
    if args.der == True:
        ort.der()



#######################################################
#
#####  MODE = 3 for cubic I, hex I, trig I, tet I and ortho
#
########################################################


if int(args.mode[0])==3:
    print("Mode=",int(args.mode[0]), ", direct calculation of the isotropic magnetoelastic constants (b)")
    
    mode3 = mode3.run(args)
    
    if args.gen == True:
        mode3.gen()
    if args.der == True:
        mode3.der()
    

        