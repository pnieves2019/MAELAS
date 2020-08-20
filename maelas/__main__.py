#!/bin/bash

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import stat

import maelas.parser   as parser
import maelas.generate as generate
import maelas.relax    as relax
import maelas.test_mae as test_mae
from   maelas.data import SymmetryData

from pymatgen import Lattice, Structure
from pymatgen.transformations.standard_transformations import ConventionalCellTransformation,DeformStructureTransformation
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar
from scipy.optimize import curve_fit
from pyfiglet import Figlet
from sklearn.metrics import r2_score


##############################################

f = Figlet(font='slant')
print(f.renderText('MAELAS'))

# ____                _                                  _   _ _
#|  _ \ __ _ _ __ ___(_)_ __   __ _    ___ _ __ ___   __| | | (_)_ __   ___
#| |_) / _` | '__/ __| | '_ \ / _` |  / __| '_ ` _ \ / _` | | | | '_ \ / _ \
#|  __/ (_| | |  \__ \ | | | | (_| | | (__| | | | | | (_| | | | | | | |  __/
#|_|   \__,_|_|  |___/_|_| |_|\__, |  \___|_| |_| |_|\__,_| |_|_|_| |_|\___|
#                             |___/
args  = parser.MAELAS_Options()()


print("MAELAS code v1.0")
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
    generator = generate.VASP(args)
    symData   = generator.poscar()
    nat        = symData.number_of_species
    structure2 = symData.structure
    sg         = symData.space_group
    pg         = symData.point_group


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

if args.der == True:
    print('--------------------------------------------------------------------------------------------------------')
    print("Derivation of anisotropic magnetostriction coefficients from the energy written in the OSZICAR files:")
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

    elif int(args.sg0[0]) > 0:

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

if args.gen == True:
    generator.incar()
    generator.kpoints()
    generator.scripts()

#######################################################
#
##### CUBIC (I) ###########  230 >= space group >= 207
#
########################################################


if 230 >= sg >= 207:
    print("Cubic (I) system")
    print("Point group =", str(pg))
    print("Number of anisotropic magnestostriction coefficients =", 2)

    if args.gen == True:

        for i in range(int(args.ndist[0])):


            strain1 = - float(args.strain[0])+2*(float(args.strain[0])/(float(args.ndist[0])-1))*i

            print("epsilon", strain1)


        #Generation POSCAR file

        #lambda_001


            a3 = 1.0 + strain1
            a1 = 1/math.sqrt(a3)
            a2 = a1
            dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
            structure3 = dd.apply_transformation(structure2)
            pos_name = "POSCAR_1_" + str(i+1)

            structure33 = Poscar(structure3)
            structure33.write_file(filename = pos_name,significant_figures=16)



        #lambda_111

            const = (4/(4-3*(strain1**2)+strain1**3))**(1/3)

            a12 = const*strain1*0.5
            a13 = a12
            a21 = a12
            a22 = a12
            a23 = a12
            a31 = a12
            a32 = a12
            a33 = a12

            a11 = const
            a22 = const
            a33 = const

            cc = DeformStructureTransformation(deformation=((a11, a12, a13), (a21, a22, a23), (a31, a32, a33)))
            structure4 = cc.apply_transformation(structure2)
            pos_name2 = "POSCAR_2_" + str(i+1)

            structure44 = Poscar(structure4)
            structure44.write_file(filename = pos_name2,significant_figures=16)



    # INCAR_1_1 m=0,0,1

        path_inc_ncl_1_1 = 'INCAR_1_1'
        inc_ncl_1_1 = open(path_inc_ncl_1_1,'w')
        inc_ncl_list_1_1 = generator.inc_ncl_list[:]
        inc_ncl_list_1_1 += ['SAXIS = 0 0 1.0\n']

        for j in range(len(inc_ncl_list_1_1)):
            inc_ncl_1_1.write(str(inc_ncl_list_1_1[j]))

        inc_ncl_1_1.close()


    # INCAR_1_2 m=1,0,0

        path_inc_ncl_1_2 = 'INCAR_1_2'
        inc_ncl_1_2 = open(path_inc_ncl_1_2,'w')
        inc_ncl_list_1_2 = generator.inc_ncl_list[:]
        inc_ncl_list_1_2 += ['SAXIS = 1.0 0 0.0\n']

        for j in range(len(inc_ncl_list_1_2)):
            inc_ncl_1_2.write(str(inc_ncl_list_1_2[j]))

        inc_ncl_1_2.close()

    # INCAR_2_1 m=1,1,1

        path_inc_ncl_2_1 = 'INCAR_2_1'
        inc_ncl_2_1 = open(path_inc_ncl_2_1,'w')
        inc_ncl_list_2_1 = generator.inc_ncl_list[:]
        inc_ncl_list_2_1 += ['SAXIS = 1.0 1.0 1.0\n']

        for j in range(len(inc_ncl_list_2_1)):
            inc_ncl_2_1.write(str(inc_ncl_list_2_1[j]))

        inc_ncl_2_1.close()


    # INCAR_2_2 m=1,0,-1

        path_inc_ncl_2_2 = 'INCAR_2_2'
        inc_ncl_2_2 = open(path_inc_ncl_2_2,'w')
        inc_ncl_list_2_2 = generator.inc_ncl_list[:]
        inc_ncl_list_2_2 += ['SAXIS = 1.0 0.0 -1.0\n']

        for j in range(len(inc_ncl_list_2_2)):
            inc_ncl_2_2.write(str(inc_ncl_list_2_2[j]))

        inc_ncl_2_2.close()


    # Derivation of magnetostriction coefficients:

    if args.der == True:

        for j in range(1,3):

            for k in range(1,3):

                path_dat = "ene_" + str(j) + "_" + str(k) + ".dat"
                dat = open(path_dat,'w')


                for i in range(int(args.ndist[0])):

                    pos_name = "POSCAR_" + str(j) + "_" + str(i+1)

                    struct = Structure.from_file(pos_name)

                    latt = struct.lattice.matrix

                    if j == 1:
                        var1 = latt[2][2]
                    elif j == 2:
                        var1 = math.sqrt((latt[0][0]+latt[1][0]+latt[2][0])**2+(latt[0][1]+latt[1][1]+latt[2][1])**2+(latt[0][2]+latt[1][2]+latt[2][2])**2)


                    path_osz = "OSZICAR_" + str(j) + "_" + str(i+1) + "_" + str(k)
                    osz = open(path_osz,'r')
                    ene0 = osz.readlines()
                    ene1 = ene0[len(ene0)-2]
                    ene2 = ene1[11:32]

                    osz.close()

                    dat.write(repr(var1))
                    dat.write('  ')
                    dat.write(str(ene2))
                    dat.write('\n')


                dat.close()


       # fitting and plot


        def K(x,a,b,c):
            return a*x*x+b*x+c


        print("")
        print("Fit of quadratic function f(x)=a*x\u00B2+b*x+c to energy vs distortion data")
        print("")
        print("-------------------------")
        print("Calculation of \u03BB001:")
        print("-------------------------")
        print(" ")
        print('Lattice distorsion along [0,0,1] direction')
        print("")

        f = open('ene_1_1.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []
        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)

        print("Fitting parameters for spin parallel to 001 (data from file ene_1_1.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_1_1.png")
            print("")

        l1 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l1)
        print("")


        plt.plot(x, y, 'bo', label='data in ene_1_1.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Lattice distorsion along [0,0,1] direction (Å)')
        plt.title('Calculation of \u03BB\u2080\u2080\u2081 (spin = [0,0,1]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)

        plt.savefig('fit_ene_1_1.png')
        plt.close()

        f = open('ene_1_2.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []
        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)
        print("Fitting parameters for spin parallel to 100 (data from file ene_1_2.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_1_2.png")
            print("")

        l2 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l2)
        print("")


        lambda001 = 2.0*(2.0/3.0)*((l1 -l2)/(l1+l2))


        plt.plot(x, y, 'bo', label='data in ene_1_2.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Lattice distorsion along [0,0,1] direction (Å)')
        plt.title('Calculation of \u03BB\u2080\u2080\u2081 (spin = [1,0,0]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)

        plt.savefig('fit_ene_1_2.png')
        plt.close()



        #make figure dE_1.png

        fig = 'dE_1.png'
        spin1 = '0,0,1'
        spin2 = '1,0,0'
        dist = '0,0,1'
        tit = "Calculation of \u03BB\u2080\u2080\u2081"
        f1 = open('ene_1_1.dat','r')
        f2 = open('ene_1_2.dat','r')


        s1 = f1.readlines()
        s2 = f2.readlines()
        f1.close
        f2.close

        x = []
        y = []
        y2 = []
        for j in s1:
            x.append(float(j.split()[0]))
            y.append(float(j.split()[1]))

        for j in s2:
                y2.append(float(j.split()[1]))

        x = np.array(x)
        y = np.array(y)
        y2 = np.array(y2)

        plt.plot(x, (y2-y)*1e6, 'o-')

        ylabel ='E[' + str(spin2) + '] - E['+ str(spin1) + '] (\u03BCeV)'
        plt.ylabel(ylabel)
        label = "Lattice distorsion along [" + str(dist) + "] direction (Å)"
        plt.xlabel(label)
        plt.title(tit)
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig(fig)
        plt.close()







        print(" ")
        print("-------------------------")
        print("Calculation of \u03BB111:")
        print("-------------------------")
        print(" ")
        print('Lattice distorsion along [1,1,1] direction')
        print("")

        f = open('ene_2_1.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []
        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)
        print("Fitting parameters for spin parallel to 111 (data from file ene_2_1.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_2_1.png")
            print("")

        l1 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l1)
        print("")


        plt.plot(x, y, 'bo', label='data in ene_2_1.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Lattice distorsion along [1,1,1] direction (Å)')
        plt.title('Calculation of \u03BB\u2081\u2081\u2081 (spin = [1,1,1]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)

        plt.savefig('fit_ene_2_1.png')
        plt.close()


        f = open('ene_2_2.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []

        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)
        print("Fitting parameters for spin parallel to 10-1 (data from file ene_2_2.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_2_2.png")
            print("")

        l2 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l2)
        print("")

        lambda111 = 2.0*(2.0/3.0)*((l1 -l2)/(l1+l2))

        lambda_s = (2.0/5.0)*lambda001 + (3.0/5.0)*lambda111


        plt.plot(x, y, 'bo', label='data in ene_2_2.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.xlabel('Lattice distorsion along [1,1,1] direction (Å)')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.title('Calculation of \u03BB\u2081\u2081\u2081 (spin = [1,0,-1]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig('fit_ene_2_2.png')
        plt.close()



        #make figure dE_2.png

        fig = 'dE_2.png'
        spin1 = '1,1,1'
        spin2 = '1,0,-1'
        dist = '1,1,1'
        tit = "Calculation of \u03BB\u2081\u2081\u2081"
        f1 = open('ene_2_1.dat','r')
        f2 = open('ene_2_2.dat','r')


        s1 = f1.readlines()
        s2 = f2.readlines()
        f1.close
        f2.close

        x = []
        y = []
        y2 = []
        for j in s1:
            x.append(float(j.split()[0]))
            y.append(float(j.split()[1]))

        for j in s2:
                y2.append(float(j.split()[1]))

        x = np.array(x)
        y = np.array(y)
        y2 = np.array(y2)

        plt.plot(x, (y2-y)*1e6, 'o-')

        ylabel ='E[' + str(spin2) + '] - E['+ str(spin1) + '] (\u03BCeV)'
        plt.ylabel(ylabel)
        label = "Lattice distorsion along [" + str(dist) + "] direction (Å)"
        plt.xlabel(label)
        plt.title(tit)
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig(fig)
        plt.close()








        print(" ")
        print("----------------------------------------------")
        print("Anisotropic magnetostriction coefficients:")
        print("----------------------------------------------")
        print(" ")
        print(" ")
        print("Using the convention in reference J.R. Cullen et al., in Materials, Science and Technology (VCH Publishings, 1994), pp.529-565:")
        print(" ")
        print("\u03BB001 =", lambda001*1e6,u'x 10\u207B\u2076')
        print(" ")
        print("\u03BB111 =", lambda111*1e6,u'x 10\u207B\u2076')
        print(" ")
        print("(Polycrystal) \u03BBs =", lambda_s*1e6,u'x 10\u207B\u2076')


        if args.delas == True:
            print(" ")
            print(" ")
            print("----------------------------------------------")
            print("Calculation of magnetoelastic constants:")
            print("----------------------------------------------")
            print(" ")
            print("Reading the elastic tensor file =", str(args.elas[0]))
            print(" ")




            elasdat = open(args.elas[0],'r')
            elasline = elasdat.readlines()
            elasline0 = elasline[2]
            elasline1 = elasline[5]
            c11 = float(elasline0[0:8])
            c12 = float(elasline0[8:16])
            c44 = float(elasline1[24:32])

            elasdat.close()


            b1 = -(3/2)*(c11-c12)*lambda001

            b2 = -3*c44*lambda111


            print("c11 =", str(c11), 'GPa')
            print(" ")
            print("c12 =", str(c12), 'GPa')
            print(" ")
            print("c44 =", str(c44), 'GPa')
            print(" ")
            print("Warning: If these elastic constants are not the same as in the input elastic tensor file", str(args.elas[0]),", then check that the format of the elastic tensor is exactly the same as in the standard output file ELADAT generated by ELAS code (see Example folder)")
            print(" ")
            print(" ")
            print("Magnetoelastic constants:")
            print(" ")
            print("b1 =", str(b1), 'GPa')
            print(" ")
            print("b2 =", str(b2), 'GPa')
            print(" ")









########################################################################

### HEXAGONAL (I) and TETRAGONAL (I) ##### SG 177 - 194   &  SG 89 - 142

########################################################################


elif (177 <= sg <= 194) or (89 <= sg <= 142):


    if 177 <= sg <= 194:
        print("Hexagonal (I) system")
        print("Point group =", str(pg))
        print("Number of anisotropic magnestostriction coefficients =", 4)



    if 89 <= sg <= 142:
        print("Tetragonal (I) system")
        print("Point group =", str(pg))
        print("Number of anisotropic magnestostriction coefficients =", 5)

    if args.gen == True:

        if 177 <= sg <= 194:
            # Convention: lattice vector a1 along x-axis
            angle = -math.pi*(60.0/180.0)
            dd = DeformStructureTransformation(deformation=((math.cos(angle), math.sin(angle), 0), (-math.sin(angle), math.cos(angle), 0), (0, 0, 1)))
            structure2b = dd.apply_transformation(structure2)
        else:
            structure2b = structure2



        for i in range(int(args.ndist[0])):


            strain1 = - float(args.strain[0])+2*(float(args.strain[0])/(float(args.ndist[0])-1))*i

            print("epsilon", strain1)


        #Generation POSCAR file

        #lambda_alpha_1_2




            a1 = 1.0 + strain1
            a2 = 1/math.sqrt(a1)
            a3 = a2
            dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
            structure3 = dd.apply_transformation(structure2b)
            pos_name = "POSCAR_1_" + str(i+1)

            structure33 = Poscar(structure3)
            structure33.write_file(filename = pos_name,significant_figures=16)

        #lambda_alpha_2_2


            a3 = 1.0 + strain1
            a1 = 1/math.sqrt(a3)
            a2 = a1
            dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
            structure3 = dd.apply_transformation(structure2b)
            pos_name2 = "POSCAR_2_" + str(i+1)

            structure33 = Poscar(structure3)
            structure33.write_file(filename = pos_name2,significant_figures=16)

        #lambda_gamma_2


            a1 = 1.0 + strain1
            a2 = 1/math.sqrt(a1)
            a3 = a2
            dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
            structure3 = dd.apply_transformation(structure2b)
            pos_name3 = "POSCAR_3_" + str(i+1)

            structure33 = Poscar(structure3)
            structure33.write_file(filename = pos_name3,significant_figures=16)

        #lambda_epsilon_2

            latt_par = structure2b.lattice.matrix
            
            const = (1/(1-(strain1*0.5)**2))**(1/3)
            
            a11 = const
            a12 = 0.0
            a13 = const*strain1*0.5*(latt_par[2][2]/latt_par[0][0])
            a21 = 0.0
            a22 = const
            a23 = 0.0
            a31 = const*strain1*0.5*(latt_par[0][0]/latt_par[2][2])
            a32 = 0.0
            a33 = const

            cc = DeformStructureTransformation(deformation=((a11, a12, a13), (a21, a22, a23), (a31, a32, a33)))
            structure4 = cc.apply_transformation(structure2b)
            pos_name4 = "POSCAR_4_" + str(i+1)

            structure44 = Poscar(structure4)
            structure44.write_file(filename = pos_name4,significant_figures=16)



            if 89 <= sg <= 142:
                # lambda delta2
                const = (1/(1-(strain1*0.5)**2))**(1/3)

                a11 = const
                a12 = const*strain1*0.5
                a13 = 0.0
                a21 = const*strain1*0.5
                a22 = const
                a23 = 0.0
                a31 = 0.0
                a32 = 0.0
                a33 = const

                cc = DeformStructureTransformation(deformation=((a11, a12, a13), (a21, a22, a23), (a31, a32, a33)))
                structure5 = cc.apply_transformation(structure2b)
                pos_name5 = "POSCAR_5_" + str(i+1)

                structure55 = Poscar(structure5)
                structure55.write_file(filename = pos_name5,significant_figures=16)

        # INCAR_1_1 m=1,1,1

        path_inc_ncl_1_1 = 'INCAR_1_1'
        inc_ncl_1_1 = open(path_inc_ncl_1_1,'w')
        inc_ncl_list_1_1 = generator.inc_ncl_list[:]
        inc_ncl_list_1_1 += ['SAXIS = 1.0 1.0 1.0\n']

        for j in range(len(inc_ncl_list_1_1)):
            inc_ncl_1_1.write(str(inc_ncl_list_1_1[j]))

        inc_ncl_1_1.close()


    # INCAR_1_2 m=1,1,0

        path_inc_ncl_1_2 = 'INCAR_1_2'
        inc_ncl_1_2 = open(path_inc_ncl_1_2,'w')
        inc_ncl_list_1_2 = generator.inc_ncl_list[:]
        inc_ncl_list_1_2 += ['SAXIS = 1.0 1.0 0.0\n']

        for j in range(len(inc_ncl_list_1_2)):
            inc_ncl_1_2.write(str(inc_ncl_list_1_2[j]))

        inc_ncl_1_2.close()



    # INCAR_2_1 m=0,0,1

        path_inc_ncl_2_1 = 'INCAR_2_1'
        inc_ncl_2_1 = open(path_inc_ncl_2_1,'w')
        inc_ncl_list_2_1 = generator.inc_ncl_list[:]
        inc_ncl_list_2_1 += ['SAXIS = 0.0 0.0 1.0\n']

        for j in range(len(inc_ncl_list_2_1)):
            inc_ncl_2_1.write(str(inc_ncl_list_2_1[j]))

        inc_ncl_2_1.close()


    # INCAR_2_2 m=1,0,0

        path_inc_ncl_2_2 = 'INCAR_2_2'
        inc_ncl_2_2 = open(path_inc_ncl_2_2,'w')
        inc_ncl_list_2_2 = generator.inc_ncl_list[:]
        inc_ncl_list_2_2 += ['SAXIS = 1.0 0.0 0.0\n']

        for j in range(len(inc_ncl_list_2_2)):
            inc_ncl_2_2.write(str(inc_ncl_list_2_2[j]))

        inc_ncl_2_2.close()


        # INCAR_3_1 m=1,0,0

        path_inc_ncl_3_1 = 'INCAR_3_1'
        inc_ncl_3_1 = open(path_inc_ncl_3_1,'w')
        inc_ncl_list_3_1 = generator.inc_ncl_list[:]
        inc_ncl_list_3_1 += ['SAXIS = 1.0 0.0 0.0\n']

        for j in range(len(inc_ncl_list_3_1)):
            inc_ncl_3_1.write(str(inc_ncl_list_3_1[j]))

        inc_ncl_3_1.close()


    # INCAR_3_2 m=0,1,0

        path_inc_ncl_3_2 = 'INCAR_3_2'
        inc_ncl_3_2 = open(path_inc_ncl_3_2,'w')
        inc_ncl_list_3_2 = generator.inc_ncl_list[:]
        inc_ncl_list_3_2 += ['SAXIS = 0.0 1.0 0.0\n']

        for j in range(len(inc_ncl_list_3_2)):
            inc_ncl_3_2.write(str(inc_ncl_list_3_2[j]))

        inc_ncl_3_2.close()


        # INCAR_4_1 m=1,0,1

        path_inc_ncl_4_1 = 'INCAR_4_1'
        inc_ncl_4_1 = open(path_inc_ncl_4_1,'w')
        inc_ncl_list_4_1 = generator.inc_ncl_list[:]
        inc_ncl_list_4_1 += ['SAXIS = 1.0 0.0 1.0\n']

        for j in range(len(inc_ncl_list_4_1)):
            inc_ncl_4_1.write(str(inc_ncl_list_4_1[j]))

        inc_ncl_4_1.close()


        # INCAR_4_2 m=-1,0,1

        path_inc_ncl_4_2 = 'INCAR_4_2'
        inc_ncl_4_2 = open(path_inc_ncl_4_2,'w')
        inc_ncl_list_4_2 = generator.inc_ncl_list[:]
        inc_ncl_list_4_2 += ['SAXIS = -1.0 0.0 1.0\n']

        for j in range(len(inc_ncl_list_4_2)):
            inc_ncl_4_2.write(str(inc_ncl_list_4_2[j]))

        inc_ncl_4_2.close()



        if 89 <= sg <= 142:

            # INCAR_5_1 m=1,1,0

            path_inc_ncl_5_1 = 'INCAR_5_1'
            inc_ncl_5_1 = open(path_inc_ncl_5_1,'w')
            inc_ncl_list_5_1 = generator.inc_ncl_list[:]
            inc_ncl_list_5_1 += ['SAXIS = 1.0 1.0 0.0\n']

            for j in range(len(inc_ncl_list_5_1)):
                inc_ncl_5_1.write(str(inc_ncl_list_5_1[j]))

            inc_ncl_5_1.close()


            # INCAR_5_2 m=-1,1,0

            path_inc_ncl_5_2 = 'INCAR_5_2'
            inc_ncl_5_2 = open(path_inc_ncl_5_2,'w')
            inc_ncl_list_5_2 = generator.inc_ncl_list[:]
            inc_ncl_list_5_2 += ['SAXIS = -1.0 1.0 0.0\n']

            for j in range(len(inc_ncl_list_5_2)):
                inc_ncl_5_2.write(str(inc_ncl_list_5_2[j]))

            inc_ncl_5_2.close()




  # Derivation of magnetostriction coefficients:

    if args.der == True:

        if (89 <= sg <= 142):
            nmax = 6
        else:
            nmax = 5


        for j in range(1,nmax):

            for k in range(1,3):

                path_dat = "ene_" + str(j) + "_" + str(k) + ".dat"
                dat = open(path_dat,'w')


                for i in range(int(args.ndist[0])):

                    pos_name = "POSCAR_" + str(j) + "_" + str(i+1)

                    struct = Structure.from_file(pos_name)

                    latt = struct.lattice.matrix

                    if j == 1:
                        var1 = latt[0][0]

                    elif j == 2:
                        var1 = latt[2][2]

                    elif j == 3:
                        var1 = latt[0][0]

                    elif j == 4:
                        var1 = math.sqrt((latt[0][0]+latt[2][0])**2+(latt[0][1]+latt[2][1])**2+(latt[0][2]+latt[2][2])**2)

                    else:
                        if 89 <= sg <= 142:
                            var1 = math.sqrt((latt[0][0]+latt[1][0])**2+(latt[0][1]+latt[1][1])**2+(latt[0][2]+latt[1][2])**2)



                    path_osz = "OSZICAR_" + str(j) + "_" + str(i+1) + "_" + str(k)
                    osz = open(path_osz,'r')
                    ene0 = osz.readlines()
                    ene1 = ene0[len(ene0)-2]
                    ene2 = ene1[11:32]

                    osz.close()

                    dat.write(repr(var1))
                    dat.write('  ')
                    dat.write(str(ene2))
                    dat.write('\n')


                dat.close()



        # fitting and plot


        def K(x,a,b,c):
            return a*x*x+b*x+c


        print("")
        print("Fit of quadratic function f(x)=a*x\u00B2+b*x+c to energy vs distortion data")

        print(" ")
        print("-------------------------")
        print('Calculation of \u03BB 1\u03B1,2:')
        print("-------------------------")
        print(" ")
        print('Lattice distorsion along [1,0,0] direction')
        print("")

        f = open('ene_1_1.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []
        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)


        print("Fitting parameters for spin parallel to 111 (data from file ene_1_1.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_1_1.png")
            print("")

        l1 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l1)
        print("")


        plt.plot(x, y, 'bo', label='data in ene_1_1.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Lattice distorsion along [1,0,0] direction (Å)')
        plt.title('Calculation of \u03BB 1\u03B1,2 (spin = [1,1,1]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)

        plt.savefig('fit_ene_1_1.png')
        plt.close()

        f = open('ene_1_2.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []
        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)
        print("Fitting parameters for spin parallel to 110 (data from file ene_1_2.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])
        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_1_2.png")
            print("")

        l2 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l2)
        print("")

        lambda_alpha_1_2 = 2.0*3.0*((l1 -l2)/(l1+l2))


        plt.plot(x, y, 'bo', label='data in ene_1_2.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Lattice distorsion along [1,0,0] direction (Å)')
        plt.title('Calculation of \u03BB 1\u03B1,2 (spin = [1,1,0]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)

        plt.savefig('fit_ene_1_2.png')
        plt.close()



        #make figure dE_1.png

        fig = 'dE_1.png'
        spin1 = '1,1,1'
        spin2 = '1,1,0'
        dist = '1,0,0'
        tit = "Calculation of \u03BB 1\u03B1,2  "
        f1 = open('ene_1_1.dat','r')
        f2 = open('ene_1_2.dat','r')


        s1 = f1.readlines()
        s2 = f2.readlines()
        f1.close
        f2.close

        x = []
        y = []
        y2 = []
        for j in s1:
            x.append(float(j.split()[0]))
            y.append(float(j.split()[1]))

        for j in s2:
                y2.append(float(j.split()[1]))

        x = np.array(x)
        y = np.array(y)
        y2 = np.array(y2)

        plt.plot(x, (y2-y)*1e6, 'o-')

        ylabel ='E[' + str(spin2) + '] - E['+ str(spin1) + '] (\u03BCeV)'
        plt.ylabel(ylabel)
        label = "Lattice distorsion along [" + str(dist) + "] direction (Å)"
        plt.xlabel(label)
        plt.title(tit)
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig(fig)
        plt.close()






        print(" ")
        print("-------------------------")
        print("Calculation of \u03BB 2\u03B1,2:")
        print("-------------------------")
        print(" ")
        print('Lattice distorsion along [0,0,1] direction')
        print("")


        f = open('ene_2_1.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []
        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)

        print("Fitting parameters for spin parallel to 001 (data from file ene_2_1.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])
        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_2_1.png")
            print("")

        l1 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l1)
        print("")


        plt.plot(x, y, 'bo', label='data in ene_2_1.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Lattice distorsion along [0,0,1] direction (Å)')
        plt.title('Calculation of \u03BB 2\u03B1,2 (spin = [0,0,1]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)

        plt.savefig('fit_ene_2_1.png')
        plt.close()


        f = open('ene_2_2.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []

        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)
        print("Fitting parameters for spin parallel to 100 (data from file ene_2_2.dat)::")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])
        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_2_2.png")
            print("")

        l2 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l2)
        print("")

        lambda_alpha_2_2 = 2.0*((l1 -l2)/(l1+l2))


        plt.plot(x, y, 'bo', label='data in ene_2_2.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.xlabel('Lattice distorsion along [0,0,1] direction (Å)')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.title('Calculation of \u03BB 2\u03B1,2 (spin = [1,0,0]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig('fit_ene_2_2.png')
        plt.close()


       #make figure dE_2.png

        fig = 'dE_2.png'
        spin1 = '0,0,1'
        spin2 = '1,0,0'
        dist = '0,0,1'
        tit = "Calculation of \u03BB 2\u03B1,2  "
        f1 = open('ene_2_1.dat','r')
        f2 = open('ene_2_2.dat','r')


        s1 = f1.readlines()
        s2 = f2.readlines()
        f1.close
        f2.close

        x = []
        y = []
        y2 = []
        for j in s1:
            x.append(float(j.split()[0]))
            y.append(float(j.split()[1]))

        for j in s2:
                y2.append(float(j.split()[1]))

        x = np.array(x)
        y = np.array(y)
        y2 = np.array(y2)

        plt.plot(x, (y2-y)*1e6, 'o-')

        ylabel ='E[' + str(spin2) + '] - E['+ str(spin1) + '] (\u03BCeV)'
        plt.ylabel(ylabel)
        label = "Lattice distorsion along [" + str(dist) + "] direction (Å)"
        plt.xlabel(label)
        plt.title(tit)
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig(fig)
        plt.close()







        print(" ")
        print("-------------------------")
        print("Calculation of \u03BB \u03B3,2:")
        print("-------------------------")
        print(" ")
        print('Lattice distorsion along [1,0,0] direction')
        print("")


        f = open('ene_3_1.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []
        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)

        print("Fitting parameters for spin parallel to 100 (data from file ene_3_1.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])
        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_3_1.png")
            print("")

        l1 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l1)
        print("")


        plt.plot(x, y, 'bo', label='data in ene_3_1.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Lattice distorsion along [1,0,0] direction (Å)')
        plt.title('Calculation of \u03BB \u03B3,2 (spin = [1,0,0]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)

        plt.savefig('fit_ene_3_1.png')
        plt.close()


        f = open('ene_3_2.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []

        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)
        print("Fitting parameters for spin parallel to 010 (data from file ene_3_2.dat:")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])
        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_3_2.png")
            print("")

        l2 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l2)
        print("")

        lambda_gamma_2 = 2.0*((l1 -l2)/(l1+l2))


        plt.plot(x, y, 'bo', label='data in ene_3_2.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.xlabel('Lattice distorsion along [1,0,0] direction (Å)')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.title('Calculation of \u03BB \u03B3,2 (spin = [0,1,0]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig('fit_ene_3_2.png')
        plt.close()



        #make figure dE_3.png

        fig = 'dE_3.png'
        spin1 = '1,0,0'
        spin2 = '0,1,0'
        dist = '1,0,0'
        tit = "Calculation of \u03BB \u03B3,2  "
        f1 = open('ene_3_1.dat','r')
        f2 = open('ene_3_2.dat','r')


        s1 = f1.readlines()
        s2 = f2.readlines()
        f1.close
        f2.close

        x = []
        y = []
        y2 = []
        for j in s1:
            x.append(float(j.split()[0]))
            y.append(float(j.split()[1]))

        for j in s2:
                y2.append(float(j.split()[1]))

        x = np.array(x)
        y = np.array(y)
        y2 = np.array(y2)

        plt.plot(x, (y2-y)*1e6, 'o-')

        ylabel ='E[' + str(spin2) + '] - E['+ str(spin1) + '] (\u03BCeV)'
        plt.ylabel(ylabel)
        label = "Lattice distorsion along [" + str(dist) + "] direction (Å)"
        plt.xlabel(label)
        plt.title(tit)
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig(fig)
        plt.close()








        print(" ")
        print("-------------------------")
        print("Calculation of \u03BB \u03B5,2:")
        print("-------------------------")
        print(" ")
        print('Lattice distorsion along [1,0,1] direction')
        print("")


        f = open('ene_4_1.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []
        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)
        print("Fitting parameters for spin parallel to 101 (data from file ene_4_1.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])
        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_4_1.png")
            print("")

        l1 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l1)
        print("")


        plt.plot(x, y, 'bo', label='data in ene_4_1.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Lattice distorsion along [1,0,1] direction (Å)')
        plt.title('Calculation of \u03BB \u03B5,2 (spin = [1,0,1]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)

        plt.savefig('fit_ene_4_1.png')
        plt.close()


        f = open('ene_4_2.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []

        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)
        print("Fitting parameters for spin parallel to -101 (data from file ene_4_2.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])
        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_4_2.png")
            print("")

        l2 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l2)
        print("")

        
        aa0 = SpacegroupAnalyzer(structure0,symprec=sym1, angle_tolerance=sym2)
        structure1 = aa0.get_conventional_standard_structure(international_monoclinic=True)

        bb0 = ConventionalCellTransformation(symprec=sym1, angle_tolerance=sym2, international_monoclinic=True)
        structure2 = bb0.apply_transformation(structure1) 
        
        if 177 <= sg <= 194:
            # Convention: lattice vector a1 along x-axis
            angle = -math.pi*(60.0/180.0)
            dd = DeformStructureTransformation(deformation=((math.cos(angle), math.sin(angle), 0), (-math.sin(angle), math.cos(angle), 0), (0, 0, 1)))
            structure2b = dd.apply_transformation(structure2)
        else:
            structure2b = structure2
        
        
        latt_par = structure2b.lattice.matrix
            
        latt_a = latt_par[0][0]
        latt_c = latt_par[2][2]
        
        
        eta_par = (2.0*latt_a*latt_c)/(latt_a**2+latt_c**2)
        
        lambda_epsilon_2 = 2.0*((l1-l2)/(l1+l2))*(1.0/eta_par)




        plt.plot(x, y, 'bo', label='data in ene_4_2.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.xlabel('Lattice distorsion along [1,0,1] direction (Å)')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.title('Calculation of \u03BB \u03B5,2 (spin = [-1,0,1]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig('fit_ene_4_2.png')
        plt.close()



        #make figure dE_4.png

        fig = 'dE_4.png'
        spin1 = '1,0,1'
        spin2 = '-1,0,1'
        dist = '1,0,1'
        tit = "Calculation of \u03BB \u03B5,2  "
        f1 = open('ene_4_1.dat','r')
        f2 = open('ene_4_2.dat','r')


        s1 = f1.readlines()
        s2 = f2.readlines()
        f1.close
        f2.close

        x = []
        y = []
        y2 = []
        for j in s1:
            x.append(float(j.split()[0]))
            y.append(float(j.split()[1]))

        for j in s2:
                y2.append(float(j.split()[1]))

        x = np.array(x)
        y = np.array(y)
        y2 = np.array(y2)

        plt.plot(x, (y2-y)*1e6, 'o-')

        ylabel ='E[' + str(spin2) + '] - E['+ str(spin1) + '] (\u03BCeV)'
        plt.ylabel(ylabel)
        label = "Lattice distorsion along [" + str(dist) + "] direction (Å)"
        plt.xlabel(label)
        plt.title(tit)
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig(fig)
        plt.close()







        if 89 <= sg <= 142:

            print(" ")
            print("-------------------------")
            print("Calculation of \u03BB \u03B4,2:")
            print("-------------------------")
            print(" ")
            print('Lattice distorsion along [1,1,0] direction')
            print("")

            f = open('ene_5_1.dat','r')
            l = f.readlines()
            f.close

            x = []
            y = []
            for i in l:
                x.append(float(i.split()[0]))
                y.append(float(i.split()[1]))

            x = np.array(x)
            y = np.array(y)

            params = curve_fit(K, x, y)
            print("Fitting parameters for spin parallel to 110 (data from file ene_5_1.dat):")
            print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])
            r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
            print("R-squared =", r_squared)
            print("")

            if r_squared < 0.98:
                print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_5_1.png")
                print("")

            l1 = -params[0][1] / (2.0 * params[0][0])

            print("X minimum = -b/(2*a) =", l1)
            print("")

            plt.plot(x, y, 'bo', label='data in ene_5_1.dat')
            popt, pcov = curve_fit(K, x, y)
            t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
            plt.plot(t, K(t, *popt), 'r--', label='fit')
            plt.ylabel('Energy (eV)')
            plt.legend()
            plt.xlabel('Lattice distorsion along [1,1,0] direction (Å)')
            plt.title('Calculation of \u03BB \u03B4,2 (spin = [1,1,0]) ')
            plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
            plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)

            plt.savefig('fit_ene_5_1.png')
            plt.close()


            f = open('ene_5_2.dat','r')
            l = f.readlines()
            f.close

            x = []
            y = []

            for i in l:
                x.append(float(i.split()[0]))
                y.append(float(i.split()[1]))

            x = np.array(x)
            y = np.array(y)

            params = curve_fit(K, x, y)
            print("Fitting parameters for spin parallel to -110 (data from file ene_5_2.dat):")
            print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])
            r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
            print("R-squared =", r_squared)
            print("")

            if r_squared < 0.98:
                print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_5_2.png")
                print("")

            l2 = -params[0][1] / (2.0 * params[0][0])

            print("X minimum = -b/(2*a) =", l2)
            print("")

            lambda_delta = 2.0*((l1 -l2)/(l1+l2))


            plt.plot(x, y, 'bo', label='data in ene_5_2.dat')
            popt, pcov = curve_fit(K, x, y)
            t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
            plt.plot(t, K(t, *popt), 'r--', label='fit')
            plt.xlabel('Lattice distorsion along [1,1,0] direction (Å)')
            plt.ylabel('Energy (eV)')
            plt.legend()
            plt.title('Calculation of \u03BB \u03B4,2 (spin = [-1,1,0])')
            plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
            plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
            plt.savefig('fit_ene_5_2.png')
            plt.close()


            #make figure dE_5.png

            fig = 'dE_5.png'
            spin1 = '1,1,0'
            spin2 = '-1,1,0'
            dist = '1,1,0'
            tit = "Calculation of \u03BB \u03B4,2  "
            f1 = open('ene_5_1.dat','r')
            f2 = open('ene_5_2.dat','r')


            s1 = f1.readlines()
            s2 = f2.readlines()
            f1.close
            f2.close

            x = []
            y = []
            y2 = []
            for j in s1:
                x.append(float(j.split()[0]))
                y.append(float(j.split()[1]))

            for j in s2:
                y2.append(float(j.split()[1]))

            x = np.array(x)
            y = np.array(y)
            y2 = np.array(y2)

            plt.plot(x, (y2-y)*1e6, 'o-')

            ylabel ='E[' + str(spin2) + '] - E['+ str(spin1) + '] (\u03BCeV)'
            plt.ylabel(ylabel)
            label = "Lattice distorsion along [" + str(dist) + "] direction (Å)"
            plt.xlabel(label)
            plt.title(tit)
            plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
            plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
            plt.savefig(fig)
            plt.close()






        print(" ")
        print("----------------------------------------------")
        print("Anisotropic magnetostriction coefficients:")
        print("----------------------------------------------")
        print(" ")

        if 177 <= sg <= 194:
            print(" ")
            print("Using the convention in reference E.A. Clark et al., Phys. Rev. 138, A216 (1965):")
            print(" ")
            print("\u03BB 1\u03B1,2 =", lambda_alpha_1_2*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB 2\u03B1,2 =", lambda_alpha_2_2*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB \u03B3,2 =", lambda_gamma_2*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB \u03B5,2 =", lambda_epsilon_2*1e6,u'x 10\u207B\u2076')

            print(" ")
            print("...............")
            print(" ")
            print("Using the convention in reference W.P. Mason, Phys. Rev. 96, 302 (1954):")
            print(" ")
            print("\u03BBA =", (-lambda_alpha_1_2+0.5*lambda_gamma_2)*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BBB =", (-lambda_alpha_1_2-0.5*lambda_gamma_2)*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BBC =", -lambda_alpha_2_2*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BBD =", 0.5*(lambda_epsilon_2+0.5*(-lambda_alpha_1_2+0.5*lambda_gamma_2-lambda_alpha_2_2))*1e6,u'x 10\u207B\u2076')

            print(" ")
            print("...............")
            print(" ")
            print("Using the convention in reference R.R. Birss, Advances in Physics 8, 252 (1959):")
            print(" ")
            print("Q2 =", (-lambda_alpha_1_2-0.5*lambda_gamma_2)*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("Q4 =", (lambda_alpha_1_2+0.5*lambda_gamma_2-lambda_alpha_2_2)*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("Q6 =", 2*lambda_epsilon_2*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("Q8 =", lambda_gamma_2*1e6,u'x 10\u207B\u2076')

            print(" ")
            print("...............")
            print(" ")
            print("Using the convention in reference E.R. Callen et al., Phys. Rev. 139, A455 (1965):")
            print(" ")
            print("\u03BB 12\u03B1 =", (2/math.sqrt(3))*(2*lambda_alpha_1_2+lambda_alpha_2_2)*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB 22\u03B1 =", (1/math.sqrt(3))*(-lambda_alpha_1_2+lambda_alpha_2_2)*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB \u03B3 =", lambda_gamma_2*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB \u03B5 =", lambda_epsilon_2*1e6,u'x 10\u207B\u2076')

            if args.delas == True:

                print(" ")
                print(" ")
                print("----------------------------------------------")
                print("Calculation of magnetoelastic constants:")
                print("----------------------------------------------")
                print(" ")
                print("Reading the elastic tensor file =", str(args.elas[0]))
                print(" ")




                elasdat = open(args.elas[0],'r')
                elasline = elasdat.readlines()
                elasline0 = elasline[2]
                elasline1 = elasline[4]
                elasline2 = elasline[5]
                c11 = float(elasline0[0:8])
                c12 = float(elasline0[8:16])
                c13 = float(elasline0[16:24])
                c33 = float(elasline1[16:24])
                c44 = float(elasline2[24:32])


                elasdat.close()


                b21 = -(c11+c12)*lambda_alpha_1_2-c13*lambda_alpha_2_2

                b22 = -2*c13*lambda_alpha_1_2-c33*lambda_alpha_2_2

                b3 = -(c11-c12)*lambda_gamma_2

                b4 = -2*c44*lambda_epsilon_2


                print("c11 =", str(c11), 'GPa')
                print(" ")
                print("c12 =", str(c12), 'GPa')
                print(" ")
                print("c13 =", str(c13), 'GPa')
                print(" ")
                print("c33 =", str(c33), 'GPa')
                print(" ")
                print("c44 =", str(c44), 'GPa')
                print(" ")
                print("Warning: If these elastic constants are not the same as in the input elastic tensor file", str(args.elas[0]),", then check that the format of the elastic tensor is exactly the same as in the standard output file ELADAT generated by ELAS code (see Example folder)")
                print(" ")
                print(" ")
                print("Magnetoelastic constants:")
                print(" ")
                print("Using the convention in reference J.R. Cullen et al., in Materials, Science and Technology (VCH Publishings, 1994), pp.529-565:")
                print(" ")
                print("b21 =", str(b21), 'GPa')
                print(" ")
                print("b22 =", str(b22), 'GPa')
                print(" ")
                print("b3 =", str(b3), 'GPa')
                print(" ")
                print("b4 =", str(b4), 'GPa')
                print(" ")




        if 89 <= sg <= 142:

            print(" ")
            print("Using the convention in reference J.R. Cullen et al., in Materials, Science and Technology (VCH Publishings, 1994), pp.529-565:")
            print(" ")
            print("\u03BB 1\u03B1,2 =", lambda_alpha_1_2*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB 2\u03B1,2 =", lambda_alpha_2_2*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB \u03B3,2 =", lambda_gamma_2*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB \u03B5,2 =", lambda_epsilon_2*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB \u03B4,2 =", lambda_delta*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("...............")
            print(" ")
            print("Using the convention in reference W.P. Mason, Phys. Rev. 96, 302 (1954):")
            print(" ")
            print("\u03BB1 =", (-lambda_alpha_1_2+0.5*lambda_gamma_2)*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB2 =", 0.5*(lambda_epsilon_2-0.5*lambda_alpha_2_2-0.5*lambda_alpha_1_2+0.25*lambda_gamma_2)*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB3 =", (0.5*lambda_delta-lambda_alpha_1_2)*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB4 =", -lambda_alpha_2_2*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB5 =", (-lambda_alpha_1_2-0.5*lambda_gamma_2)*1e6,u'x 10\u207B\u2076')


            if args.delas == True:

                print(" ")
                print(" ")
                print("----------------------------------------------")
                print("Calculation of magnetoelastic constants:")
                print("----------------------------------------------")
                print(" ")
                print("Reading the elastic tensor file =", str(args.elas[0]))
                print(" ")




                elasdat = open(args.elas[0],'r')
                elasline = elasdat.readlines()
                elasline0 = elasline[2]
                elasline1 = elasline[4]
                elasline2 = elasline[5]
                elasline3 = elasline[7]
                c11 = float(elasline0[0:8])
                c12 = float(elasline0[8:16])
                c13 = float(elasline0[16:24])
                c33 = float(elasline1[16:24])
                c44 = float(elasline2[24:32])
                c66 = float(elasline3[40:48])

                elasdat.close()


                b21 = -(c11+c12)*lambda_alpha_1_2-c13*lambda_alpha_2_2

                b22 = -2*c13*lambda_alpha_1_2-c33*lambda_alpha_2_2

                b3 = -(c11-c12)*lambda_gamma_2

                b4 = -2*c44*lambda_epsilon_2

                b3p = -2*c66*lambda_delta


                print("c11 =", str(c11), 'GPa')
                print(" ")
                print("c12 =", str(c12), 'GPa')
                print(" ")
                print("c13 =", str(c13), 'GPa')
                print(" ")
                print("c33 =", str(c33), 'GPa')
                print(" ")
                print("c44 =", str(c44), 'GPa')
                print(" ")
                print("c66 =", str(c66), 'GPa')
                print(" ")
                print("Warning: If these elastic constants are not the same as in the input elastic tensor file", str(args.elas[0]),", then check that the format of the elastic tensor is exactly the same as in the standard output file ELADAT generated by ELAS code (see Example folder)")
                print(" ")
                print(" ")
                print("Magnetoelastic constants:")
                print(" ")
                print("Using the convention in reference J.R. Cullen et al., in Materials, Science and Technology (VCH Publishings, 1994), pp.529-565:")
                print(" ")
                print("b21 =", str(b21), 'GPa')
                print(" ")
                print("b22 =", str(b22), 'GPa')
                print(" ")
                print("b3 =", str(b3), 'GPa')
                print(" ")
                print("b'3 =", str(b3p), 'GPa')
                print(" ")
                print("b4 =", str(b4), 'GPa')
                print(" ")

                print("The equation of the magnetoelastic energy can be found in the User Manual")








#################################################################

##### TRIGONAL (I) ##### SG 149 - 167

#################################################################






elif 149 <= sg <= 167:
    print("Trigonal system")
    print("Point group =", str(pg))
    print("Number of anisotropic magnestostriction coefficients =", 6)

    if args.gen == True:

        # Convention: lattice vector a1 along x-axis
        angle = -math.pi*(60.0/180.0)
        dd = DeformStructureTransformation(deformation=((math.cos(angle), math.sin(angle), 0), (-math.sin(angle), math.cos(angle), 0), (0, 0, 1)))
        structure2b = dd.apply_transformation(structure2)


        for i in range(int(args.ndist[0])):


            strain1 = - float(args.strain[0])+2*(float(args.strain[0])/(float(args.ndist[0])-1))*i

            print("epsilon", strain1)


        #Generation POSCAR file

        #lambda_alpha_1_2


            a1 = 1.0 + strain1
            a2 = 1/math.sqrt(a1)
            a3 = a2
            dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
            structure3 = dd.apply_transformation(structure2b)
            pos_name3 = "POSCAR_1_" + str(i+1)

            structure33 = Poscar(structure3)
            structure33.write_file(filename = pos_name3,significant_figures=16)

        #lambda_alpha_2_2


            a3 = 1.0 + strain1
            a2 = 1/math.sqrt(a3)
            a1 = a2
            dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
            structure4 = dd.apply_transformation(structure2b)
            pos_name4 = "POSCAR_2_" + str(i+1)

            structure44 = Poscar(structure4)
            structure44.write_file(filename = pos_name4,significant_figures=16)

        #lambda_gamma_1

            a1 = 1.0 + strain1
            a2 = 1/math.sqrt(a1)
            a3 = a2
            dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
            structure5 = dd.apply_transformation(structure2b)
            pos_name5 = "POSCAR_3_" + str(i+1)

            structure55 = Poscar(structure5)
            structure55.write_file(filename = pos_name5,significant_figures=16)


        #lambda_gamma_2

            latt_par = structure2b.lattice.matrix
            
            const = (1/(1-(strain1*0.5)**2))**(1/3)
            
            a11 = const
            a12 = 0.0
            a13 = const*strain1*0.5*(latt_par[2][2]/latt_par[0][0])
            a21 = 0.0
            a22 = const
            a23 = 0.0
            a31 = const*strain1*0.5*(latt_par[0][0]/latt_par[2][2])
            a32 = 0.0
            a33 = const
            

            cc = DeformStructureTransformation(deformation=((a11, a12, a13), (a21, a22, a23), (a31, a32, a33)))
            structure6 = cc.apply_transformation(structure2b)
            pos_name6 = "POSCAR_4_" + str(i+1)

            structure66 = Poscar(structure6)
            structure66.write_file(filename = pos_name6,significant_figures=16)

        #lambda_1_2

            pos_name7 = "POSCAR_5_" + str(i+1)

            structure77 = Poscar(structure6)
            structure77.write_file(filename = pos_name7,significant_figures=16)


        #lambda_2_1

            pos_name8 = "POSCAR_6_" + str(i+1)

            structure88 = Poscar(structure6)
            structure88.write_file(filename = pos_name8,significant_figures=16)


    # INCAR_1_1 m=0,0,1

        path_inc_ncl_1_1 = 'INCAR_1_1'
        inc_ncl_1_1 = open(path_inc_ncl_1_1,'w')
        inc_ncl_list_1_1 = generator.inc_ncl_list[:]
        inc_ncl_list_1_1 += ['SAXIS = 0 0 1.0\n']

        for j in range(len(inc_ncl_list_1_1)):
            inc_ncl_1_1.write(str(inc_ncl_list_1_1[j]))

        inc_ncl_1_1.close()


    # INCAR_1_2 m=1,1,0

        path_inc_ncl_1_2 = 'INCAR_1_2'
        inc_ncl_1_2 = open(path_inc_ncl_1_2,'w')
        inc_ncl_list_1_2 = generator.inc_ncl_list[:]
        inc_ncl_list_1_2 += ['SAXIS = 1.0 1.0 0.0\n']

        for j in range(len(inc_ncl_list_1_2)):
            inc_ncl_1_2.write(str(inc_ncl_list_1_2[j]))

        inc_ncl_1_2.close()


    # INCAR_2_1 m=0,0,1

        path_inc_ncl_2_1 = 'INCAR_2_1'
        inc_ncl_2_1 = open(path_inc_ncl_2_1,'w')
        inc_ncl_list_2_1 = generator.inc_ncl_list[:]
        inc_ncl_list_2_1 += ['SAXIS = 0 0 1.0\n']

        for j in range(len(inc_ncl_list_2_1)):
            inc_ncl_2_1.write(str(inc_ncl_list_2_1[j]))

        inc_ncl_2_1.close()


    # INCAR_2_2 m=1,0,0

        path_inc_ncl_2_2 = 'INCAR_2_2'
        inc_ncl_2_2 = open(path_inc_ncl_2_2,'w')
        inc_ncl_list_2_2 = generator.inc_ncl_list[:]
        inc_ncl_list_2_2 += ['SAXIS = 1.0 0.0 0.0\n']

        for j in range(len(inc_ncl_list_2_2)):
            inc_ncl_2_2.write(str(inc_ncl_list_2_2[j]))

        inc_ncl_2_2.close()


    # INCAR_3_1 m=1,0,0

        path_inc_ncl_3_1 = 'INCAR_3_1'
        inc_ncl_3_1 = open(path_inc_ncl_3_1,'w')
        inc_ncl_list_3_1 = generator.inc_ncl_list[:]
        inc_ncl_list_3_1 += ['SAXIS = 1.0 0.0 0.0\n']

        for j in range(len(inc_ncl_list_3_1)):
            inc_ncl_3_1.write(str(inc_ncl_list_3_1[j]))

        inc_ncl_3_1.close()


    # INCAR_3_2 m=0,1,0

        path_inc_ncl_3_2 = 'INCAR_3_2'
        inc_ncl_3_2 = open(path_inc_ncl_3_2,'w')
        inc_ncl_list_3_2 = generator.inc_ncl_list[:]
        inc_ncl_list_3_2 += ['SAXIS = 0.0 1.0 0.0\n']

        for j in range(len(inc_ncl_list_3_2)):
            inc_ncl_3_2.write(str(inc_ncl_list_3_2[j]))

        inc_ncl_3_2.close()



     # INCAR_4_1 m=1,0,1

        path_inc_ncl_4_1 = 'INCAR_4_1'
        inc_ncl_4_1 = open(path_inc_ncl_4_1,'w')
        inc_ncl_list_4_1 = generator.inc_ncl_list[:]
        inc_ncl_list_4_1 += ['SAXIS = 1.0 0.0 1.0\n']

        for j in range(len(inc_ncl_list_4_1)):
            inc_ncl_4_1.write(str(inc_ncl_list_4_1[j]))

        inc_ncl_4_1.close()


    # INCAR_4_2 m=1,0,-1

        path_inc_ncl_4_2 = 'INCAR_4_2'
        inc_ncl_4_2 = open(path_inc_ncl_4_2,'w')
        inc_ncl_list_4_2 = generator.inc_ncl_list[:]
        inc_ncl_list_4_2 += ['SAXIS = 1.0 0.0 -1.0\n']

        for j in range(len(inc_ncl_list_4_2)):
            inc_ncl_4_2.write(str(inc_ncl_list_4_2[j]))

        inc_ncl_4_2.close()



     # INCAR_5_1 m=0,1,1

        path_inc_ncl_5_1 = 'INCAR_5_1'
        inc_ncl_5_1 = open(path_inc_ncl_5_1,'w')
        inc_ncl_list_5_1 = generator.inc_ncl_list[:]
        inc_ncl_list_5_1 += ['SAXIS = 0 1.0 1.0\n']

        for j in range(len(inc_ncl_list_5_1)):
            inc_ncl_5_1.write(str(inc_ncl_list_5_1[j]))

        inc_ncl_5_1.close()


    # INCAR_5_2 m=0,1,-1

        path_inc_ncl_5_2 = 'INCAR_5_2'
        inc_ncl_5_2 = open(path_inc_ncl_5_2,'w')
        inc_ncl_list_5_2 = generator.inc_ncl_list[:]
        inc_ncl_list_5_2 += ['SAXIS = 0 1.0 -1.0\n']

        for j in range(len(inc_ncl_list_5_2)):
            inc_ncl_5_2.write(str(inc_ncl_list_5_2[j]))

        inc_ncl_5_2.close()


     # INCAR_6_1 m=1,1,0

        path_inc_ncl_6_1 = 'INCAR_6_1'
        inc_ncl_6_1 = open(path_inc_ncl_6_1,'w')
        inc_ncl_list_6_1 = generator.inc_ncl_list[:]
        inc_ncl_list_6_1 += ['SAXIS = 1.0 1.0 0\n']

        for j in range(len(inc_ncl_list_6_1)):
            inc_ncl_6_1.write(str(inc_ncl_list_6_1[j]))

        inc_ncl_6_1.close()


    # INCAR_6_2 m=1,-1,0

        path_inc_ncl_6_2 = 'INCAR_6_2'
        inc_ncl_6_2 = open(path_inc_ncl_6_2,'w')
        inc_ncl_list_6_2 = generator.inc_ncl_list[:]
        inc_ncl_list_6_2 += ['SAXIS = 1.0 -1.0 0\n']

        for j in range(len(inc_ncl_list_6_2)):
            inc_ncl_6_2.write(str(inc_ncl_list_6_2[j]))

        inc_ncl_6_2.close()








    # Derivation of magnetostriction coefficients:

    if args.der == True:

        for j in range(1,7):

            for k in range(1,3):

                path_dat = "ene_" + str(j) + "_" + str(k) + ".dat"
                dat = open(path_dat,'w')


                for i in range(int(args.ndist[0])):

                    pos_name = "POSCAR_" + str(j) + "_" + str(i+1)

                    struct = Structure.from_file(pos_name)

                    latt = struct.lattice.matrix

                    if j == 1:
                        var1 = latt[0][0]
                    elif j == 2:
                        var1 = latt[2][2]
                    elif j == 3:
                        var1 = latt[0][0]
                    elif j == 4:
                        var1 = math.sqrt((latt[0][0]+latt[2][0])**2+(latt[0][1]+latt[2][1])**2+(latt[0][2]+latt[2][2])**2)
                    elif j == 5:
                        var1 = math.sqrt((latt[0][0]+latt[2][0])**2+(latt[0][1]+latt[2][1])**2+(latt[0][2]+latt[2][2])**2)
                    elif j == 6:
                        var1 = math.sqrt((latt[0][0]+latt[2][0])**2+(latt[0][1]+latt[2][1])**2+(latt[0][2]+latt[2][2])**2)

                    path_osz = "OSZICAR_" + str(j) + "_" + str(i+1) + "_" + str(k)
                    osz = open(path_osz,'r')
                    ene0 = osz.readlines()
                    ene1 = ene0[len(ene0)-2]
                    ene2 = ene1[11:32]

                    osz.close()

                    dat.write(repr(var1))
                    dat.write('  ')
                    dat.write(str(ene2))
                    dat.write('\n')


                dat.close()


       # fitting and plot


        def K(x,a,b,c):
            return a*x*x+b*x+c


        print("")
        print("Fit of quadratic function f(x)=a*x\u00B2+b*x+c to energy vs distortion data")
        print("")
        print("-------------------------")
        print("Calculation of \u03BB \u03B11,2:")
        print("-------------------------")
        print(" ")
        print('Lattice distorsion along [1,0,0] direction')
        print("")

        f = open('ene_1_1.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []
        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)

        print("Fitting parameters for spin parallel to 001 (data from file ene_1_1.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_1_1.png")
            print("")

        l1 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l1)
        print("")


        plt.plot(x, y, 'bo', label='data in ene_1_1.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Lattice distorsion along [1,0,0] direction (Å)')
        plt.title('Calculation of \u03BB \u03B11,2 (spin = [0,0,1]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)

        plt.savefig('fit_ene_1_1.png')
        plt.close()

        f = open('ene_1_2.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []
        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)
        print("Fitting parameters for spin parallel to 110 (data from file ene_1_2.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_1_2.png")
            print("")

        l2 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l2)
        print("")


        lambda_alpha_1_2 = 2.0*((l1 -l2)/(l1+l2))


        plt.plot(x, y, 'bo', label='data in ene_1_2.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Lattice distorsion along [0,0,1] direction (Å)')
        plt.title('Calculation of \u03BB \u03B11,2 (spin = [1,1,0]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)

        plt.savefig('fit_ene_1_2.png')
        plt.close()


        #make figure dE_1.png

        fig = 'dE_1.png'
        spin1 = '0,0,1'
        spin2 = '1,1,0'
        dist = '0,0,1'
        tit = "Calculation of \u03BB \u03B11,2  "
        f1 = open('ene_1_1.dat','r')
        f2 = open('ene_1_2.dat','r')


        s1 = f1.readlines()
        s2 = f2.readlines()
        f1.close
        f2.close

        x = []
        y = []
        y2 = []
        for j in s1:
            x.append(float(j.split()[0]))
            y.append(float(j.split()[1]))

        for j in s2:
                y2.append(float(j.split()[1]))

        x = np.array(x)
        y = np.array(y)
        y2 = np.array(y2)

        plt.plot(x, (y2-y)*1e6, 'o-')

        ylabel ='E[' + str(spin2) + '] - E['+ str(spin1) + '] (\u03BCeV)'
        plt.ylabel(ylabel)
        label = "Lattice distorsion along [" + str(dist) + "] direction (Å)"
        plt.xlabel(label)
        plt.title(tit)
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig(fig)
        plt.close()




        print(" ")
        print("-------------------------")
        print("Calculation of \u03BB \u03B12,2:")
        print("-------------------------")
        print(" ")
        print('Lattice distorsion along [0,0,1] direction')
        print("")

        f = open('ene_2_1.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []
        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)
        print("Fitting parameters for spin parallel to 001 (data from file ene_2_1.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_2_1.png")
            print("")

        l1 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l1)
        print("")


        plt.plot(x, y, 'bo', label='data in ene_2_1.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Lattice distorsion along [0,0,1] direction (Å)')
        plt.title('Calculation of \u03BB \u03B12,2 (spin = [0,0,1]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)

        plt.savefig('fit_ene_2_1.png')
        plt.close()


        f = open('ene_2_2.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []

        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)
        print("Fitting parameters for spin parallel to 100 (data from file ene_2_2.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_2_2.png")
            print("")

        l2 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l2)
        print("")

        lambda_alpha_2_2 = 2.0*((l1 -l2)/(l1+l2))


        plt.plot(x, y, 'bo', label='data in ene_2_2.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.xlabel('Lattice distorsion along [0,0,1] direction (Å)')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.title('Calculation of \u03BB \u03B12,2 (spin = [1,0,0]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig('fit_ene_2_2.png')
        plt.close()



        #make figure dE_2.png

        fig = 'dE_2.png'
        spin1 = '0,0,1'
        spin2 = '1,0,0'
        dist = '0,0,1'
        tit = "Calculation of \u03BB \u03B12,2  "
        f1 = open('ene_2_1.dat','r')
        f2 = open('ene_2_2.dat','r')


        s1 = f1.readlines()
        s2 = f2.readlines()
        f1.close
        f2.close

        x = []
        y = []
        y2 = []
        for j in s1:
            x.append(float(j.split()[0]))
            y.append(float(j.split()[1]))

        for j in s2:
                y2.append(float(j.split()[1]))

        x = np.array(x)
        y = np.array(y)
        y2 = np.array(y2)

        plt.plot(x, (y2-y)*1e6, 'o-')

        ylabel ='E[' + str(spin2) + '] - E['+ str(spin1) + '] (\u03BCeV)'
        plt.ylabel(ylabel)
        label = "Lattice distorsion along [" + str(dist) + "] direction (Å)"
        plt.xlabel(label)
        plt.title(tit)
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig(fig)
        plt.close()




        print(" ")
        print("-------------------------")
        print("Calculation of \u03BB \u02631:")
        print("-------------------------")
        print(" ")
        print('Lattice distorsion along [1,0,0] direction')
        print("")

        f = open('ene_3_1.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []
        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)
        print("Fitting parameters for spin parallel to 100 (data from file ene_3_1.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_3_1.png")
            print("")

        l1 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l1)
        print("")


        plt.plot(x, y, 'bo', label='data in ene_3_1.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Lattice distorsion along [1,0,0] direction (Å)')
        plt.title('Calculation of \u03BB \u02631 (spin = [1,0,0]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)

        plt.savefig('fit_ene_3_1.png')
        plt.close()


        f = open('ene_3_2.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []

        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)
        print("Fitting parameters for spin parallel to 010 (data from file ene_3_2.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_3_2.png")
            print("")

        l2 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l2)
        print("")

        lambda_gamma_1 = 2.0*((l1 -l2)/(l1+l2))


        plt.plot(x, y, 'bo', label='data in ene_3_2.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.xlabel('Lattice distorsion along [1,0,0] direction (Å)')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.title('Calculation of \u03BB \u02631 (spin = [0,1,0]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig('fit_ene_3_2.png')
        plt.close()


        #make figure dE_3.png

        fig = 'dE_3.png'
        spin1 = '1,0,0'
        spin2 = '0,1,0'
        dist = '1,0,0'
        tit = "Calculation of \u03BB \u02631  "
        f1 = open('ene_3_1.dat','r')
        f2 = open('ene_3_2.dat','r')


        s1 = f1.readlines()
        s2 = f2.readlines()
        f1.close
        f2.close

        x = []
        y = []
        y2 = []
        for j in s1:
            x.append(float(j.split()[0]))
            y.append(float(j.split()[1]))

        for j in s2:
                y2.append(float(j.split()[1]))

        x = np.array(x)
        y = np.array(y)
        y2 = np.array(y2)

        plt.plot(x, (y2-y)*1e6, 'o-')

        ylabel ='E[' + str(spin2) + '] - E['+ str(spin1) + '] (\u03BCeV)'
        plt.ylabel(ylabel)
        label = "Lattice distorsion along [" + str(dist) + "] direction (Å)"
        plt.xlabel(label)
        plt.title(tit)
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig(fig)
        plt.close()







        print(" ")
        print("-------------------------")
        print("Calculation of \u03BB \u02632:")
        print("-------------------------")
        print(" ")
        print('Lattice distorsion along [1,0,1] direction')
        print("")

        f = open('ene_4_1.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []
        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)
        print("Fitting parameters for spin parallel to 101 (data from file ene_4_1.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_4_1.png")
            print("")

        l1 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l1)
        print("")


        plt.plot(x, y, 'bo', label='data in ene_4_1.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Lattice distorsion along [1,0,1] direction (Å)')
        plt.title('Calculation of \u03BB \u02632 (spin = [1,0,1]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)

        plt.savefig('fit_ene_4_1.png')
        plt.close()


        f = open('ene_4_2.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []

        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)
        print("Fitting parameters for spin parallel to 10-1 (data from file ene_4_2.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_4_2.png")
            print("")

        l2 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l2)
        print("")

        
        aa0 = SpacegroupAnalyzer(structure0,symprec=sym1, angle_tolerance=sym2)
        structure1 = aa0.get_conventional_standard_structure(international_monoclinic=True)

        bb0 = ConventionalCellTransformation(symprec=sym1, angle_tolerance=sym2, international_monoclinic=True)
        structure2 = bb0.apply_transformation(structure1) 
        
        # Convention: lattice vector a1 along x-axis
        angle = -math.pi*(60.0/180.0)
        dd = DeformStructureTransformation(deformation=((math.cos(angle), math.sin(angle), 0), (-math.sin(angle), math.cos(angle), 0), (0, 0, 1)))
        structure2b = dd.apply_transformation(structure2)
        
        
        latt_par = structure2b.lattice.matrix
            
        latt_a = latt_par[0][0]
        latt_c = latt_par[2][2]
        
        
        eta_par = (latt_a*latt_c)/(latt_a**2+latt_c**2)
        
        
        lambda_gamma_2 = 2.0*((l1 -l2)/(l1+l2))*(1.0/eta_par)


        plt.plot(x, y, 'bo', label='data in ene_4_2.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.xlabel('Lattice distorsion along [1,0,1] direction (Å)')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.title('Calculation of \u03BB \u02632 (spin = [1,0,-1]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig('fit_ene_4_2.png')
        plt.close()


       #make figure dE_4.png

        fig = 'dE_4.png'
        spin1 = '1,0,1'
        spin2 = '1,0,-1'
        dist = '1,0,1'
        tit = "Calculation of \u03BB \u02632  "
        f1 = open('ene_4_1.dat','r')
        f2 = open('ene_4_2.dat','r')


        s1 = f1.readlines()
        s2 = f2.readlines()
        f1.close
        f2.close

        x = []
        y = []
        y2 = []
        for j in s1:
            x.append(float(j.split()[0]))
            y.append(float(j.split()[1]))

        for j in s2:
                y2.append(float(j.split()[1]))

        x = np.array(x)
        y = np.array(y)
        y2 = np.array(y2)

        plt.plot(x, (y2-y)*1e6, 'o-')

        ylabel ='E[' + str(spin2) + '] - E['+ str(spin1) + '] (\u03BCeV)'
        plt.ylabel(ylabel)
        label = "Lattice distorsion along [" + str(dist) + "] direction (Å)"
        plt.xlabel(label)
        plt.title(tit)
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig(fig)
        plt.close()






        print(" ")
        print("-------------------------")
        print("Calculation of \u03BB 12:")
        print("-------------------------")
        print(" ")
        print('Lattice distorsion along [1,0,1] direction')
        print("")

        f = open('ene_5_1.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []
        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)
        print("Fitting parameters for spin parallel to 011 (data from file ene_5_1.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_5_1.png")
            print("")

        l1 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l1)
        print("")


        plt.plot(x, y, 'bo', label='data in ene_5_1.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Lattice distorsion along [1,0,1] direction (Å)')
        plt.title('Calculation of \u03BB 12 (spin = [0,1,1]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)

        plt.savefig('fit_ene_5_1.png')
        plt.close()


        f = open('ene_5_2.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []

        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)
        print("Fitting parameters for spin parallel to 01-1 (data from file ene_5_2.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_5_2.png")
            print("")

        l2 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l2)
        print("")

        eta_par = (latt_a*latt_c)/(2.0*(latt_a**2+latt_c**2))
        
        
        lambda_1_2 = 2.0*((l1 -l2)/(l1+l2))*(1.0/eta_par)


        plt.plot(x, y, 'bo', label='data in ene_5_2.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.xlabel('Lattice distorsion along [1,0,1] direction (Å)')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.title('Calculation of \u03BB 12 (spin = [0,1,-1]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig('fit_ene_5_2.png')
        plt.close()



        #make figure dE_5.png

        fig = 'dE_5.png'
        spin1 = '0,1,1'
        spin2 = '0,1,-1'
        dist = '1,0,1'
        tit = "Calculation of \u03BB 12  "
        f1 = open('ene_5_1.dat','r')
        f2 = open('ene_5_2.dat','r')


        s1 = f1.readlines()
        s2 = f2.readlines()
        f1.close
        f2.close

        x = []
        y = []
        y2 = []
        for j in s1:
            x.append(float(j.split()[0]))
            y.append(float(j.split()[1]))

        for j in s2:
                y2.append(float(j.split()[1]))

        x = np.array(x)
        y = np.array(y)
        y2 = np.array(y2)

        plt.plot(x, (y2-y)*1e6, 'o-')

        ylabel ='E[' + str(spin2) + '] - E['+ str(spin1) + '] (\u03BCeV)'
        plt.ylabel(ylabel)
        label = "Lattice distorsion along [" + str(dist) + "] direction (Å)"
        plt.xlabel(label)
        plt.title(tit)
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig(fig)
        plt.close()



        print(" ")
        print("-------------------------")
        print("Calculation of \u03BB 21:")
        print("-------------------------")
        print(" ")
        print('Lattice distorsion along [1,0,1] direction')
        print("")

        f = open('ene_6_1.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []
        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)
        print("Fitting parameters for spin parallel to 110 (data from file ene_6_1.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_6_1.png")
            print("")

        l1 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l1)
        print("")


        plt.plot(x, y, 'bo', label='data in ene_6_1.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Lattice distorsion along [1,0,1] direction (Å)')
        plt.title('Calculation of \u03BB 21 (spin = [1,1,0]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)

        plt.savefig('fit_ene_6_1.png')
        plt.close()


        f = open('ene_6_2.dat','r')
        l = f.readlines()
        f.close

        x = []
        y = []

        for i in l:
            x.append(float(i.split()[0]))
            y.append(float(i.split()[1]))

        x = np.array(x)
        y = np.array(y)

        params = curve_fit(K, x, y)
        print("Fitting parameters for spin parallel to 1-10 (data from file ene_6_2.dat):")
        print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_6_2.png")
            print("")

        l2 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -b/(2*a) =", l2)
        print("")

        eta_par = (latt_a*latt_c)/(latt_a**2+latt_c**2)
        
        lambda_2_1 = 2.0*((l1 -l2)/(l1+l2))*(1.0/eta_par)


        plt.plot(x, y, 'bo', label='data in ene_6_2.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.xlabel('Lattice distorsion along [1,0,1] direction (Å)')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.title('Calculation of \u03BB 21 (spin = [1,-1,0]) ')
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig('fit_ene_6_2.png')
        plt.close()



        #make figure dE_6.png

        fig = 'dE_6.png'
        spin1 = '1,1,0'
        spin2 = '1,-1,0'
        dist = '1,0,1'
        tit = "Calculation of \u03BB 21 "
        f1 = open('ene_6_1.dat','r')
        f2 = open('ene_6_2.dat','r')


        s1 = f1.readlines()
        s2 = f2.readlines()
        f1.close
        f2.close

        x = []
        y = []
        y2 = []
        for j in s1:
            x.append(float(j.split()[0]))
            y.append(float(j.split()[1]))

        for j in s2:
                y2.append(float(j.split()[1]))

        x = np.array(x)
        y = np.array(y)
        y2 = np.array(y2)

        plt.plot(x, (y2-y)*1e6, 'o-')

        ylabel ='E[' + str(spin2) + '] - E['+ str(spin1) + '] (\u03BCeV)'
        plt.ylabel(ylabel)
        label = "Lattice distorsion along [" + str(dist) + "] direction (Å)"
        plt.xlabel(label)
        plt.title(tit)
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig(fig)
        plt.close()



        print(" ")
        print("----------------------------------------------")
        print("Anisotropic magnetostriction coefficients:")
        print("----------------------------------------------")
        print(" ")
        print("Using the convention in reference J.R. Cullen et al., in Materials, Science and Technology (VCH Publishings, 1994), pp.529-565:")
        print(" ")
        print("\u03BB \u03B11,2  =", lambda_alpha_1_2*1e6,u'x 10\u207B\u2076')
        print(" ")
        print("\u03BB \u03B12,2 =", lambda_alpha_2_2*1e6,u'x 10\u207B\u2076')
        print(" ")
        print("\u03BB \u02631 =", lambda_gamma_1*1e6,u'x 10\u207B\u2076')
        print(" ")
        print("\u03BB \u02632 =", lambda_gamma_2*1e6,u'x 10\u207B\u2076')
        print(" ")
        print("\u03BB 12 =", lambda_1_2*1e6,u'x 10\u207B\u2076')
        print(" ")
        print("\u03BB 21 =", lambda_2_1*1e6,u'x 10\u207B\u2076')


        if args.delas == True:

                print(" ")
                print(" ")
                print("----------------------------------------------")
                print("Calculation of magnetoelastic constants:")
                print("----------------------------------------------")
                print(" ")
                print("Reading the elastic tensor file =", str(args.elas[0]))
                print(" ")




                elasdat = open(args.elas[0],'r')
                elasline = elasdat.readlines()
                elasline0 = elasline[2]
                elasline1 = elasline[4]
                elasline2 = elasline[5]
                elasline3 = elasline[7]
                c11 = float(elasline0[0:8])
                c12 = float(elasline0[8:16])
                c13 = float(elasline0[16:24])
                c14 = float(elasline0[24:32])
                c33 = float(elasline1[16:24])
                c44 = float(elasline2[24:32])
                c66 = float(elasline3[40:48])

                elasdat.close()


                b21 = -(c11+c12)*lambda_alpha_1_2-c13*lambda_alpha_2_2

                b22 = -2*c13*lambda_alpha_1_2-c33*lambda_alpha_2_2

                b3 = c14*lambda_2_1+0.5*(-c11+c12)*lambda_gamma_1

                b4 = -c14*lambda_1_2+c44*lambda_gamma_2

                b14 =  c44*lambda_2_1-c14*lambda_gamma_1

                b34 = 0.5*(-c11+c12)*lambda_1_2+c14*lambda_gamma_2



                print("c11 =", str(c11), 'GPa')
                print(" ")
                print("c12 =", str(c12), 'GPa')
                print(" ")
                print("c13 =", str(c13), 'GPa')
                print(" ")
                print("c14 =", str(c14), 'GPa')
                print(" ")
                print("c33 =", str(c33), 'GPa')
                print(" ")
                print("c44 =", str(c44), 'GPa')
                print(" ")
                print("c66 =", str(c66), 'GPa')
                print(" ")
                print("Warning: If these elastic constants are not the same as in the input elastic tensor file", str(args.elas[0]),", then check that the format of the elastic tensor is exactly the same as in the standard output file ELADAT generated by ELAS code (see Example folder)")
                print(" ")
                print(" ")
                print("Magnetoelastic constants:")
                print(" ")
                print("Using the convention in reference J.R. Cullen et al., in Materials, Science and Technology (VCH Publishings, 1994), pp.529-565:")
                print(" ")
                print("b21 =", str(b21), 'GPa')
                print(" ")
                print("b22 =", str(b22), 'GPa')
                print(" ")
                print("b3 =", str(b3), 'GPa')
                print(" ")
                print("b4 =", str(b4), 'GPa')
                print(" ")
                print("b14 =", str(b14), 'GPa')
                print(" ")
                print("b34 =", str(b34), 'GPa')
                print(" ")

                print("The equation of the magnetoelastic energy can be found in the User Manual")






#################################################################

##### ORTHORHOMBIC #### SG 16 - 74

#################################################################


elif 16 <= sg <= 74:
    print("Orthorhombic system")
    print("Point group =", str(pg))
    print("Number of anisotropic magnestostriction coefficients =", 9)

    if args.gen == True:


        # AELAS and IEEE lattice convention: c<a<b


        latt0 = structure2.lattice.matrix
        coordsnew = np.zeros((len(structure2.species), 3))

        for i in range(len(structure2.species)):
            coordsnew[i][0] = float(structure2.frac_coords[i][1])
            coordsnew[i][1] = float(structure2.frac_coords[i][2])
            coordsnew[i][2] = float(structure2.frac_coords[i][0])


        lattice = Lattice.from_parameters(a=latt0[1][1], b=latt0[2][2], c=latt0[0][0], alpha=90, beta=90, gamma=90)
        structure2b = Structure(lattice, structure2.species, coordsnew)


        for i in range(int(args.ndist[0])):


            strain1 = - float(args.strain[0])+2*(float(args.strain[0])/(float(args.ndist[0])-1))*i

            print("epsilon", strain1)


        #Generation POSCAR file

        #lambda_1


            a1 = 1.0 + strain1
            a2 = 1/math.sqrt(a1)
            a3 = a2
            dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
            structure3 = dd.apply_transformation(structure2b)
            pos_name3 = "POSCAR_1_" + str(i+1)

            structure33 = Poscar(structure3)
            structure33.write_file(filename = pos_name3,significant_figures=16)

        #lambda_2


            a1 = 1.0 + strain1
            a2 = 1/math.sqrt(a1)
            a3 = a2
            dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
            structure4 = dd.apply_transformation(structure2b)
            pos_name4 = "POSCAR_2_" + str(i+1)

            structure44 = Poscar(structure4)
            structure44.write_file(filename = pos_name4,significant_figures=16)

        #lambda_3

            a2 = 1.0 + strain1
            a1 = 1/math.sqrt(a2)
            a3 = a1
            dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
            structure5 = dd.apply_transformation(structure2b)
            pos_name5 = "POSCAR_3_" + str(i+1)

            structure55 = Poscar(structure5)
            structure55.write_file(filename = pos_name5,significant_figures=16)

        #lambda_4

            a2 = 1.0 + strain1
            a1 = 1/math.sqrt(a2)
            a3 = a1
            dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
            structure6 = dd.apply_transformation(structure2b)
            pos_name6 = "POSCAR_4_" + str(i+1)

            structure66 = Poscar(structure6)
            structure66.write_file(filename = pos_name6,significant_figures=16)


       #lambda_5

            a3 = 1.0 + strain1
            a1 = 1/math.sqrt(a3)
            a2 = a1
            dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
            structure7 = dd.apply_transformation(structure2b)
            pos_name7 = "POSCAR_5_" + str(i+1)

            structure77 = Poscar(structure7)
            structure77.write_file(filename = pos_name7,significant_figures=16)

       #lambda_6

            a3 = 1.0 + strain1
            a1 = 1/math.sqrt(a3)
            a2 = a1
            dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
            structure8 = dd.apply_transformation(structure2b)
            pos_name8 = "POSCAR_6_" + str(i+1)

            structure88 = Poscar(structure8)
            structure88.write_file(filename = pos_name8,significant_figures=16)

        #lambda_7

            latt_par = structure2b.lattice.matrix
            
            const = (1/(1-(strain1*0.5)**2))**(1/3)
            
            a11 = const
            a12 = const*strain1*0.5*(latt_par[1][1]/latt_par[0][0])
            a13 = 0.0
            a21 = const*strain1*0.5*(latt_par[0][0]/latt_par[1][1])
            a22 = const
            a23 = 0.0
            a31 = 0.0
            a32 = 0.0
            a33 = const 
            

            cc = DeformStructureTransformation(deformation=((a11, a12, a13), (a21, a22, a23), (a31, a32, a33)))
            structure9 = cc.apply_transformation(structure2b)
            pos_name9 = "POSCAR_7_" + str(i+1)

            structure99 = Poscar(structure9)
            structure99.write_file(filename = pos_name9,significant_figures=16)

        #lambda_8
                 
            
            const = (1/(1-(strain1*0.5)**2))**(1/3)

            a11 = const
            a12 = 0.0
            a13 = const*strain1*0.5*(latt_par[2][2]/latt_par[0][0])
            a21 = 0.0
            a22 = const
            a23 = 0.0
            a31 = const*strain1*0.5*(latt_par[0][0]/latt_par[2][2])
            a32 = 0.0
            a33 = const

            cc = DeformStructureTransformation(deformation=((a11, a12, a13), (a21, a22, a23), (a31, a32, a33)))
            structure10 = cc.apply_transformation(structure2b)
            pos_name10 = "POSCAR_8_" + str(i+1)

            structure1010 = Poscar(structure10)
            structure1010.write_file(filename = pos_name10,significant_figures=16)


        #lambda_9

            const = (1/(1-(strain1*0.5)**2))**(1/3)

            a11 = const
            a12 = 0.0
            a13 = 0.0
            a21 = 0.0
            a22 = const
            a23 = const*strain1*0.5*(latt_par[2][2]/latt_par[1][1])
            a31 = 0.0
            a32 = const*strain1*0.5*(latt_par[1][1]/latt_par[2][2])
            a33 = const

            cc = DeformStructureTransformation(deformation=((a11, a12, a13), (a21, a22, a23), (a31, a32, a33)))
            structure11 = cc.apply_transformation(structure2b)
            pos_name11 = "POSCAR_9_" + str(i+1)

            structure1111 = Poscar(structure11)
            structure1111.write_file(filename = pos_name11,significant_figures=16)



    # INCAR_1_1 m=1,0,0

        path_inc_ncl_1_1 = 'INCAR_1_1'
        inc_ncl_1_1 = open(path_inc_ncl_1_1,'w')
        inc_ncl_list_1_1 = generator.inc_ncl_list[:]
        inc_ncl_list_1_1 += ['SAXIS = 1.0 0 0\n']

        for j in range(len(inc_ncl_list_1_1)):
            inc_ncl_1_1.write(str(inc_ncl_list_1_1[j]))

        inc_ncl_1_1.close()


    # INCAR_1_2 m=0,0,1

        path_inc_ncl_1_2 = 'INCAR_1_2'
        inc_ncl_1_2 = open(path_inc_ncl_1_2,'w')
        inc_ncl_list_1_2 = generator.inc_ncl_list[:]
        inc_ncl_list_1_2 += ['SAXIS = 0.0 0.0 1.0\n']

        for j in range(len(inc_ncl_list_1_2)):
            inc_ncl_1_2.write(str(inc_ncl_list_1_2[j]))

        inc_ncl_1_2.close()


    # INCAR_2_1 m=0,1,0

        path_inc_ncl_2_1 = 'INCAR_2_1'
        inc_ncl_2_1 = open(path_inc_ncl_2_1,'w')
        inc_ncl_list_2_1 = generator.inc_ncl_list[:]
        inc_ncl_list_2_1 += ['SAXIS = 0 1.0 0\n']

        for j in range(len(inc_ncl_list_2_1)):
            inc_ncl_2_1.write(str(inc_ncl_list_2_1[j]))

        inc_ncl_2_1.close()


    # INCAR_2_2 m=0,0,1

        path_inc_ncl_2_2 = 'INCAR_2_2'
        inc_ncl_2_2 = open(path_inc_ncl_2_2,'w')
        inc_ncl_list_2_2 = generator.inc_ncl_list[:]
        inc_ncl_list_2_2 += ['SAXIS = 0.0 0.0 1.0\n']

        for j in range(len(inc_ncl_list_2_2)):
            inc_ncl_2_2.write(str(inc_ncl_list_2_2[j]))

        inc_ncl_2_2.close()


    # INCAR_3_1 m=1,0,0

        path_inc_ncl_3_1 = 'INCAR_3_1'
        inc_ncl_3_1 = open(path_inc_ncl_3_1,'w')
        inc_ncl_list_3_1 = generator.inc_ncl_list[:]
        inc_ncl_list_3_1 += ['SAXIS = 1.0 0.0 0.0\n']

        for j in range(len(inc_ncl_list_3_1)):
            inc_ncl_3_1.write(str(inc_ncl_list_3_1[j]))

        inc_ncl_3_1.close()


    # INCAR_3_2 m=0,0,1

        path_inc_ncl_3_2 = 'INCAR_3_2'
        inc_ncl_3_2 = open(path_inc_ncl_3_2,'w')
        inc_ncl_list_3_2 = generator.inc_ncl_list[:]
        inc_ncl_list_3_2 += ['SAXIS = 0.0 0.0 1.0\n']

        for j in range(len(inc_ncl_list_3_2)):
            inc_ncl_3_2.write(str(inc_ncl_list_3_2[j]))

        inc_ncl_3_2.close()



     # INCAR_4_1 m=0,1,0

        path_inc_ncl_4_1 = 'INCAR_4_1'
        inc_ncl_4_1 = open(path_inc_ncl_4_1,'w')
        inc_ncl_list_4_1 = generator.inc_ncl_list[:]
        inc_ncl_list_4_1 += ['SAXIS = 0 1.0 0\n']

        for j in range(len(inc_ncl_list_4_1)):
            inc_ncl_4_1.write(str(inc_ncl_list_4_1[j]))

        inc_ncl_4_1.close()


    # INCAR_4_2 m=0,0,1

        path_inc_ncl_4_2 = 'INCAR_4_2'
        inc_ncl_4_2 = open(path_inc_ncl_4_2,'w')
        inc_ncl_list_4_2 = generator.inc_ncl_list[:]
        inc_ncl_list_4_2 += ['SAXIS = 0.0 0.0 1.0\n']

        for j in range(len(inc_ncl_list_4_2)):
            inc_ncl_4_2.write(str(inc_ncl_list_4_2[j]))

        inc_ncl_4_2.close()



     # INCAR_5_1 m=1,0,0

        path_inc_ncl_5_1 = 'INCAR_5_1'
        inc_ncl_5_1 = open(path_inc_ncl_5_1,'w')
        inc_ncl_list_5_1 = generator.inc_ncl_list[:]
        inc_ncl_list_5_1 += ['SAXIS = 1.0 0.0 0.0\n']

        for j in range(len(inc_ncl_list_5_1)):
            inc_ncl_5_1.write(str(inc_ncl_list_5_1[j]))

        inc_ncl_5_1.close()


    # INCAR_5_2 m=0,0,1

        path_inc_ncl_5_2 = 'INCAR_5_2'
        inc_ncl_5_2 = open(path_inc_ncl_5_2,'w')
        inc_ncl_list_5_2 = generator.inc_ncl_list[:]
        inc_ncl_list_5_2 += ['SAXIS = 0 0.0 1.0\n']

        for j in range(len(inc_ncl_list_5_2)):
            inc_ncl_5_2.write(str(inc_ncl_list_5_2[j]))

        inc_ncl_5_2.close()


     # INCAR_6_1 m=0,1,0

        path_inc_ncl_6_1 = 'INCAR_6_1'
        inc_ncl_6_1 = open(path_inc_ncl_6_1,'w')
        inc_ncl_list_6_1 = generator.inc_ncl_list[:]
        inc_ncl_list_6_1 += ['SAXIS = 0.0 1.0 0\n']

        for j in range(len(inc_ncl_list_6_1)):
            inc_ncl_6_1.write(str(inc_ncl_list_6_1[j]))

        inc_ncl_6_1.close()


    # INCAR_6_2 m=0,0,1

        path_inc_ncl_6_2 = 'INCAR_6_2'
        inc_ncl_6_2 = open(path_inc_ncl_6_2,'w')
        inc_ncl_list_6_2 = generator.inc_ncl_list[:]
        inc_ncl_list_6_2 += ['SAXIS = 0.0 0.0 1.0\n']

        for j in range(len(inc_ncl_list_6_2)):
            inc_ncl_6_2.write(str(inc_ncl_list_6_2[j]))

        inc_ncl_6_2.close()





     # INCAR_7_1 m=1,1,0

        path_inc_ncl_7_1 = 'INCAR_7_1'
        inc_ncl_7_1 = open(path_inc_ncl_7_1,'w')
        inc_ncl_list_7_1 = generator.inc_ncl_list[:]
        inc_ncl_list_7_1 += ['SAXIS = 1.0 1.0 0\n']

        for j in range(len(inc_ncl_list_7_1)):
            inc_ncl_7_1.write(str(inc_ncl_list_7_1[j]))

        inc_ncl_7_1.close()


    # INCAR_7_2 m=0,0,1

        path_inc_ncl_7_2 = 'INCAR_7_2'
        inc_ncl_7_2 = open(path_inc_ncl_7_2,'w')
        inc_ncl_list_7_2 = generator.inc_ncl_list[:]
        inc_ncl_list_7_2 += ['SAXIS = 0.0 0.0 1.0\n']

        for j in range(len(inc_ncl_list_7_2)):
            inc_ncl_7_2.write(str(inc_ncl_list_7_2[j]))

        inc_ncl_7_2.close()



     # INCAR_8_1 m=1,0,1

        path_inc_ncl_8_1 = 'INCAR_8_1'
        inc_ncl_8_1 = open(path_inc_ncl_8_1,'w')
        inc_ncl_list_8_1 = generator.inc_ncl_list[:]
        inc_ncl_list_8_1 += ['SAXIS = 1.0 0.0 1.0\n']

        for j in range(len(inc_ncl_list_8_1)):
            inc_ncl_8_1.write(str(inc_ncl_list_8_1[j]))

        inc_ncl_8_1.close()


    # INCAR_8_2 m=0,0,1

        path_inc_ncl_8_2 = 'INCAR_8_2'
        inc_ncl_8_2 = open(path_inc_ncl_8_2,'w')
        inc_ncl_list_8_2 = generator.inc_ncl_list[:]
        inc_ncl_list_8_2 += ['SAXIS = 0.0 0.0 1.0\n']

        for j in range(len(inc_ncl_list_8_2)):
            inc_ncl_8_2.write(str(inc_ncl_list_8_2[j]))

        inc_ncl_8_2.close()


     # INCAR_9_1 m=0,1,1

        path_inc_ncl_9_1 = 'INCAR_9_1'
        inc_ncl_9_1 = open(path_inc_ncl_9_1,'w')
        inc_ncl_list_9_1 = generator.inc_ncl_list[:]
        inc_ncl_list_9_1 += ['SAXIS = 0.0 1.0 1.0\n']

        for j in range(len(inc_ncl_list_9_1)):
            inc_ncl_9_1.write(str(inc_ncl_list_9_1[j]))

        inc_ncl_9_1.close()


    # INCAR_9_2 m=0,0,1

        path_inc_ncl_9_2 = 'INCAR_9_2'
        inc_ncl_9_2 = open(path_inc_ncl_9_2,'w')
        inc_ncl_list_9_2 = generator.inc_ncl_list[:]
        inc_ncl_list_9_2 += ['SAXIS = 0.0 0.0 1.0\n']

        for j in range(len(inc_ncl_list_9_2)):
            inc_ncl_9_2.write(str(inc_ncl_list_9_2[j]))

        inc_ncl_9_2.close()






    # Derivation of magnetostriction coefficients:

    if args.der == True:

        for j in range(1,10):

            for k in range(1,3):

                path_dat = "ene_" + str(j) + "_" + str(k) + ".dat"
                dat = open(path_dat,'w')


                for i in range(int(args.ndist[0])):

                    pos_name = "POSCAR_" + str(j) + "_" + str(i+1)

                    struct = Structure.from_file(pos_name)

                    latt = struct.lattice.matrix

                    if j == 1:
                        var1 = latt[0][0]
                    elif j == 2:
                        var1 = latt[0][0]
                    elif j == 3:
                        var1 = latt[1][1]
                    elif j == 4:
                        var1 = latt[1][1]
                    elif j == 5:
                        var1 = latt[2][2]
                    elif j == 6:
                        var1 = latt[2][2]
                    elif j == 7:
                        var1 = math.sqrt((latt[0][0]+latt[1][0])**2+(latt[0][1]+latt[1][1])**2+(latt[0][2]+latt[1][2])**2)
                    elif j == 8:
                        var1 = math.sqrt((latt[0][0]+latt[2][0])**2+(latt[0][1]+latt[2][1])**2+(latt[0][2]+latt[2][2])**2)
                    elif j == 9:
                        var1 = math.sqrt((latt[1][0]+latt[2][0])**2+(latt[1][1]+latt[2][1])**2+(latt[1][2]+latt[2][2])**2)


                    path_osz = "OSZICAR_" + str(j) + "_" + str(i+1) + "_" + str(k)
                    osz = open(path_osz,'r')
                    ene0 = osz.readlines()
                    ene1 = ene0[len(ene0)-2]
                    ene2 = ene1[11:32]

                    osz.close()

                    dat.write(repr(var1))
                    dat.write('  ')
                    dat.write(str(ene2))
                    dat.write('\n')


                dat.close()


       # fitting and plot


        def K(x,a,b,c):
            return a*x*x+b*x+c


        print("")
        print("Fit of quadratic function f(x)=a*x\u00B2+b*x+c to energy vs distortion data")


        lambda_ortho = []


        list_spin = ['1,0,0','0,0,1','0,1,0','0,0,1','1,0,0','0,0,1','0,1,0','0,0,1','1,0,0','0,0,1','0,1,0','0,0,1','1,1,0','0,0,1','1,0,1','0,0,1','0,1,1','0,0,1']
        list_dist = ['1,0,0','1,0,0','0,1,0','0,1,0','0,0,1','0,0,1','1,1,0','1,0,1','0,1,1']

        # AELAS and IEEE lattice convention: c<a<b


        aa0 = SpacegroupAnalyzer(structure0,symprec=sym1, angle_tolerance=sym2)
        structure1 = aa0.get_conventional_standard_structure(international_monoclinic=True)

        bb0 = ConventionalCellTransformation(symprec=sym1, angle_tolerance=sym2, international_monoclinic=True)
        structure2 = bb0.apply_transformation(structure1)
           
        
        latt0 = structure2.lattice.matrix
        coordsnew = np.zeros((len(structure2.species), 3))

        for i in range(len(structure2.species)):
            coordsnew[i][0] = float(structure2.frac_coords[i][1])
            coordsnew[i][1] = float(structure2.frac_coords[i][2])
            coordsnew[i][2] = float(structure2.frac_coords[i][0])


        lattice = Lattice.from_parameters(a=latt0[1][1], b=latt0[2][2], c=latt0[0][0], alpha=90, beta=90, gamma=90)
        structure2b = Structure(lattice, structure2.species, coordsnew)
        
        latt_par = structure2b.lattice.matrix
        
        latt_a = latt_par[0][0]
        latt_b = latt_par[1][1]
        latt_c = latt_par[2][2]
        
        
        for i in range(1,10):


            ene_dat1 = "ene_" + str(i) + "_1.dat"
            ene_dat2 = "ene_" + str(i) + "_2.dat"
            spin1 = str(list_spin[2*i-2])
            spin2 = str(list_spin[2*i-1])
            dist = str(list_dist[i-1])
            fig1 = 'fit_ene_' + str(i) + '_1.png'
            fig2 = 'fit_ene_' + str(i) + '_2.png'


            print("")
            print("-------------------------")
            print("Calculation of \u03BB", i,":")
            print("-------------------------")
            print(" ")
            print('Lattice distorsion along [', dist ,'] direction')
            print("")



            f = open(ene_dat1,'r')
            l = f.readlines()
            f.close

            x = []
            y = []
            for j in l:
                x.append(float(j.split()[0]))
                y.append(float(j.split()[1]))

            x = np.array(x)
            y = np.array(y)

            params = curve_fit(K, x, y)


            print('Fitting parameters for spin parallel to [', spin1 ,'] data from file ',ene_dat1,')')
            print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])

            r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
            print("R-squared =", r_squared)
            print("")

            if r_squared < 0.98:
                print("WARNING!! R-squared is lower than 0.98. Check figure ", fig1)
                print("")

            l1 = -params[0][1] / (2.0 * params[0][0])

            print("X minimum = -b/(2*a) =", l1)
            print("")


            plt.plot(x, y, 'bo', label=ene_dat1 )
            popt, pcov = curve_fit(K, x, y)
            t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
            plt.plot(t, K(t, *popt), 'r--', label='fit')
            plt.ylabel('Energy (eV)')
            plt.legend()
            label = "Lattice distorsion along [" + str(dist) + "] direction (Å)"
            tit = 'Calculation of \u03BB' + str(i) + ', spin = [' + str(spin1) + '] '
            plt.xlabel(label)
            plt.title(tit)
            plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
            plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)


            plt.savefig(fig1)
            plt.close()




            f = open(ene_dat2,'r')
            l = f.readlines()
            f.close

            x = []
            y = []
            for j in l:
                x.append(float(j.split()[0]))
                y.append(float(j.split()[1]))

            x = np.array(x)
            y = np.array(y)

            params = curve_fit(K, x, y)

            print('Fitting parameters for spin parallel to [', spin2 ,'] data from file ',ene_dat2,')')
            print("a =", params[0][0], ", b =", params[0][1], ", c =", params[0][2])

            r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
            print("R-squared =", r_squared)
            print("")

            if r_squared < 0.98:
                print("WARNING!! R-squared is lower than 0.98. Check figure ", fig2)
                print("")

            l2 = -params[0][1] / (2.0 * params[0][0])

            print("X minimum = -b/(2*a) =", l2)
            print("")


            plt.plot(x, y, 'bo', label=ene_dat2 )
            popt, pcov = curve_fit(K, x, y)
            t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
            plt.plot(t, K(t, *popt), 'r--', label='fit')
            plt.ylabel('Energy (eV)')
            plt.legend()
            label = "Lattice distorsion along [" + str(dist) + "] direction (Å)"
            tit = "Calculation of \u03BB" + str(i) + ", spin = [" + str(spin2) + "] "
            plt.xlabel(label)
            plt.title(tit)
            plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
            plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)


            plt.savefig(fig2)
            plt.close()


            #make figure dE_X.png

            f1 = open(ene_dat1,'r')
            f2 = open(ene_dat2,'r')
            s1 = f1.readlines()
            s2 = f2.readlines()
            f1.close
            f2.close

            x = []
            y = []
            y2 = []
            for j in s1:
                x.append(float(j.split()[0]))
                y.append(float(j.split()[1]))

            for j in s2:
                y2.append(float(j.split()[1]))

            x = np.array(x)
            y = np.array(y)
            y2 = np.array(y2)


            plt.plot(x, (y2-y)*1e6, 'o-')

            ylabel ='E[' + str(spin2) + '] - E['+ str(spin1) + '] (\u03BCeV)'
            plt.ylabel(ylabel)
            label = "Lattice distorsion along [" + str(dist) + "] direction (Å)"
            tit = "Calculation of \u03BB" + str(i)
            plt.xlabel(label)
            plt.title(tit)
            plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
            plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)

            fig3 = 'dE_' + str(i) + '.png'
            plt.savefig(fig3)
            plt.close()



            if i == 7:
                
                lmb_1 = lambda_ortho[0]
                
                lmb_2 = lambda_ortho[1]
                
                lmb_3 = lambda_ortho[2]
                
                lmb_4 = lambda_ortho[3]
                
                lmb_7_1 = ((latt_a**2+latt_b**2)*(l1 -l2))/(latt_a*latt_b*(l1+l2))
                
                lmb_7_2 =((latt_a-latt_b)*(latt_a*(lmb_1+lmb_2)-latt_b*(lmb_3+lmb_4)))/(4*latt_a*latt_b)
                
                lmb_7 = lmb_7_1 - lmb_7_2
                
                lambda_ortho += [lmb_7]
            
            elif i == 8:
                
                lmb_1 = lambda_ortho[0]
                
                lmb_5 = lambda_ortho[4]
                
                lmb_8_1 = ((latt_a**2+latt_c**2)*(l1 -l2))/(latt_a*latt_c*(l1+l2))
                
                lmb_8_2 =((latt_a-latt_c)*(latt_a*lmb_1-latt_c*lmb_5))/(4*latt_a*latt_c)
                
                lmb_8 = lmb_8_1 - lmb_8_2
                
                lambda_ortho += [lmb_8]
                   
                
            elif i == 9:
                
                lmb_4 = lambda_ortho[3]
                
                lmb_6 = lambda_ortho[5]
                
                lmb_9_1 = ((latt_b**2+latt_c**2)*(l1 -l2))/(latt_b*latt_c*(l1+l2))
                
                lmb_9_2 =((latt_b-latt_c)*(latt_b*lmb_4-latt_c*lmb_6))/(4*latt_b*latt_c)
                
                lmb_9 = lmb_9_1 - lmb_9_2
                
                lambda_ortho += [lmb_9]
                
                
                
            else:
                
                lambda_ortho += [2.0*((l1 -l2)/(l1+l2))] 


        print(" ")
        print("----------------------------------------------")
        print("Anisotropic magnetostriction coefficients:")
        print("----------------------------------------------")
        print(" ")
        print("Using the convention in reference W.P. Mason, Phys. Rev. 96, 302 (1954):")
        print(" ")

        for i in range(len(lambda_ortho)):

            print("\u03BB",i+1," =", lambda_ortho[i]*1e6,u'x 10\u207B\u2076')
            print(" ")
            
        

        if args.delas == True:

                print(" ")
                print(" ")
                print("----------------------------------------------")
                print("Calculation of magnetoelastic constants:")
                print("----------------------------------------------")
                print(" ")
                print("Reading the elastic tensor file =", str(args.elas[0]))
                print(" ")




                elasdat = open(args.elas[0],'r')
                elasline = elasdat.readlines()
                elasline0 = elasline[2]
                elasline1 = elasline[3]
                elasline2 = elasline[4]
                elasline3 = elasline[5]
                elasline4 = elasline[6]
                elasline5 = elasline[7]
                c11 = float(elasline0[0:8])
                c12 = float(elasline0[8:16])
                c13 = float(elasline0[16:24])
                c23 = float(elasline1[16:24])
                c22 = float(elasline1[8:16])
                c33 = float(elasline2[16:24])
                c44 = float(elasline3[24:32])
                c55 = float(elasline4[32:40])
                c66 = float(elasline5[40:48])

                elasdat.close()


                b1 = -c11*lambda_ortho[0]-c12*lambda_ortho[2]-c13*lambda_ortho[4]

                b2 = -c11*lambda_ortho[1]-c12*lambda_ortho[3]-c13*lambda_ortho[5]

                b3 = -c12*lambda_ortho[0]-c22*lambda_ortho[2]-c23*lambda_ortho[4]

                b4 = -c12*lambda_ortho[1]-c22*lambda_ortho[3]-c23*lambda_ortho[5]

                b5 = -c13*lambda_ortho[0]-c23*lambda_ortho[2]-c33*lambda_ortho[4]

                b6 = -c13*lambda_ortho[1]-c23*lambda_ortho[3]-c33*lambda_ortho[5]

                b7 = c66*(lambda_ortho[0]+lambda_ortho[1]+lambda_ortho[2]+lambda_ortho[3]-4*lambda_ortho[6])

                b8 = c55*(lambda_ortho[0]+lambda_ortho[4]-4*lambda_ortho[7])

                b9 = c44*(lambda_ortho[3]+lambda_ortho[5]-4*lambda_ortho[8])


                print("c11 =", str(c11), 'GPa')
                print(" ")
                print("c12 =", str(c12), 'GPa')
                print(" ")
                print("c13 =", str(c13), 'GPa')
                print(" ")
                print("c23 =", str(c23), 'GPa')
                print(" ")
                print("c22 =", str(c22), 'GPa')
                print(" ")
                print("c33 =", str(c33), 'GPa')
                print(" ")
                print("c44 =", str(c44), 'GPa')
                print(" ")
                print("c55 =", str(c55), 'GPa')
                print(" ")
                print("c66 =", str(c66), 'GPa')
                print(" ")
                print("Warning: If these elastic constants are not the same as in the input elastic tensor file", str(args.elas[0]),", then check that the format of the elastic tensor is exactly the same as in the standard output file ELADAT generated by ELAS code (see Example folder)")
                print(" ")
                print(" ")
                print("Magnetoelastic constants:")
                print(" ")
                print("b1 =", str(b1), 'GPa')
                print(" ")
                print("b2 =", str(b2), 'GPa')
                print(" ")
                print("b3 =", str(b3), 'GPa')
                print(" ")
                print("b4 =", str(b4), 'GPa')
                print(" ")
                print("b5 =", str(b5), 'GPa')
                print(" ")
                print("b6 =", str(b6), 'GPa')
                print(" ")
                print("b7 =", str(b7), 'GPa')
                print(" ")
                print("b8 =", str(b8), 'GPa')
                print(" ")
                print("b9 =", str(b9), 'GPa')
                print(" ")
                print("The equation of the magnetoelastic energy can be found in the User Manual")
                print(" ")
