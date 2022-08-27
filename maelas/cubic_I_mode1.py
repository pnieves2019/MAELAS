#!/bin/bash

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import stat

import maelas.parser   as parser
import maelas.generate as generate


from pymatgen import Lattice, Structure
from pymatgen.transformations.standard_transformations import ConventionalCellTransformation,DeformStructureTransformation
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

class run:

  def __init__(self,args):
        self.args = args
        
  
  
  def gen(self):

        
        generator = generate.VASP(self.args)
        structure0 = Structure.from_file(self.args.pos[0])
    
        sym1 = float(self.args.sympre[0])
        sym2 = float(self.args.symang[0])
    
        aa = SpacegroupAnalyzer(structure0,symprec=sym1, angle_tolerance=sym2)
        
        if self.args.noconv == False:
          structure1 = aa.get_conventional_standard_structure(international_monoclinic=True)
          bb = ConventionalCellTransformation(symprec=sym1, angle_tolerance=sym2, international_monoclinic=True)
          structure2 = bb.apply_transformation(structure1)
          
        
        if self.args.noconv == True:
          structure2 = structure0

        
        for i in range(int(self.args.ndist[0])):


            strain1 = - float(self.args.strain[0])+2*(float(self.args.strain[0])/(float(self.args.ndist[0])-1))*i

            print("strain", strain1)


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

  def der(self):
  
        structure0 = Structure.from_file(self.args.pos[0])
        nat = len(structure0.species)
        vol = float(structure0.volume)
        
        for j in range(1,3):

            for k in range(1,3):

                path_dat = "ene_" + str(j) + "_" + str(k) + ".dat"
                dat = open(path_dat,'w')


                for i in range(int(self.args.ndist[0])):

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
        print("Fit of quadratic function f(x)=A*x\u00B2+B*x+C to energy vs cell length")
        print("")
        print("-------------------------")
        print("Calculation of \u03BB001:")
        print("-------------------------")
        print(" ")
        print('Unit cell length along [0,0,1] direction')
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
        print("A =", params[0][0], ", B =", params[0][1], ", C =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_1_1.png")
            print("")

        l1 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -B/(2*A) =", l1)
        print("")


        plt.plot(x, y, 'bo', label='data in ene_1_1.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Unit cell length along [0,0,1] direction (\u212B)')
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
        print("A =", params[0][0], ", B =", params[0][1], ", C =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_1_2.png")
            print("")

        l2 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -B/(2*A) =", l2)
        print("")

        nn = int(self.args.ndist[0])+1
        if nn % 2 == 0:
            lli = int((nn-2)/2)
            mae001 = y[lli]

        lambda001 = 2.0*(2.0/3.0)*((l1 -l2)/(l1+l2))


        plt.plot(x, y, 'bo', label='data in ene_1_2.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Unit cell length along [0,0,1] direction (\u212B)')
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
        label = "Unit cell length along [" + str(dist) + "] direction (\u212B)"
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
        print('Unit cell length along [1,1,1] direction')
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
        print("A =", params[0][0], ", B =", params[0][1], ", C =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_2_1.png")
            print("")

        l1 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -B/(2*A) =", l1)
        print("")
        
        if nn % 2 == 0:
            lli = int((nn-2)/2)
            mae111 = y[lli]


        plt.plot(x, y, 'bo', label='data in ene_2_1.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Unit cell length along [1,1,1] direction (\u212B)')
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
        print("A =", params[0][0], ", B =", params[0][1], ", C =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_2_2.png")
            print("")

        l2 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -B/(2*A) =", l2)
        print("")
        
        if nn % 2 == 0:
            lli = int((nn-2)/2)
            mae101 = y[lli]

        lambda111 = 2.0*(2.0/3.0)*((l1 -l2)/(l1+l2))

        lambda_s = (2.0/5.0)*lambda001 + (3.0/5.0)*lambda111


        plt.plot(x, y, 'bo', label='data in ene_2_2.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.xlabel('Unit cell length along [1,1,1] direction (\u212B)')
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
        label = "Unit cell length along [" + str(dist) + "] direction (\u212B)"
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
        print("Using the convention in reference E. D. T. de Lacheisserie, Magnetostriction: Theory and Application of Magnetoelasticity (CRC Press, Boca Raton, FL, 1993):")
        print(" ")
        print("\u03BB \u03B3,2 =", (3.0/2.0)*lambda001*1e6,u'x 10\u207B\u2076')
        print(" ")
        print("\u03BB \u03B5,2 =", (3.0/2.0)*lambda111*1e6,u'x 10\u207B\u2076')
        print(" ")
        print("Polycrystal (uniform stress approximation):")
        print(" ")
        print("\u03BBs =", lambda_s*1e6,u'x 10\u207B\u2076')
        
        
        if nn % 2 == 0:         
            print("----------------------------------------------")
            print("Magnetocrystalline anisotropy energy:")
            print("----------------------------------------------")
            print(" ")
            print("These energies correspond to the central points in the data files ene_1_1.dat, ene_2_1.dat, and ene_2_2.dat:")
            print(" ")
            print("E(0,0,1) = ",mae001," eV")
            print(" ")
            print("E(1,1,1) = ",mae111," eV")
            print(" ")
            print("E(1,0,-1) = ",mae101," eV")
            print(" ")
            print("E(1,1,1) - E(0,0,1) = ",(mae111 - mae001)*1e6,u'x 10\u207B\u2076 eV')
            print(" ")
            print("[E(1,1,1) - E(0,0,1)]/Natom = ",((mae111 - mae001)/nat)*1e6,u'x 10\u207B\u2076 eV/atom')
            print(" ")
            print("E(1,0,-1) - E(0,0,1) = ",(mae101 - mae001)*1e6,u'x 10\u207B\u2076 eV')
            print(" ")
            print("[E(1,0,-1) - E(0,0,1)]/Natom = ",((mae101 - mae001)/nat)*1e6,u'x 10\u207B\u2076 eV/atom')
            print(" ")


        if self.args.delas == True:
            print(" ")
            print(" ")
            print("----------------------------------------------")
            print("Calculation of magnetoelastic constants:")
            print("----------------------------------------------")
            print(" ")
            print("Reading the elastic tensor file =", str(self.args.elas[0]))
            print(" ")




            elasdat = open(self.args.elas[0],'r')
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
            print("Warning: If these elastic constants are not the same as in the input elastic tensor file", str(self.args.elas[0]),", then check that the format of the elastic tensor is exactly the same as in the standard output file ELADAT generated by AELAS code (see Example folder)")
            print(" ")
            print(" ")
            print("Magnetoelastic constants:")
            print(" ")
            print("Using the convention in reference P. Nieves et al., Comput. Phys. Commun. 264, 107964 (2021):")
            print(" ")
            print("b1 =", str(b1), 'GPa')
            print(" ")
            print("b2 =", str(b2), 'GPa')
            print(" ")
            print("Using the convention in reference E. D. T. de Lacheisserie, Magnetostriction: Theory and Application of Magnetoelasticity (CRC Press, Boca Raton, FL, 1993):")
            print(" ")
            print("b \u03B3,2 =", str(b1), 'GPa')
            print(" ")
            print("b \u03B5,2 =", str(b2), 'GPa')
            print(" ")




  
      