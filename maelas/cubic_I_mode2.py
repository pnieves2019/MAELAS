#!/bin/bash

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import stat

import maelas.parser   as parser
import maelas.generate as generate

from pymatgen.core import Lattice, Structure
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
        generator.poscar2()
        generator.incar()
    
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

        #b_1


            a1 = 1.0 + strain1
            a2 = 1.0
            a3 = 1.0
            dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
            structure3 = dd.apply_transformation(structure2)
            pos_name = "POSCAR_1_" + str(i+1)

            structure33 = Poscar(structure3)
            structure33.write_file(filename = pos_name,significant_figures=16)



        #b_2

            

            a12 = strain1
            a13 = 0.0
            a21 = strain1
            a23 = 0.0
            a31 = 0.0
            a32 = 0.0

            a11 = 1.0
            a22 = 1.0
            a33 = 1.0

            cc = DeformStructureTransformation(deformation=((a11, a12, a13), (a21, a22, a23), (a31, a32, a33)))
            structure4 = cc.apply_transformation(structure2)
            pos_name2 = "POSCAR_2_" + str(i+1)

            structure44 = Poscar(structure4)
            structure44.write_file(filename = pos_name2,significant_figures=16)



    # INCAR_1_1 m=1,0,0

        path_inc_ncl_1_1 = 'INCAR_1_1'
        inc_ncl_1_1 = open(path_inc_ncl_1_1,'w')
        inc_ncl_list_1_1 = generator.inc_ncl_list[:]
        inc_ncl_list_1_1 += ['SAXIS = 1.0 0 0.0\n']

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

    # INCAR_2_1 m=1,1,0

        path_inc_ncl_2_1 = 'INCAR_2_1'
        inc_ncl_2_1 = open(path_inc_ncl_2_1,'w')
        inc_ncl_list_2_1 = generator.inc_ncl_list[:]
        inc_ncl_list_2_1 += ['SAXIS = 1.0 1.0 0.0\n']

        for j in range(len(inc_ncl_list_2_1)):
            inc_ncl_2_1.write(str(inc_ncl_list_2_1[j]))

        inc_ncl_2_1.close()


    # INCAR_2_2 m=-1,1,0

        path_inc_ncl_2_2 = 'INCAR_2_2'
        inc_ncl_2_2 = open(path_inc_ncl_2_2,'w')
        inc_ncl_list_2_2 = generator.inc_ncl_list[:]
        inc_ncl_list_2_2 += ['SAXIS = -1.0 1.0 0.0\n']

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
                    
                    strain1 = - float(self.args.strain[0])+2*(float(self.args.strain[0])/(float(self.args.ndist[0])-1))*i

                    print("strain", strain1)


                    path_osz = "OSZICAR_" + str(j) + "_" + str(i+1) + "_" + str(k)
                    osz = open(path_osz,'r')
                    ene0 = osz.readlines()
                    ene1 = ene0[len(ene0)-2]
                    ene2 = ene1[11:32]

                    osz.close()

                    dat.write(str(strain1))
                    dat.write('  ')
                    dat.write(str(ene2))
                    dat.write('\n')


                dat.close()


       # fitting and plot


        def K(x,a,b):
            return a*x+b


        print("")
        print("Fit of linear function f(x)=A*x+B to energy vs strain")
        print("")
        print("-------------------------")
        print("Calculation of b_1:")
        print("-------------------------")
        print(" ")
        print('Deformation x-direction (\u03B5xx), spin1=[1,0,0] and spin2=[1,1,0]')
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
        
        
        f2 = open('ene_1_2.dat','r')
        l2 = f2.readlines()
        f2.close

        x2 = []
        y2 = []
        for i in l2:
            x2.append(float(i.split()[0]))
            y2.append(float(i.split()[1]))

        x2 = np.array(x2)
        y2 = np.array(y2)
        

        params = curve_fit(K, x, y-y2)

        print("Fitting parameters:")
        print("A =", params[0][0], ", B =", params[0][1])

        r_squared = r2_score(y-y2, K(x,params[0][0],params[0][1]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure dE_1.png")
            print("")

        b1 = 2.0*params[0][0]


        nn = int(self.args.ndist[0])+1
        if nn % 2 == 0:
            lli = int((nn-2)/2)
            mae100 = y[lli]
            mae110 = y2[lli]




        #make figure dE_1.png

        fig = 'dE_1.png'
        spin1 = '1,0,0'
        spin2 = '1,1,0'
        dist = '1,0,0'
        tit = "Calculation of b1"

        plt.plot(x, (y-y2)*1e6, 'o', label='data')
        plt.plot(x, K(x, params[0][0],params[0][1])*1e6, 'r-', label='fit')
        plt.legend()
        ax = plt.gca()
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))

        ylabel ='E[' + str(spin1) + '] - E['+ str(spin2) + '] (\u03BCeV)'
        plt.ylabel(ylabel)
        label = "Strain \u03B5xx"
        plt.xlabel(label)
        plt.title(tit)
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig(fig)
        plt.close()




        print("")
        print("Fit of linear function f(x)=A*x+B to energy vs strain")
        print("")
        print("-------------------------")
        print("Calculation of b_2:")
        print("-------------------------")
        print(" ")
        print('Deformation xy-direction  (\u03B5xy), spin1=[1,1,0] and spin2=[-1,1,0]')
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
        
        
        f2 = open('ene_2_2.dat','r')
        l2 = f2.readlines()
        f2.close

        x2 = []
        y2 = []
        for i in l2:
            x2.append(float(i.split()[0]))
            y2.append(float(i.split()[1]))

        x2 = np.array(x2)
        y2 = np.array(y2)
        

        params = curve_fit(K, x, y-y2)

        print("Fitting parameters:")
        print("A =", params[0][0], ", B =", params[0][1])

        r_squared = r2_score(y-y2, K(x,params[0][0],params[0][1]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure dE_2.png")
            print("")

        b2 = params[0][0]/2.0






        #make figure dE_2.png

        fig = 'dE_2.png'
        spin1 = '1,1,0'
        spin2 = '-1,1,0'
        dist = '1,1,0'
        tit = "Calculation of b2" 

        plt.plot(x, (y-y2)*1e6, 'o', label='data')
        plt.plot(x, K(x, params[0][0],params[0][1])*1e6, 'r-', label='fit')
        plt.legend()
        ax = plt.gca()
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))

        ylabel ='E[' + str(spin1) + '] - E['+ str(spin2) + '] (\u03BCeV)'
        plt.ylabel(ylabel)
        label = "Strain \u03B5xy"
        plt.xlabel(label)
        plt.title(tit)
        plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
        plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
        plt.savefig(fig)
        plt.close()










        print(" ")
        print("----------------------------------------------")
        print("Anisotropic magnetoelastic constants:")
        print("----------------------------------------------")
        print(" ")
        print(" ")
        print("Using the convention in reference P. Nieves et al., Comput. Phys. Commun. 264, 107964 (2021):")
        print(" ")
        print("b1 =", b1,u' eV')
        print(" ")
        print("b2 =", b2,u' eV')
        print(" ")
        print("b1 =", b1/nat,u' eV/atom')
        print(" ")
        print("b2 =", b2/nat,u' eV/atom')
        print(" ")
        print("b1 =", b1/vol,u' eV/A^3')
        print(" ")
        print("b2 =", b2/vol,u' eV/A^3')
        print(" ")
        print("b1 =", (b1/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
        print(" ")
        print("b2 =", (b2/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
        print(" ")
        print("Using the convention in reference E. D. T. de Lacheisserie, Magnetostriction: Theory and Application of Magnetoelasticity (CRC Press, Boca Raton, FL, 1993).")
        print(" ")
        print("b \u03B3,2 =", (b1/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
        print(" ")
        print("b \u03B5,2 =", (b2/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
        print(" ")

        bani1 = (b1/vol)*1.602176565e-19*1e30*1e-6
        bani2 = (b2/vol)*1.602176565e-19*1e30*1e-6
        
        path_inc_std = 'MAGANI'
        inc_std = open(path_inc_std,'w')
        inc_std.write(str(bani1))
        inc_std.write('\n')
        inc_std.write(str(bani2))
        inc_std.close()

        
        
        if nn % 2 == 0:         
            print("----------------------------------------------")
            print("Magnetocrystalline anisotropy energy:")
            print("----------------------------------------------")
            print(" ")
            print("These energies correspond to the central points in the data files ene_1_1.dat and ene_1_2.dat:")
            print(" ")
            print("E(1,0,0) = ",mae100," eV")
            print(" ")
            print("E(1,1,0) = ",mae110," eV")
            print(" ")
            print("E(1,1,0) - E(1,0,0) = ",(mae110 - mae100)*1e6,u'x 10\u207B\u2076 eV')
            print(" ")
            print("[E(1,1,0) - E(1,0,0)]/Natom = ",((mae110 - mae100)/nat)*1e6,u'x 10\u207B\u2076 eV/atom')
            print(" ")



        if self.args.delas == True:
            print(" ")
            print(" ")
            print("----------------------------------------------")
            print("Calculation of magnetostrictive coefficients:")
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


            bb1 = (b1/vol)*1.602176565e-19*1e30*1e-9  #GPa

            bb2 = (b2/vol)*1.602176565e-19*1e30*1e-9  #GPa
            
            lambda001 = -((2.0*bb1)/3.0)/(c11-c12)
            
            lambda111 = -bb2/(3.0*c44)
            
            lambda_s = (2.0/5.0)*lambda001 + (3.0/5.0)*lambda111


            print("c11 =", str(c11), 'GPa')
            print(" ")
            print("c12 =", str(c12), 'GPa')
            print(" ")
            print("c44 =", str(c44), 'GPa')
            print(" ")
            print("Warning: If these elastic constants are not the same as in the input elastic tensor file", str(self.args.elas[0]),", then check that the format of the elastic tensor is exactly the same as in the standard output file ELADAT generated by AELAS code (see Example folder)")
            print(" ")
            print(" ")
            print("Magnetostrictive coefficients:")
            print(" ")
            print("Using the convention in reference P. Nieves et al., Comput. Phys. Commun. 264, 107964 (2021):")
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
            print("----------------------------------------------")
            print(" ")
            print("Spontaneous volume magnetostriction (anisotropic contribution):")
            print(" ")
            print("Cubic crystal at ferromagnetic state with equal domains along all equivalent easy directions: ")
            print(" ")
            print("\u03C9_s^ani = 0.0")
            print(" ")
            print("----------------------------------------------")
            print(" ")
            print("Polycrystal (uniform stress approximation):")
            print(" ")
            print("\u03BBs =", lambda_s*1e6,u'x 10\u207B\u2076')




        
        