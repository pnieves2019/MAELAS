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

        # Convention: lattice vector a1 along x-axis
        angle = -math.pi*(60.0/180.0)
        dd = DeformStructureTransformation(deformation=((math.cos(angle), math.sin(angle), 0), (-math.sin(angle), math.cos(angle), 0), (0, 0, 1)))
        structure2b = dd.apply_transformation(structure2)


        for i in range(int(self.args.ndist[0])):


            strain1 = - float(self.args.strain[0])+2*(float(self.args.strain[0])/(float(self.args.ndist[0])-1))*i

            print("strain", strain1)


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

######################### OLD deformation v1.0
            
           # pos_name7 = "POSCAR_5_" + str(i+1)

           # structure77 = Poscar(structure6)
           # structure77.write_file(filename = pos_name7,significant_figures=16)
           
           
######################## NEW deformation v2.0


            a1 = 1.0 + strain1
            a2 = 1/math.sqrt(a1)
            a3 = a2
            dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
            structure7 = dd.apply_transformation(structure2b)
            pos_name7 = "POSCAR_5_" + str(i+1)

            structure77 = Poscar(structure7)
            structure77.write_file(filename = pos_name7,significant_figures=16)

#######################################           

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

  def der(self):
  
        structure0 = Structure.from_file(self.args.pos[0])
        nat = len(structure0.species)
        vol = float(structure0.volume)
        
        sym1 = float(self.args.sympre[0])
        sym2 = float(self.args.symang[0])
    
        aa = SpacegroupAnalyzer(structure0,symprec=sym1, angle_tolerance=sym2)
        
        
        for j in range(1,7):

            for k in range(1,3):

                path_dat = "ene_" + str(j) + "_" + str(k) + ".dat"
                dat = open(path_dat,'w')


                for i in range(int(self.args.ndist[0])):

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
                        # OLD length direction v1.0
                        #var1 = math.sqrt((latt[0][0]+latt[2][0])**2+(latt[0][1]+latt[2][1])**2+(latt[0][2]+latt[2][2])**2)
                        # NEW length direction v2.0
                        var1 = latt[0][0]
                        
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
        print("Fit of quadratic function f(x)=A*x\u00B2+B*x+C to energy vs cell length")
        print("")
        print("-------------------------")
        print("Calculation of \u03BB \u03B11,2:")
        print("-------------------------")
        print(" ")
        print('Unit cell length along [1,0,0] direction')
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
        plt.xlabel('Unit cell length along [1,0,0] direction (\u212B)')
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


        lambda_alpha_1_2 = 2.0*((l1 -l2)/(l1+l2))


        plt.plot(x, y, 'bo', label='data in ene_1_2.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Unit cell length along [0,0,1] direction (\u212B)')
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
        label = "Unit cell length along [" + str(dist) + "] direction (\u212B)"
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
        print('Unit cell length along [0,0,1] direction')
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
        
        nn = int(self.args.ndist[0])+1
        if nn % 2 == 0:
            lli = int((nn-2)/2)
            mae001 = y[lli]


        plt.plot(x, y, 'bo', label='data in ene_2_1.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Unit cell length along [0,0,1] direction (\u212B)')
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
            mae100 = y[lli]

        lambda_alpha_2_2 = 2.0*((l1 -l2)/(l1+l2))


        plt.plot(x, y, 'bo', label='data in ene_2_2.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.xlabel('Unit cell length along [0,0,1] direction (\u212B)')
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
        label = "Unit cell length along [" + str(dist) + "] direction (\u212B)"
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
        print('Unit cell length along [1,0,0] direction')
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
        print("A =", params[0][0], ", B =", params[0][1], ", C =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_3_1.png")
            print("")

        l1 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -B/(2*A) =", l1)
        print("")


        plt.plot(x, y, 'bo', label='data in ene_3_1.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Unit cell length along [1,0,0] direction (\u212B)')
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
        print("A =", params[0][0], ", B =", params[0][1], ", C =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_3_2.png")
            print("")

        l2 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -B/(2*A) =", l2)
        print("")

        lambda_gamma_1 = 2.0*((l1 -l2)/(l1+l2))


        plt.plot(x, y, 'bo', label='data in ene_3_2.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.xlabel('Unit cell length along [1,0,0] direction (\u212B)')
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
        label = "Unit cell length along [" + str(dist) + "] direction (\u212B)"
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
        print('Unit cell length along [a,0,c] direction')
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
        print("A =", params[0][0], ", B =", params[0][1], ", C =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_4_1.png")
            print("")

        l1 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -B/(2*A) =", l1)
        print("")


        plt.plot(x, y, 'bo', label='data in ene_4_1.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Unit cell length along [a,0,c] direction (\u212B)')
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
        print("A =", params[0][0], ", B =", params[0][1], ", C =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_4_2.png")
            print("")

        l2 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -B/(2*A) =", l2)
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
        plt.xlabel('Unit cell length along [a,0,c] direction (\u212B)')
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
        dist = 'a,0,c'
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
        label = "Unit cell length along [" + str(dist) + "] direction (\u212B)"
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
        print('Unit cell length along [1,0,0] direction')
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
        print("A =", params[0][0], ", B =", params[0][1], ", C =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_5_1.png")
            print("")

        l1 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -B/(2*A) =", l1)
        print("")


        plt.plot(x, y, 'bo', label='data in ene_5_1.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Unit cell length along [1,0,0] direction (\u212B)')
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
        print("A =", params[0][0], ", B =", params[0][1], ", C =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_5_2.png")
            print("")

        l2 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -B/(2*A) =", l2)
        print("")

        # OLD calculation v1.0
        
        #eta_par = (latt_a*latt_c)/(2.0*(latt_a**2+latt_c**2))
        #lambda_1_2 = 2.0*((l1 -l2)/(l1+l2))*(1.0/eta_par)
        
        # New calculation v2.0
        lambda_1_2 = 2.0*2.0*((l1 -l2)/(l1+l2))

        plt.plot(x, y, 'bo', label='data in ene_5_2.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.xlabel('Unit cell length along [1,0,0] direction (\u212B)')
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
        dist = '1,0,0'
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
        label = "Unit cell length along [" + str(dist) + "] direction (\u212B)"
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
        print('Unit cell length along [a,0,c] direction')
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
        print("A =", params[0][0], ", B =", params[0][1], ", C =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_6_1.png")
            print("")

        l1 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -B/(2*A) =", l1)
        print("")


        plt.plot(x, y, 'bo', label='data in ene_6_1.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Unit cell length along [a,0,c] direction (\u212B)')
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
        print("A =", params[0][0], ", B =", params[0][1], ", C =", params[0][2])

        r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
        print("R-squared =", r_squared)
        print("")

        if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_6_2.png")
            print("")

        l2 = -params[0][1] / (2.0 * params[0][0])

        print("X minimum = -B/(2*A) =", l2)
        print("")

        eta_par = (latt_a*latt_c)/(latt_a**2+latt_c**2)
        
        lambda_2_1 = 2.0*((l1 -l2)/(l1+l2))*(1.0/eta_par)


        plt.plot(x, y, 'bo', label='data in ene_6_2.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.xlabel('Unit cell length along [a,0,c] direction (\u212B)')
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
        dist = 'a,0,c'
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
        print(" ")
        print("Warning!: The results calculated with -mode 1 can have significant errors, especially for non-cubic crystal symmetries. We strongly recommned to verify these results through -mode 2.")
        print(" ")
        print(" ")
        print("/////////////////////")
        print("Spontaneous volume magnetostriction (anisotropic contribution):")
        print("/////////////////////")
        print(" ")
        print("Trigonal crystal with easy axis (ferromagnetic state with equal domains along all equivalent easy directions): ")
        print(" ")
        print("\u03C9_s^ani =", ((4.0/3.0)*lambda_alpha_1_2+(2.0/3.0)*lambda_alpha_2_2)*1e6,u'x 10\u207B\u2076')
        print(" ")
        print("......................... ")
        print(" ")
        print("Trigonal crystal with easy plane (ferromagnetic state with equal domains along all equivalent easy directions): ")
        print(" ")
        print("\u03C9_s^ani =", (-(2.0/3.0)*lambda_alpha_1_2-(1.0/3.0)*lambda_alpha_2_2)*1e6,u'x 10\u207B\u2076')
        print(" ")
        print("/////////////////////")
        print("Polycrystal (uniform stress approximation):")
        print("/////////////////////")
        print(" ")
        print(" ")
        print("Trigonal crystal with easy axis (reference demagnetized state with equal domains along all easy directions): ")
        print(" ")
        print("\u03BE =", 1e6*((4.0/15.0)*lambda_alpha_1_2+(1.0/15.0)*lambda_alpha_2_2-(2.0/15.0)*lambda_gamma_1-(1.0/15.0)*lambda_gamma_2-(1.0/3.0)*(2.0*lambda_alpha_1_2+lambda_alpha_2_2)),u'x 10\u207B\u2076')
        print(" ")
        print("\u03B7 =", 1e6*(-(2.0/15.0)*lambda_alpha_1_2+(2.0/15.0)*lambda_alpha_2_2+(2.0/5.0)*lambda_gamma_1+(1.0/5.0)*lambda_gamma_2),u'x 10\u207B\u2076')
        print(" ")
        print("......................... ")
        print(" ")
        print("Trigonal crystal with easy plane (reference demagnetized state with equal domains along all easy directions): ")
        print(" ")
        print("\u03BE =", 1e6*((4.0/15.0)*lambda_alpha_1_2+(1.0/15.0)*lambda_alpha_2_2-(2.0/15.0)*lambda_gamma_1-(1.0/15.0)*lambda_gamma_2),u'x 10\u207B\u2076')
        print(" ")
        print("\u03B7 =", 1e6*(-(2.0/15.0)*lambda_alpha_1_2+(2.0/15.0)*lambda_alpha_2_2+(2.0/5.0)*lambda_gamma_1+(1.0/5.0)*lambda_gamma_2),u'x 10\u207B\u2076')
        print(" ")
        print(" ") 
        
        
        
        if nn % 2 == 0:         
                print("----------------------------------------------")
                print("Magnetocrystalline anisotropy energy:")
                print("----------------------------------------------")
                print(" ")
                print("These energies correspond to the central points in the data files ene_2_1.dat and ene_2_2.dat:")
                print(" ")
                print("E(0,0,1) = ",mae001," eV")
                print(" ")
                print("E(1,0,0) = ",mae100," eV")
                print(" ")
                print("E(1,0,0) - E(0,0,1) = ",(mae100 - mae001)*1e6,u'x 10\u207B\u2076 eV')
                print(" ")
                print("[E(1,0,0) - E(0,0,1)]/Natom = ",((mae100 - mae001)/nat)*1e6,u'x 10\u207B\u2076 eV/atom')


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


                #OLD b(lmb) v1.0 

                #b21 = -(c11+c12)*lambda_alpha_1_2-c13*lambda_alpha_2_2

                #b22 = -2*c13*lambda_alpha_1_2-c33*lambda_alpha_2_2

                #b3 = c14*lambda_2_1+0.5*(-c11+c12)*lambda_gamma_1

                #b4 = -c14*lambda_1_2+c44*lambda_gamma_2

                #b14 =  c44*lambda_2_1-c14*lambda_gamma_1

                #b34 = 0.5*(-c11+c12)*lambda_1_2+c14*lambda_gamma_2
                
                #NEW b(lmb) v2.0 
                          
                b21 = -(c11+c12)*lambda_alpha_1_2-c13*lambda_alpha_2_2

                b22 = -2.0*c13*lambda_alpha_1_2-c33*lambda_alpha_2_2
                
                b3 = -c14*lambda_2_1+(c12-c11)*lambda_gamma_1
                
                b4 = -c14*lambda_1_2-c44*lambda_gamma_2
                
                b14 = -c44*lambda_2_1-2.0*c14*lambda_gamma_1
                
                b34 = 0.5*(c12-c11)*lambda_1_2-c14*lambda_gamma_2



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
                print("Warning: If these elastic constants are not the same as in the input elastic tensor file", str(self.args.elas[0]),", then check that the format of the elastic tensor is exactly the same as in the standard output file ELADAT generated by AELAS code (see Example folder)")
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


        
        
        