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
          
        
        if int(self.args.sg0[0]) == 0:
            sg = aa.get_space_group_number()
        elif 0 < int(self.args.sg0[0]) <= 230:
            sg = int(self.args.sg0[0])   
        else:
            print("Space group number must be in the range 1-230")
            exit(-1)

        if 177 <= sg <= 194:
            # Convention: lattice vector a1 along x-axis
            angle = -math.pi*(60.0/180.0)
            dd = DeformStructureTransformation(deformation=((math.cos(angle), math.sin(angle), 0), (-math.sin(angle), math.cos(angle), 0), (0, 0, 1)))
            structure2b = dd.apply_transformation(structure2)
        else:
            structure2b = structure2



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

  def der(self):
  
        generator = generate.VASP(self.args)
        structure0 = Structure.from_file(self.args.pos[0])
        generator.poscar2()
    
        sym1 = float(self.args.sympre[0])
        sym2 = float(self.args.symang[0])
    
        aa = SpacegroupAnalyzer(structure0,symprec=sym1, angle_tolerance=sym2)
        
        if self.args.noconv == False:
          structure1 = aa.get_conventional_standard_structure(international_monoclinic=True)
          bb = ConventionalCellTransformation(symprec=sym1, angle_tolerance=sym2, international_monoclinic=True)
          structure2 = bb.apply_transformation(structure1)
          
        
        if self.args.noconv == True:
          structure2 = structure0
        
        nat = len(structure2.species)
        vol = float(structure2.volume)
        
        
        if int(self.args.sg0[0]) == 0:
            sg = aa.get_space_group_number()
        elif 0 < int(self.args.sg0[0]) <= 230:
            sg = int(self.args.sg0[0])   
        else:
            print("Space group number must be in the range 1-230")
            exit(-1)
        
        
        if (89 <= sg <= 142):
            nmax = 6
        else:
            nmax = 5


        for j in range(1,nmax):

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
        print("Fit of quadratic function f(x)=A*x\u00B2+B*x+C to energy vs cell length")

        print(" ")
        print("-------------------------")
        print('Calculation of \u03BB 1\u03B1,2:')
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


        print("Fitting parameters for spin parallel to 111 (data from file ene_1_1.dat):")
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

        lambda_alpha_1_2 = 2.0*3.0*((l1 -l2)/(l1+l2))


        plt.plot(x, y, 'bo', label='data in ene_1_2.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.xlabel('Unit cell length along [1,0,0] direction (\u212B)')
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
        label = "Unit cell length along [" + str(dist) + "] direction (\u212B)"
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
        label = "Unit cell length along [" + str(dist) + "] direction (\u212B)"
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

        lambda_gamma_2 = 2.0*((l1 -l2)/(l1+l2))


        plt.plot(x, y, 'bo', label='data in ene_3_2.dat')
        popt, pcov = curve_fit(K, x, y)
        t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
        plt.plot(t, K(t, *popt), 'r--', label='fit')
        plt.xlabel('Unit cell length along [1,0,0] direction (\u212B)')
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
        label = "Unit cell length along [" + str(dist) + "] direction (\u212B)"
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
        plt.xlabel('Unit cell length along [a,0,c] direction (\u212B)')
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
        dist = 'a,0,c'
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
        label = "Unit cell length along [" + str(dist) + "] direction (\u212B)"
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
            print('Unit cell length along [1,1,0] direction')
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
            plt.xlabel('Unit cell length along [1,1,0] direction (\u212B)')
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

            lambda_delta = 2.0*((l1 -l2)/(l1+l2))


            plt.plot(x, y, 'bo', label='data in ene_5_2.dat')
            popt, pcov = curve_fit(K, x, y)
            t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
            plt.plot(t, K(t, *popt), 'r--', label='fit')
            plt.xlabel('Unit cell length along [1,1,0] direction (\u212B)')
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
        
            q2 = (-lambda_alpha_1_2-0.5*lambda_gamma_2)*1e6
            q4 = (lambda_alpha_1_2+0.5*lambda_gamma_2-lambda_alpha_2_2)*1e6
            q6 = 2*lambda_epsilon_2*1e6
            q8 = lambda_gamma_2*1e6
        
            print("Q2 =", q2,u'x 10\u207B\u2076')
            print(" ")
            print("Q4 =", q4,u'x 10\u207B\u2076')
            print(" ")
            print("Q6 =", q6,u'x 10\u207B\u2076')
            print(" ")
            print("Q8 =", q8,u'x 10\u207B\u2076')

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
            
            print(" ")
            print("...............")
            print(" ")        
            print("Using the convention in reference E. D. T. de Lacheisserie, Magnetostriction: Theory and Application of Magnetoelasticity (CRC Press, Boca Raton, FL, 1993).")
            print(" ")
            print("\u03BB 1\u03B1,2 =", (2*lambda_alpha_1_2+lambda_alpha_2_2)*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB 2\u03B1,2 =", (2/3)*(-lambda_alpha_1_2+lambda_alpha_2_2)*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB \u03B5 =", lambda_gamma_2*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB \u03B6 =", lambda_epsilon_2*1e6,u'x 10\u207B\u2076')
            
        
            print(" ")
            print(" ")
            print(" ")
            print("/////////////////////")
            print("Spontaneous volume magnetostriction (anisotropic contribution):")
            print("/////////////////////")
            print(" ")
            print("Hexagonal crystal with easy axis (ferromagnetic state with equal domains along all equivalent easy directions): ")
            print(" ")
            print("\u03C9_s^ani =", ((4.0/3.0)*lambda_alpha_1_2+(2.0/3.0)*lambda_alpha_2_2)*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("......................... ")
            print(" ")
            print("Hexagonal crystal with easy plane (ferromagnetic state with equal domains along all equivalent easy directions): ")
            print(" ")
            print("\u03C9_s^ani =", (-(2.0/3.0)*lambda_alpha_1_2-(1.0/3.0)*lambda_alpha_2_2)*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("/////////////////////")
            print("Polycrystal (uniform stress approximation):")
            print("/////////////////////")
            print(" ")
            print("Using the convention in reference R.R. Birss, Advances in Physics 8, 252 (1959):")
            print(" ")
            print("Hexagonal crystal with easy axis (reference demagnetized state with equal domains along all easy directions): ")
            print(" ")
            print("\u03BE =", (2.0/3.0)*q2+(4.0/15.0)*q4-(1.0/15.0)*q6+(1.0/15.0)*q8,u'x 10\u207B\u2076')
            print(" ")
            print("\u03B7 =", -(2.0/15.0)*q4+(1.0/15.0)*q6+(7.0/15.0)*q8,u'x 10\u207B\u2076')
            print(" ")
            print("......................... ")
            print(" ")
            print("Hexagonal crystal with easy plane (reference demagnetized state with equal domains along all easy directions): ")
            print(" ")
            print("\u03BE =", -(1.0/3.0)*q2-(1.0/15.0)*q4-(1.0/15.0)*q6-(4.0/15.0)*q8,u'x 10\u207B\u2076')
            print(" ")
            print("\u03B7 =", -(2.0/15.0)*q4+(1.0/15.0)*q6+(7.0/15.0)*q8,u'x 10\u207B\u2076')
            print(" ")
            print(" ")
            print("Warning!: The results calculated with -mode 1 can have significant errors, especially for non-cubic crystal symmetries. We strongly recommned to verify these results through -mode 2.")
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
                print("Warning: If these elastic constants are not the same as in the input elastic tensor file", str(self.args.elas[0]),", then check that the format of the elastic tensor is exactly the same as in the standard output file ELADAT generated by AELAS code (see Example folder)")
                print(" ")
                print(" ")
                print("Magnetoelastic constants:")
                print(" ")
                print("Using the convention in reference P. Nieves et al., Comput. Phys. Commun. 264, 107964 (2021):")
                print(" ")
                print("b21 =", str(b21), 'GPa')
                print(" ")
                print("b22 =", str(b22), 'GPa')
                print(" ")
                print("b3 =", str(b3), 'GPa')
                print(" ")
                print("b4 =", str(b4), 'GPa')
                print(" ")
                print("Using the convention in reference E. D. T. de Lacheisserie, Magnetostriction: Theory and Application of Magnetoelasticity (CRC Press, Boca Raton, FL, 1993).")
                print(" ")             
                
                b1a2 = (1/3)*(2.0*math.sqrt(2.0)*b21+math.sqrt(2.0)*b22)
                
                b2a2 = -(2/3)*(b21-b22)
                
                print("b 1\u03B1,2 =", str(b1a2), 'GPa')
                print(" ")
                print("b 2\u03B1,2 =", str(b2a2), 'GPa')
                print(" ")
                print("b \u03B5,2 =", str(b3), 'GPa')
                print(" ")
                print("b \u03B6,2 =", str(b4), 'GPa')
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
            print(" ")
            print("...............")
            print(" ")
            print("Using the convention in reference E. D. T. de Lacheisserie, Magnetostriction: Theory and Application of Magnetoelasticity (CRC Press, Boca Raton, FL, 1993).")
            print(" ")
            print(" ")
            print("\u03BB 1\u03B1,2 =", (2/3)*(2*lambda_alpha_1_2+lambda_alpha_2_2)*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB 2\u03B1,2 =", (2/3)*(-lambda_alpha_1_2+lambda_alpha_2_2)*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB \u03B3,2 =", lambda_gamma_2*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB \u03B5,2 =", lambda_epsilon_2*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB \u03B4,2 =", lambda_delta*1e6,u'x 10\u207B\u2076')
            print(" ")
            
            print(" ")
            print("Warning!: The results calculated with -mode 1 can have significant errors, especially for non-cubic crystal symmetries. We strongly recommned to verify these results through -mode 2.")
            print(" ")
            print(" ")
            print(" ")
            print("/////////////////////")
            print("Spontaneous volume magnetostriction (anisotropic contribution):")
            print("/////////////////////")
            print(" ")
            print("Tetragonal crystal with easy axis (ferromagnetic state with equal domains along all equivalent easy directions): ")
            print(" ")
            print("\u03C9_s^ani =", ((4.0/3.0)*lambda_alpha_1_2+(2.0/3.0)*lambda_alpha_2_2)*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("......................... ")
            print(" ")
            print("Tetragonal crystal with easy plane (ferromagnetic state with equal domains along all equivalent easy directions): ")
            print(" ")
            print("\u03C9_s^ani =", (-(2.0/3.0)*lambda_alpha_1_2-(1.0/3.0)*lambda_alpha_2_2)*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("/////////////////////")
            print("Polycrystal (uniform stress approximation):")
            print("/////////////////////")
            print(" ")
            print(" ")
            print("Tetragonal crystal with easy axis (reference demagnetized state with equal domains along all easy directions): ")
            print(" ")
            print("\u03BE =", 1e6*((4.0/15.0)*lambda_alpha_1_2+(1.0/15.0)*lambda_alpha_2_2-(2.0/15.0)*lambda_epsilon_2-(1.0/15.0)*lambda_gamma_2-(1.0/15.0)*lambda_delta-(1.0/3.0)*(2.0*lambda_alpha_1_2+lambda_alpha_2_2)),u'x 10\u207B\u2076')
            print(" ")
            print("\u03B7 =", 1e6*(-(2.0/15.0)*lambda_alpha_1_2+(2.0/15.0)*lambda_alpha_2_2+(2.0/5.0)*lambda_epsilon_2+(1.0/5.0)*lambda_gamma_2+(1.0/5.0)*lambda_delta),u'x 10\u207B\u2076')
            print(" ")
            print("......................... ")
            print(" ")
            print("Tetragonal crystal with easy plane (reference demagnetized state with equal domains along all easy directions): ")
            print(" ")
            print("\u03BE =", 1e6*((4.0/15.0)*lambda_alpha_1_2+(1.0/15.0)*lambda_alpha_2_2-(2.0/15.0)*lambda_epsilon_2-(1.0/15.0)*lambda_gamma_2-(1.0/15.0)*lambda_delta),u'x 10\u207B\u2076')
            print(" ")
            print("\u03B7 =", 1e6*(-(2.0/15.0)*lambda_alpha_1_2+(2.0/15.0)*lambda_alpha_2_2+(2.0/5.0)*lambda_epsilon_2+(1.0/5.0)*lambda_gamma_2+(1.0/5.0)*lambda_delta),u'x 10\u207B\u2076')
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
                print("b'3 =", str(b3p), 'GPa')
                print(" ")
                print("b4 =", str(b4), 'GPa')
                print(" ")

                print("The equation of the magnetoelastic energy can be found in the User Manual")
                print(" ")
                print("Using the convention in reference E. D. T. de Lacheisserie, Magnetostriction: Theory and Application of Magnetoelasticity (CRC Press, Boca Raton, FL, 1993):")
                print(" ")
                
                b1a2 = (1/3)*(2.0*math.sqrt(2.0)*b21+math.sqrt(2.0)*b22)
                
                b2a2 = -(2/3)*(b21-b22)
                
                print("b 1\u03B1,2 =", str(b1a2), 'GPa')
                print(" ")
                print("b 2\u03B1,2 =", str(b2a2), 'GPa')
                print(" ")
                print("b \u03B3,2 =", str(b3), 'GPa')
                print(" ")
                print("b \u03B4,2 =", str(b3p), 'GPa')
                print(" ")
                print("b \u03B5,2 =", str(b4), 'GPa')
                print(" ")

                



        
        
        