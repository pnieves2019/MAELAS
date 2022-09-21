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

        
        # AELAS and IEEE lattice convention: c<a<b


        latt0 = structure2.lattice.matrix
        coordsnew = np.zeros((len(structure2.species), 3))

        Listlatt0 = [latt0[0][0],latt0[1][1],latt0[2][2]]
        Listlattnew = sorted(Listlatt0)
        
        for ii in range(len(Listlattnew)):
            if Listlattnew[0] == Listlatt0[ii]:
                indmin = ii
            if Listlattnew[1] == Listlatt0[ii]:
                indmid = ii
            if Listlattnew[2] == Listlatt0[ii]:
                indmax = ii
        
        
        
        for i in range(len(structure2.species)):
            coordsnew[i][0] = float(structure2.frac_coords[i][indmid])
            coordsnew[i][1] = float(structure2.frac_coords[i][indmax])
            coordsnew[i][2] = float(structure2.frac_coords[i][indmin])


        lattice = Lattice.from_parameters(a=latt0[indmid][indmid], b=latt0[indmax][indmax], c=latt0[indmin][indmin], alpha=90, beta=90, gamma=90)
        structure2b = Structure(lattice, structure2.species, coordsnew)


        for i in range(int(self.args.ndist[0])):


            strain1 = - float(self.args.strain[0])+2*(float(self.args.strain[0])/(float(self.args.ndist[0])-1))*i

            print("strain", strain1)


        #Generation POSCAR file

        #lambda_1


            a1 = 1.0 + strain1
            a2 = 1.0
            a3 = a2
            dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
            structure3 = dd.apply_transformation(structure2b)
            pos_name3 = "POSCAR_1_" + str(i+1)

            structure33 = Poscar(structure3)
            structure33.write_file(filename = pos_name3,significant_figures=16)

        #lambda_2


            a1 = 1.0 + strain1
            a2 = 1.0
            a3 = a2
            dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
            structure4 = dd.apply_transformation(structure2b)
            pos_name4 = "POSCAR_2_" + str(i+1)

            structure44 = Poscar(structure4)
            structure44.write_file(filename = pos_name4,significant_figures=16)

        #lambda_3

            a2 = 1.0 + strain1
            a1 = 1.0
            a3 = a1
            dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
            structure5 = dd.apply_transformation(structure2b)
            pos_name5 = "POSCAR_3_" + str(i+1)

            structure55 = Poscar(structure5)
            structure55.write_file(filename = pos_name5,significant_figures=16)

        #lambda_4

            a2 = 1.0 + strain1
            a1 = 1.0
            a3 = a1
            dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
            structure6 = dd.apply_transformation(structure2b)
            pos_name6 = "POSCAR_4_" + str(i+1)

            structure66 = Poscar(structure6)
            structure66.write_file(filename = pos_name6,significant_figures=16)


       #lambda_5

            a3 = 1.0 + strain1
            a1 = 1.0
            a2 = a1
            dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
            structure7 = dd.apply_transformation(structure2b)
            pos_name7 = "POSCAR_5_" + str(i+1)

            structure77 = Poscar(structure7)
            structure77.write_file(filename = pos_name7,significant_figures=16)

       #lambda_6

            a3 = 1.0 + strain1
            a1 = 1.0
            a2 = a1
            dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
            structure8 = dd.apply_transformation(structure2b)
            pos_name8 = "POSCAR_6_" + str(i+1)

            structure88 = Poscar(structure8)
            structure88.write_file(filename = pos_name8,significant_figures=16)

        #lambda_7

            
            a11 = 1.0
            a12 = strain1
            a13 = 0.0
            a21 = strain1
            a22 = 1.0
            a23 = 0.0
            a31 = 0.0
            a32 = 0.0
            a33 = 1.0
            

            cc = DeformStructureTransformation(deformation=((a11, a12, a13), (a21, a22, a23), (a31, a32, a33)))
            structure9 = cc.apply_transformation(structure2b)
            pos_name9 = "POSCAR_7_" + str(i+1)

            structure99 = Poscar(structure9)
            structure99.write_file(filename = pos_name9,significant_figures=16)

        #lambda_8
                 
            
            a11 = 1.0
            a12 = 0.0
            a13 = strain1
            a21 = 0.0
            a22 = 1.0
            a23 = 0.0
            a31 = strain1
            a32 = 0.0
            a33 = 1.0

            cc = DeformStructureTransformation(deformation=((a11, a12, a13), (a21, a22, a23), (a31, a32, a33)))
            structure10 = cc.apply_transformation(structure2b)
            pos_name10 = "POSCAR_8_" + str(i+1)

            structure1010 = Poscar(structure10)
            structure1010.write_file(filename = pos_name10,significant_figures=16)


        #lambda_9

            a11 = 1.0
            a12 = 0.0
            a13 = 0.0
            a21 = 0.0
            a22 = 1.0
            a23 = strain1
            a31 = 0.0
            a32 = strain1
            a33 = 1.0

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

  def der(self):
  
        structure0 = Structure.from_file(self.args.pos[0])
        nat = len(structure0.species)
        vol = float(structure0.volume)
        
        sym1 = float(self.args.sympre[0])
        sym2 = float(self.args.symang[0])
    
        aa = SpacegroupAnalyzer(structure0,symprec=sym1, angle_tolerance=sym2)
        
        
        for j in range(1,10):

            for k in range(1,3):

                path_dat = "ene_" + str(j) + "_" + str(k) + ".dat"
                dat = open(path_dat,'w')


                for i in range(int(self.args.ndist[0])):
                    
                    strain1 = - float(self.args.strain[0])+2*(float(self.args.strain[0])/(float(self.args.ndist[0])-1))*i

                    #print("strain", strain1)

                    pos_name = "POSCAR_" + str(j) + "_" + str(i+1)

                    struct = Structure.from_file(pos_name)

                    latt = struct.lattice.matrix


                    path_osz = "OSZICAR_" + str(j) + "_" + str(i+1) + "_" + str(k)
                    osz = open(path_osz,'r')
                    ene0 = osz.readlines()
                    ene1 = ene0[len(ene0)-2]
                    ene2 = ene1[11:32]

                    osz.close()

                    dat.write(repr(strain1))
                    dat.write('  ')
                    dat.write(str(ene2))
                    dat.write('\n')


                dat.close()


       # fitting and plot


        def K(x,a,b):
            return a*x+b


        print("")
        print("Fit of linear function f(x)=A*x+B to energy vs strain")


        lambda_ortho = []
        b = []


        list_spin = ['1,0,0','0,0,1','0,1,0','0,0,1','1,0,0','0,0,1','0,1,0','0,0,1','1,0,0','0,0,1','0,1,0','0,0,1','1,1,0','0,0,1','1,0,1','0,0,1','0,1,1','0,0,1']
        list_dist = ['\u03B5xx','\u03B5xx','\u03B5yy','\u03B5yy','\u03B5zz','\u03B5zz','\u03B5xy','\u03B5xz','\u03B5yz']

        nn = int(self.args.ndist[0])+1

        
        
        for i in range(1,10):


            ene_dat1 = "ene_" + str(i) + "_1.dat"
            ene_dat2 = "ene_" + str(i) + "_2.dat"
            spin1 = str(list_spin[2*i-2])
            spin2 = str(list_spin[2*i-1])
            dist = str(list_dist[i-1])
            fig = 'dE_' + str(i) + '.png'
            
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
              
            params = curve_fit(K, x, y-y2)

            print("Fitting parameters for the calculation of b",i," :")
            print("A =", params[0][0], ", B =", params[0][1])

            r_squared = r2_score(y-y2, K(x,params[0][0],params[0][1]))
            print("R-squared =", r_squared)
            print("")

            if r_squared < 0.98:
                print("WARNING!! R-squared is lower than 0.98. Check figure dE_",i,".png")
                print("")
  
            b += [params[0][0]]
            
            if i == 1:
                if nn % 2 == 0:
                    lli = int((nn-2)/2)
                    mae100 = y[lli]
                    mae001 = y2[lli]
            elif i == 2:
                if nn % 2 == 0:
                    lli = int((nn-2)/2)
                    mae010 = y[lli]
                    
            
            
            #make figure dE_i.png

            tit = 'Calculation of b' + str(i) 

            plt.plot(x, (y-y2)*1e6, 'o', label='data')
            plt.plot(x, K(x, params[0][0],params[0][1])*1e6, 'r-', label='fit')
            plt.legend()
            ax4 = plt.gca()
            ax4.xaxis.set_major_locator(plt.MaxNLocator(5))

            ylabel ='E[' + str(spin1) + '] - E['+ str(spin2) + '] (\u03BCeV)'
            plt.ylabel(ylabel)
            label = "Strain " + str(dist)
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
        print("Using the convention in reference P.Nieves et al., Comput. Phys. Commun., 264, 107964  (2021):")
        print(" ")

        for i in range(len(b)):

            print("b",i+1," =", b[i],u' eV')
            print(" ")
        print("........... ")    
        for i in range(len(b)):

            print("b",i+1," =", b[i]/nat,u' eV/atom')
            print(" ")
        print("........... ")     
        for i in range(len(b)):

            print("b",i+1," =", b[i]/vol,u' eV/A^3')
            print(" ")
        print("........... ")    
        for i in range(len(b)):

            print("b",i+1," =", (b[i]/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
            print(" ")
            
            
        print(" ")
        print("Using the convention in reference E. D. T. de Lacheisserie, Magnetostriction: Theory and Application of Magnetoelasticity (CRC Press, Boca Raton, FL, 1993):")
        print(" ")
        
        b1 = b[0]
        b2 = b[1]
        b3 = b[2]
        b4 = b[3]
        b5 = b[4]
        b6 = b[5]
        b7 = b[6]
        b8 = b[7]
        b9 = b[8]
        
        b1a2 = -(b1+b2+b3+b4+b5+b6)/(3.0*math.sqrt(2))
        b2a2 = (1.0/6.0)*(b1+b2+b3+b4-2.0*(b5+b6))
        b3a2 = (-b1-b2+b3+b4)/(math.sqrt(6))
        b1a2p = (b1-b2+b3-b4+b5-b6)/(math.sqrt(6))
        b2a2p = (-b1+b2-b3+b4+2.0*b5-2.0*b6)/(2.0*math.sqrt(3)) 
        b3a2p = (1.0/2.0)*(b1-b2-b3+b4)
        bb2 = b7
        bg2 = b9
        bd2 = b8
        print("b1\u03B1,2 =", (b1a2/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
        print(" ")
        print("b2\u03B1,2 =", (b2a2/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
        print(" ")
        print("b3\u03B1,2 =", (b3a2/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
        print(" ")
        print("b1\u03B1,2' =", (b1a2p/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
        print(" ")
        print("b2\u03B1,2' =", (b2a2p/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
        print(" ")
        print("b3\u03B1,2' =", (b3a2p/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
        print(" ")
        print("b\u03B2,2 =", (bb2/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
        print(" ")
        print("b\u03B3,2 =", (bg2/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
        print(" ")
        print("b\u03B4,2 =", (bd2/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
        print(" ")
                      
        
        bani1 = (b1a2/vol)*1.602176565e-19*1e30*1e-6
        bani2 = (b2a2/vol)*1.602176565e-19*1e30*1e-6
        bani3 = (b3a2/vol)*1.602176565e-19*1e30*1e-6
        bani4 = (b1a2p/vol)*1.602176565e-19*1e30*1e-6
        bani5 = (b2a2p/vol)*1.602176565e-19*1e30*1e-6
        bani6 = (b3a2p/vol)*1.602176565e-19*1e30*1e-6
        bani7 = (bb2/vol)*1.602176565e-19*1e30*1e-6
        bani8 = (bg2/vol)*1.602176565e-19*1e30*1e-6
        bani9 = (bd2/vol)*1.602176565e-19*1e30*1e-6
        
        path_inc_std = 'MAGANI'
        inc_std = open(path_inc_std,'w')
        inc_std.write(str(bani1))
        inc_std.write('\n')
        inc_std.write(str(bani2))
        inc_std.write('\n')
        inc_std.write(str(bani3))
        inc_std.write('\n')
        inc_std.write(str(bani4))
        inc_std.write('\n')
        inc_std.write(str(bani5))
        inc_std.write('\n')
        inc_std.write(str(bani6))
        inc_std.write('\n')
        inc_std.write(str(bani7))
        inc_std.write('\n')
        inc_std.write(str(bani8))
        inc_std.write('\n')
        inc_std.write(str(bani9))
        inc_std.close()
        
            
        
        if nn % 2 == 0:         
            print("----------------------------------------------")
            print("Magnetocrystalline anisotropy energy:")
            print("----------------------------------------------")
            print(" ")
            print("These energies correspond to the central points in the data files ene_1_1.dat, ene_2_1.dat, and ene_1_2.dat:")
            print(" ")
            print("E(1,0,0) = ",mae100," eV")
            print(" ")
            print("E(0,1,0) = ",mae010," eV")
            print(" ")
            print("E(0,0,1) = ",mae001," eV")
            print(" ")
            print("E(1,0,0) - E(0,0,1) = ",(mae100 - mae001)*1e6,u'x 10\u207B\u2076 eV')
            print(" ")
            print("[E(1,0,0) - E(0,0,1)]/Natom = ",((mae100 - mae001)/nat)*1e6,u'x 10\u207B\u2076 eV/atom')
            print(" ")
            print("E(0,1,0) - E(0,0,1) = ",(mae010 - mae001)*1e6,u'x 10\u207B\u2076 eV')
            print(" ")
            print("[E(0,1,0) - E(0,0,1)]/Natom = ",((mae010 - mae001)/nat)*1e6,u'x 10\u207B\u2076 eV/atom')
            print(" ")
        

        if self.args.delas == True:

                print(" ")
                print(" ")
                print("----------------------------------------------")
                print("Calculation of anisotropic magnetostrictive coefficients:")
                print("----------------------------------------------")
                print(" ")
                print("Reading the elastic tensor file =", str(self.args.elas[0]))
                print(" ")

                elasdat = open(self.args.elas[0],'r')
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


                

                b1 = (b[0]/vol)*1.602176565e-19*1e30*1e-9
                b2 = (b[1]/vol)*1.602176565e-19*1e30*1e-9
                b3 = (b[2]/vol)*1.602176565e-19*1e30*1e-9
                b4 = (b[3]/vol)*1.602176565e-19*1e30*1e-9
                b5 = (b[4]/vol)*1.602176565e-19*1e30*1e-9
                b6 = (b[5]/vol)*1.602176565e-19*1e30*1e-9
                b7 = (b[6]/vol)*1.602176565e-19*1e30*1e-9
                b8 = (b[7]/vol)*1.602176565e-19*1e30*1e-9
                b9 = (b[8]/vol)*1.602176565e-19*1e30*1e-9
                
                

                denom = c13**2.0*c22-2.0*c12*c13*c23+c12**2.0*c33+c11*(c23**2.0-c22*c33)

                lambda_ortho1 = (-b5*c13*c22 + b5*c12*c23 + b3*c13*c23 - b1*c23**2.0 -b3*c12*c33 + b1*c22*c33)/denom

                lambda_ortho2 = (-b6*c13*c22 + b6*c12*c23 + b4*c13*c23 - b2*c23**2.0 - b4*c12*c33 + b2*c22*c33)/denom

                lambda_ortho3 = (b5*c12*c13 - b3*c13**2.0 - b5*c11*c23 + b1*c13*c23 + b3*c11*c33 - b1*c12*c33)/denom

                lambda_ortho4 = (b6*c12*c13 - b4*c13**2.0 - b6*c11*c23 + b2*c13*c23 + b4*c11*c33 - b2*c12*c33)/denom

                lambda_ortho5 = (-b5*c12**2.0 + b3*c12*c13 + b5*c11*c22 - b1*c13*c22 - b3*c11*c23 + b1*c12*c23)/denom

                lambda_ortho6 = (-b6*c12**2.0 + b4*c12*c13 + b6*c11*c22 - b2*c13*c22 - b4*c11*c23 + b2*c12*c23)/denom

                lambda_ortho7 = ((c13-c23)*((b3+b4)*c13 - (b1+b2)*c23) + b5*(c13*c22+c11*c23-c12*(c13+c23)) + b6*(c13*c22+c11*c23-c12*(c13+c23))+(-(b3+b4)*(c11-c12)+(b1+b2)*(c12-c22))*c33)/(-4.0*denom)

                lambda_ortho7 = lambda_ortho7 - (b7/(4.0*c66))

                lambda_ortho8 = (b5*(c11-c13)*c22-b3*c11*c23+b1*(c12-c23)*c23+b5*c12*(-c12+c23)+b3*c13*(c12+c23)-b3*c12*c33+b1*c22*(-c13+c33))/(4.0*denom)

                lambda_ortho8 = lambda_ortho8 - (b8/(4.0*c55))

                lambda_ortho9 = (b4*(c12-c13)*c13+b6*c12*(-c12+c13)+b6*c11*(c22-c23)+b2*c13*(-c22+c23)+b2*c12*(c23-c33)+b4*c11*(-c23+c33))/(4.0*denom)

                lambda_ortho9 = lambda_ortho9 - (b9/(4.0*c44))

                


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
                print("Warning: If these elastic constants are not the same as in the input elastic tensor file", str(self.args.elas[0]),", then check that the format of the elastic tensor is exactly the same as in the standard output file ELADAT generated by AELAS code (see Example folder)")
                print(" ")
                print(" ")
                print("Anisotropic magnetostrictive coefficients:")
                print(" ")
                print("Using the convention in reference P. Nieves et al., Comput. Phys. Commun. 264, 107964 (2021):")
                print(" ")
                print("\u03BB1 =", lambda_ortho1*1e6,u'x 10\u207B\u2076')
                print(" ")
                print("\u03BB2 =", lambda_ortho2*1e6,u'x 10\u207B\u2076')
                print(" ")
                print("\u03BB3 =", lambda_ortho3*1e6,u'x 10\u207B\u2076')
                print(" ")
                print("\u03BB4 =", lambda_ortho4*1e6,u'x 10\u207B\u2076')
                print(" ")
                print("\u03BB5 =", lambda_ortho5*1e6,u'x 10\u207B\u2076')
                print(" ")
                print("\u03BB6 =", lambda_ortho6*1e6,u'x 10\u207B\u2076')
                print(" ")
                print("\u03BB7 =", lambda_ortho7*1e6,u'x 10\u207B\u2076')
                print(" ")
                print("\u03BB8 =", lambda_ortho8*1e6,u'x 10\u207B\u2076')
                print(" ")
                print("\u03BB9 =", lambda_ortho9*1e6,u'x 10\u207B\u2076')
                print(" ")
                print("Using the convention in reference E. D. T. de Lacheisserie, Magnetostriction: Theory and Application of Magnetoelasticity (CRC Press, Boca Raton, FL, 1993):")
                print(" ")
                
                
                lmb1 = lambda_ortho1
                lmb2 = lambda_ortho2
                lmb3 = lambda_ortho3
                lmb4 = lambda_ortho4
                lmb5 = lambda_ortho5
                lmb6 = lambda_ortho6
                lmb7 = lambda_ortho7
                lmb8 = lambda_ortho8
                lmb9 = lambda_ortho9
        

                print("\u03BB 1\u03B1,2 =", (1.0/3.0)*(-lmb1-lmb2-lmb3-lmb4-lmb5-lmb6)*1e6,u'x 10\u207B\u2076')
                print(" ")
                print("\u03BB 2\u03B1,2 =", (1.0/6.0)*(lmb1+lmb2+lmb3+lmb4-2.0*(lmb5+lmb6))*1e6,u'x 10\u207B\u2076')
                print(" ")
                print("\u03BB 3\u03B1,2 =", (1.0/6.0)*(-lmb1-lmb2+lmb3+lmb4)*1e6,u'x 10\u207B\u2076')
                print(" ")
                print("\u03BB 1\u03B1,2' =", (lmb1-lmb2+lmb3-lmb4+lmb5-lmb6)*1e6,u'x 10\u207B\u2076')
                print(" ")
                print("\u03BB 2\u03B1,2' =", (1.0/2.0)*(-lmb1+lmb2-lmb3+lmb4+2.0*lmb5-2.0*lmb6)*1e6,u'x 10\u207B\u2076')
                print(" ")
                print("\u03BB 3\u03B1,2' =", (1.0/2.0)*(lmb1-lmb2-lmb3+lmb4)*1e6,u'x 10\u207B\u2076')
                print(" ")
                print("\u03BB \u03B2,2 =", (1.0/2.0)*(-lmb1-lmb2-lmb3-lmb4+4.0*lmb7)*1e6,u'x 10\u207B\u2076')
                print(" ")
                print("\u03BB \u03B3,2 =", (1.0/2.0)*(-lmb4-lmb6+4.0*lmb9)*1e6,u'x 10\u207B\u2076')
                print(" ")
                print("\u03BB \u03B4,2 =", (1.0/2.0)*(-lmb1-lmb5+4.0*lmb8)*1e6,u'x 10\u207B\u2076')
                print(" ")
                
                
                print(" ")
                print("/////////////////////")
                print("Spontaneous volume magnetostriction (anisotropic contribution):")
                print("/////////////////////")
                print(" ")
        
                lmb1a2 = (1.0/3.0)*(-lmb1-lmb2-lmb3-lmb4-lmb5-lmb6)
                lmb1a2p = lmb1-lmb2+lmb3-lmb4+lmb5-lmb6
        
                print("Orthorhombic crystal with easy axis along lattice vector c (ferromagnetic state with equal domains along all equivalent easy directions): ")
                print(" ")
                print("\u03C9_s^ani =", lmb1a2*1e6,u'x 10\u207B\u2076')
                print(" ")
                print("......................... ")
                print(" ")
                print("Orthorhombic crystal with easy axis along lattice vector a (ferromagnetic state with equal domains along all equivalent easy directions): ")
                print(" ")
                print("\u03C9_s^ani =", (-0.5*lmb1a2+0.5*lmb1a2p)*1e6,u'x 10\u207B\u2076')
                print(" ")
                print("......................... ")
                print(" ")
                print("Orthorhombic crystal with easy axis along lattice vector b (ferromagnetic state with equal domains along all equivalent easy directions): ")
                print(" ")
                print("\u03C9_s^ani =", (-0.5*lmb1a2-0.5*lmb1a2p)*1e6,u'x 10\u207B\u2076')
                print(" ")
                print(" ")
                print("/////////////////////")
                print("Polycrystal (uniform stress approximation):")
                print("/////////////////////")
                print(" ")
                print(" ")
                print("Orthorhombic crystal with easy axis along lattice vector c (reference demagnetized state with equal domains along all easy directions): ")
                print(" ")
                print("\u03BE =", 1e6*((2.0/15.0)*(lambda_ortho1+lambda_ortho4-lambda_ortho7-lambda_ortho8-lambda_ortho9)+(1.0/6.0)*(lambda_ortho2+lambda_ortho3+lambda_ortho5+lambda_ortho6)),u'x 10\u207B\u2076')
                print(" ")
                print("\u03B7 =", 1e6*(-(1.0/15.0)*lambda_ortho1-(1.0/15.0)*lambda_ortho4-(1.0/6.0)*(lambda_ortho2+lambda_ortho3+lambda_ortho5+lambda_ortho6)+(2.0/5.0)*(lambda_ortho7+lambda_ortho8+lambda_ortho9)),u'x 10\u207B\u2076')
                print(" ")
                print("......................... ")
                print(" ")
                print("Orthorhombic crystal with easy axis along lattice vector a (reference demagnetized state with equal domains along all easy directions): ")
                print(" ")
                print("\u03BE =", 1e6*((2.0/15.0)*(-(3.0/2.0)*lambda_ortho1+lambda_ortho4-lambda_ortho7-lambda_ortho8-lambda_ortho9)+(1.0/6.0)*(lambda_ortho2-lambda_ortho3-lambda_ortho5+lambda_ortho6)),u'x 10\u207B\u2076')
                print(" ")
                print("\u03B7 =", 1e6*(-(1.0/15.0)*lambda_ortho1-(1.0/15.0)*lambda_ortho4-(1.0/6.0)*(lambda_ortho2+lambda_ortho3+lambda_ortho5+lambda_ortho6)+(2.0/5.0)*(lambda_ortho7+lambda_ortho8+lambda_ortho9)),u'x 10\u207B\u2076')
                print(" ")
                print("......................... ")
                print(" ")
                print("Orthorhombic crystal with easy axis along lattice vector b (reference demagnetized state with equal domains along all easy directions): ")
                print(" ")
                print("\u03BE =", 1e6*((2.0/15.0)*(lambda_ortho1-(3.0/2.0)*lambda_ortho4-lambda_ortho7-lambda_ortho8-lambda_ortho9)-(1.0/6.0)*(lambda_ortho2-lambda_ortho3-lambda_ortho5+lambda_ortho6)),u'x 10\u207B\u2076')
                print(" ")
                print("\u03B7 =", 1e6*(-(1.0/15.0)*lambda_ortho1-(1.0/15.0)*lambda_ortho4-(1.0/6.0)*(lambda_ortho2+lambda_ortho3+lambda_ortho5+lambda_ortho6)+(2.0/5.0)*(lambda_ortho7+lambda_ortho8+lambda_ortho9)),u'x 10\u207B\u2076')
                print(" ")
                print(" ")
                print(" ")


        
        
        