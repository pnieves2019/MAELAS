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

  def der(self):
  
        structure0 = Structure.from_file(self.args.pos[0])
        nat = len(structure0.species)
        vol = float(structure0.volume)
        
        sym1 = float(self.args.sympre[0])
        sym2 = float(self.args.symang[0])
    
        
        
        for j in range(1,10):

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
        print("Fit of quadratic function f(x)=A*x\u00B2+B*x+C to energy vs cell length")


        lambda_ortho = []


        list_spin = ['1,0,0','0,0,1','0,1,0','0,0,1','1,0,0','0,0,1','0,1,0','0,0,1','1,0,0','0,0,1','0,1,0','0,0,1','1,1,0','0,0,1','1,0,1','0,0,1','0,1,1','0,0,1']
        list_dist = ['1,0,0','1,0,0','0,1,0','0,1,0','0,0,1','0,0,1','a,b,0','a,0,c','0,b,c']

        nn = int(self.args.ndist[0])+1

        # AELAS and IEEE lattice convention: c<a<b


        aa0 = SpacegroupAnalyzer(structure0,symprec=sym1, angle_tolerance=sym2)
        structure1 = aa0.get_conventional_standard_structure(international_monoclinic=True)

        bb0 = ConventionalCellTransformation(symprec=sym1, angle_tolerance=sym2, international_monoclinic=True)
        structure2 = bb0.apply_transformation(structure1)
           
        
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
            print('Unit cell length along [', dist ,'] direction')
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
            print("A =", params[0][0], ", B =", params[0][1], ", C =", params[0][2])

            r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
            print("R-squared =", r_squared)
            print("")

            if r_squared < 0.98:
                print("WARNING!! R-squared is lower than 0.98. Check figure ", fig1)
                print("")

            l1 = -params[0][1] / (2.0 * params[0][0])

            print("X minimum = -B/(2*A) =", l1)
            print("")
            
            if i == 1:
                if nn % 2 == 0:
                    lli = int((nn-2)/2)
                    mae100 = y[lli]
            elif i == 2:
                if nn % 2 == 0:
                    lli = int((nn-2)/2)
                    mae010 = y[lli]


            plt.plot(x, y, 'bo', label=ene_dat1 )
            popt, pcov = curve_fit(K, x, y)
            t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
            plt.plot(t, K(t, *popt), 'r--', label='fit')
            plt.ylabel('Energy (eV)')
            plt.legend()
            label = "Unit cell length along [" + str(dist) + "] direction (\u212B)"
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
            print("A =", params[0][0], ", B =", params[0][1], ", C =", params[0][2])

            r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
            print("R-squared =", r_squared)
            print("")

            if r_squared < 0.98:
                print("WARNING!! R-squared is lower than 0.98. Check figure ", fig2)
                print("")

            l2 = -params[0][1] / (2.0 * params[0][0])

            print("X minimum = -B/(2*A) =", l2)
            print("")
            
            if i == 1:
                if nn % 2 == 0:
                    lli = int((nn-2)/2)
                    mae001 = y[lli]


            plt.plot(x, y, 'bo', label=ene_dat2 )
            popt, pcov = curve_fit(K, x, y)
            t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
            plt.plot(t, K(t, *popt), 'r--', label='fit')
            plt.ylabel('Energy (eV)')
            plt.legend()
            label = "Unit cell length along [" + str(dist) + "] direction (\u212B)"
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
            label = "Unit cell length along [" + str(dist) + "] direction (\u212B)"
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
                
                lmb_7_1 = ((latt_a**2.0+latt_b**2.0)*(l1 -l2))/(latt_a*latt_b*(l1+l2))
                
                lmb_7_2 =((latt_a-latt_b)*(latt_a*(lmb_1+lmb_2)-latt_b*(lmb_3+lmb_4)))/(4.0*latt_a*latt_b)
                
                lmb_7 = lmb_7_1 - lmb_7_2
                
                lambda_ortho += [lmb_7]
            
            elif i == 8:
                
                lmb_1 = lambda_ortho[0]
                
                lmb_5 = lambda_ortho[4]
                
                lmb_8_1 = ((latt_a**2.0+latt_c**2.0)*(l1 -l2))/(latt_a*latt_c*(l1+l2))
                
                lmb_8_2 =((latt_a-latt_c)*(latt_a*lmb_1-latt_c*lmb_5))/(4.0*latt_a*latt_c)
                
                lmb_8 = lmb_8_1 - lmb_8_2
                
                lambda_ortho += [lmb_8]
                   
                
            elif i == 9:
                
                lmb_4 = lambda_ortho[3]
                
                lmb_6 = lambda_ortho[5]
                
                lmb_9_1 = ((latt_b**2.0+latt_c**2.0)*(l1 -l2))/(latt_b*latt_c*(l1+l2))
                
                lmb_9_2 =((latt_b-latt_c)*(latt_b*lmb_4-latt_c*lmb_6))/(4.0*latt_b*latt_c)
                
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
            
        
        lmb1 = lambda_ortho[0]
        lmb2 = lambda_ortho[1]
        lmb3 = lambda_ortho[2]
        lmb4 = lambda_ortho[3]
        lmb5 = lambda_ortho[4]
        lmb6 = lambda_ortho[5]
        lmb7 = lambda_ortho[6]
        lmb8 = lambda_ortho[7]
        lmb9 = lambda_ortho[8]
        
        print("Using the convention in reference E. D. T. de Lacheisserie, Magnetostriction: Theory and Application of Magnetoelasticity (CRC Press, Boca Raton, FL, 1993):")
        print(" ")
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
        
        print("Warning!: The results calculated with -mode 1 can have significant errors, especially for non-cubic crystal symmetries. We strongly recommned to verify these results through -mode 2.")
        print(" ")
        print(" ")
        print("/////////////////////")
        print("Polycrystal (uniform stress approximation):")
        print("/////////////////////")
        print(" ")
        print(" ")
        print("Orthorhombic crystal with easy axis along lattice vector c (reference demagnetized state with equal domains along all easy directions): ")
        print(" ")
        print("\u03BE =", 1e6*((2.0/15.0)*(lambda_ortho[0]+lambda_ortho[3]-lambda_ortho[6]-lambda_ortho[7]-lambda_ortho[8])+(1.0/6.0)*(lambda_ortho[1]+lambda_ortho[2]+lambda_ortho[4]+lambda_ortho[5])),u'x 10\u207B\u2076')
        print(" ")
        print("\u03B7 =", 1e6*(-(1.0/15.0)*lambda_ortho[0]-(1.0/15.0)*lambda_ortho[3]-(1.0/6.0)*(lambda_ortho[1]+lambda_ortho[2]+lambda_ortho[4]+lambda_ortho[5])+(2.0/5.0)*(lambda_ortho[6]+lambda_ortho[7]+lambda_ortho[8])),u'x 10\u207B\u2076')
        print(" ")
        print("......................... ")
        print(" ")
        print("Orthorhombic crystal with easy axis along lattice vector a (reference demagnetized state with equal domains along all easy directions): ")
        print(" ")
        print("\u03BE =", 1e6*((2.0/15.0)*(-(3.0/2.0)*lambda_ortho[0]+lambda_ortho[3]-lambda_ortho[6]-lambda_ortho[7]-lambda_ortho[8])+(1.0/6.0)*(lambda_ortho[1]-lambda_ortho[2]-lambda_ortho[4]+lambda_ortho[5])),u'x 10\u207B\u2076')
        print(" ")
        print("\u03B7 =", 1e6*(-(1.0/15.0)*lambda_ortho[0]-(1.0/15.0)*lambda_ortho[3]-(1.0/6.0)*(lambda_ortho[1]+lambda_ortho[2]+lambda_ortho[4]+lambda_ortho[5])+(2.0/5.0)*(lambda_ortho[6]+lambda_ortho[7]+lambda_ortho[8])),u'x 10\u207B\u2076')
        print(" ")
        print("......................... ")
        print(" ")
        print("Orthorhombic crystal with easy axis along lattice vector b (reference demagnetized state with equal domains along all easy directions): ")
        print(" ")
        print("\u03BE =", 1e6*((2.0/15.0)*(lambda_ortho[0]-(3.0/2.0)*lambda_ortho[3]-lambda_ortho[6]-lambda_ortho[7]-lambda_ortho[8])-(1.0/6.0)*(lambda_ortho[1]-lambda_ortho[2]-lambda_ortho[4]+lambda_ortho[5])),u'x 10\u207B\u2076')
        print(" ")
        print("\u03B7 =", 1e6*(-(1.0/15.0)*lambda_ortho[0]-(1.0/15.0)*lambda_ortho[3]-(1.0/6.0)*(lambda_ortho[1]+lambda_ortho[2]+lambda_ortho[4]+lambda_ortho[5])+(2.0/5.0)*(lambda_ortho[6]+lambda_ortho[7]+lambda_ortho[8])),u'x 10\u207B\u2076')
        print(" ")
        print(" ")
        print(" ")
        
            
        
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
                print("Calculation of magnetoelastic constants:")
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
                print("Using the convention in reference E. D. T. de Lacheisserie, Magnetostriction: Theory and Application of Magnetoelasticity (CRC Press, Boca Raton, FL, 1993):")
                print(" ")
                b1a2 = -(b1+b2+b3+b4+b5+b6)/(3.0*math.sqrt(2))
                b2a2 = (1.0/6.0)*(b1+b2+b3+b4-2.0*(b5+b6))
                b3a2 = (-b1-b2+b3+b4)/(math.sqrt(6))
                b1a2p = (b1-b2+b3-b4+b5-b6)/(math.sqrt(6))
                b2a2p = (-b1+b2-b3+b4+2.0*b5-2.0*b6)/(2.0*math.sqrt(3)) 
                b3a2p = (1.0/2.0)*(b1-b2-b3+b4)
                bb2 = b7
                bg2 = b9
                bd2 = b8
                print("b1\u03B1,2 =", str(b1a2), 'GPa')
                print(" ")
                print("b2\u03B1,2 =", str(b2a2), 'GPa')
                print(" ")
                print("b3\u03B1,2 =", str(b3a2), 'GPa')
                print(" ")
                print("b1\u03B1,2' =", str(b1a2p), 'GPa')
                print(" ")
                print("b2\u03B1,2' =", str(b2a2p), 'GPa')
                print(" ")
                print("b3\u03B1,2' =", str(b3a2p), 'GPa')
                print(" ")
                print("b\u03B2,2 =", str(bb2), 'GPa')
                print(" ")
                print("b\u03B3,2 =", str(bg2), 'GPa')
                print(" ")
                print("b\u03B4,2 =", str(bd2), 'GPa')
                print(" ")

