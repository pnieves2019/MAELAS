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
        
        sg = aa.get_space_group_number()
        
        if self.args.noconv == False:
          structure1 = aa.get_conventional_standard_structure(international_monoclinic=True)
          bb = ConventionalCellTransformation(symprec=sym1, angle_tolerance=sym2, international_monoclinic=True)
          structure2 = bb.apply_transformation(structure1)
          
        
        if self.args.noconv == True:
          structure2 = structure0
          
          
        for i in range(int(self.args.ndist[0])):


            strain1 = - float(self.args.strain[0])+2*(float(self.args.strain[0])/(float(self.args.ndist[0])-1))*i

            print("strain", strain1)


        #Generation POSCAR file Cubic I
        
            if 230 >= sg >= 207:

        #b_alpha,2


              a1 = 1.0 + strain1
              a2 = 1.0 + strain1
              a3 = 1.0 + strain1
              dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
              structure3 = dd.apply_transformation(structure2)
              pos_name = "POSCAR_1_" + str(i+1)

              structure33 = Poscar(structure3)
              structure33.write_file(filename = pos_name,significant_figures=16)
              
        #Generation POSCAR file Hexagonal I, Trigonal I or Tetragonal I:
        
            elif (177 <= sg <= 194) or (149 <= sg <= 167) or (89 <= sg <= 142):
            
              if (177 <= sg <= 194) or (149 <= sg <= 167):
            # Convention: lattice vector a1 along x-axis
                angle = -math.pi*(60.0/180.0)
                dd = DeformStructureTransformation(deformation=((math.cos(angle), math.sin(angle), 0), (-math.sin(angle), math.cos(angle), 0), (0, 0, 1)))
                structure2b = dd.apply_transformation(structure2)
              else:
                structure2b = structure2

        #b_11


              a1 = 1.0 + strain1
              a2 = 1.0 + strain1
              a3 = 1.0
              dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
              structure3 = dd.apply_transformation(structure2b)
              pos_name = "POSCAR_1_" + str(i+1)

              structure33 = Poscar(structure3)
              structure33.write_file(filename = pos_name,significant_figures=16)
              
        #b_12


              a1 = 1.0
              a2 = 1.0
              a3 = 1.0 + strain1
              dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
              structure4 = dd.apply_transformation(structure2b)
              pos_name = "POSCAR_2_" + str(i+1)

              structure44 = Poscar(structure4)
              structure44.write_file(filename = pos_name,significant_figures=16)



        #Generation POSCAR file Orthorhombic:
        
            elif 16 <= sg <= 74:
            
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
        
        
        
              for iii in range(len(structure2.species)):
                coordsnew[iii][0] = float(structure2.frac_coords[iii][indmid])
                coordsnew[iii][1] = float(structure2.frac_coords[iii][indmax])
                coordsnew[iii][2] = float(structure2.frac_coords[iii][indmin])


              lattice = Lattice.from_parameters(a=latt0[indmid][indmid], b=latt0[indmax][indmax], c=latt0[indmin][indmin], alpha=90, beta=90, gamma=90)
              structure2b = Structure(lattice, structure2.species, coordsnew)
            

        #b_1alpha,0


              a1 = 1.0 + strain1
              a2 = 1.0 + strain1
              a3 = 1.0 + strain1
              dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
              structure3 = dd.apply_transformation(structure2b)
              pos_name = "POSCAR_1_" + str(i+1)

              structure33 = Poscar(structure3)
              structure33.write_file(filename = pos_name,significant_figures=16)
              
        #b_2alpha,0


              a1 = 1.0 - 0.5*strain1
              a2 = 1.0 - 0.5*strain1
              a3 = 1.0 + strain1
              dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
              structure4 = dd.apply_transformation(structure2b)
              pos_name = "POSCAR_2_" + str(i+1)

              structure44 = Poscar(structure4)
              structure44.write_file(filename = pos_name,significant_figures=16)

        #b_3alpha,0


              a1 = 1.0 + strain1
              a2 = 1.0 - strain1
              a3 = 1.0 
              dd = DeformStructureTransformation(deformation=((a1, 0, 0), (0, a2, 0), (0, 0, a3)))
              structure5 = dd.apply_transformation(structure2b)
              pos_name = "POSCAR_3_" + str(i+1)

              structure55 = Poscar(structure4)
              structure55.write_file(filename = pos_name,significant_figures=16)


        
        

 # Derivation of magnetostriction coefficients:       

  def der(self):
  
        structure0 = Structure.from_file(self.args.pos[0])
        nat = len(structure0.species)
        vol = float(structure0.volume)
        
        sym1 = float(self.args.sympre[0])
        sym2 = float(self.args.symang[0])
    
        aa = SpacegroupAnalyzer(structure0,symprec=sym1, angle_tolerance=sym2)
        
        
        sg = aa.get_space_group_number()
        pg = aa.get_point_group_symbol()
        
        if 230 >= sg >= 207:
          print("Cubic (I) system")
          print("Number of isotropic magnetoelastic constants =", 1)
          jval=2
        elif 177 <= sg <= 194:
          print("Hexagonal (I) system")
          print("Number of isotropic magnetoelastic constants =", 2)
          jval=3
        elif 149 <= sg <= 167:
          print("Trigonal (I) system")
          print("Number of isotropic magnetoelastic constants =", 2)
          jval=3
        elif 89 <= sg <= 142:
          print("Tetragonal (I) system")
          print("Number of isotropic magnetoelastic constants =", 2)
          jval=3
        elif 16 <= sg <= 74:
          print("Orthorhombic system")
          print("Number of isotropic magnetoelastic constants =", 3)
          jval=4
    
        print("Point group =", str(pg))
        
        for j in range(1,int(jval)):

            for k in range(1,2):

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


        def K(x,a,b,c):
            return a*x**2+b*x+c
            

        #Cubic I
        
        if 230 >= sg >= 207:

          print("")
          print("Fit of quadractic function f(x)=A*x^2+B*x+C to energy vs strain")
          print("")
          print("-------------------------")
          print("Calculation of b\u03B1,2:")
          print("-------------------------")
          print(" ")
          print('Deformation x-direction (\u03B5xx),  y-direction (\u03B5yy) and  z-direction (\u03B5zz) , spin1=[0,0,1]')
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

          print("Fitting parameters:")
          print("A =", params[0][0], ", B =", params[0][1], ", C =", params[0][2])

          r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
          print("R-squared =", r_squared)
          print("")

          if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_1_1.png")
            print("")

          balpha2 = params[0][1]



        #make figure dE_1.png

          fig = 'fit_ene_1_1.png'
          tit = "Calculation of b\u03B1,2"

          plt.plot(x, y, 'o', label='data in fit_ene_1_1.dat')
          t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
          plt.plot(t, K(t, params[0][0],params[0][1],params[0][2]), 'r-', label='fit')
          plt.legend()
          ax = plt.gca()
          ax.xaxis.set_major_locator(plt.MaxNLocator(5))

          ylabel ='Energy (eV)'
          plt.ylabel(ylabel)
          label = "Strain \u03B5xx , \u03B5yy, \u03B5zz"
          plt.xlabel(label)
          plt.title(tit)
          plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
          plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
          plt.savefig(fig)
          plt.close()



          print(" ")
          print("----------------------------------------------")
          print("Isotropic magnetoelastic constant:")
          print("----------------------------------------------")
          print(" ")
          print(" ")
          print("Using the convention in reference E. D. T. de Lacheisserie, Magnetostriction: Theory and Application of Magnetoelasticity (CRC Press, Boca Raton, FL, 1993).")
          print(" ")
          print("b \u03B1,2 =", balpha2,u' eV')
          print(" ")
          print("b \u03B1,2 =", balpha2/nat,u' eV/atom')
          print(" ")
          print("b \u03B1,2 =", balpha2/vol,u' eV/A^3')
          print(" ")
          print("b \u03B1,2 =", (balpha2/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
          print(" ")



          if self.args.delas == True:
            print(" ")
            print(" ")
            print("----------------------------------------------")
            print("Calculation of isotropic magnetostrictive coefficient:")
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


            bbalpha2 = (balpha2/vol)*1.602176565e-19*1e30*1e-9  #GPa

            
            lambdaalpha2 = -(bbalpha2)/(c11+2*c12)
            

            


            print("c11 =", str(c11), 'GPa')
            print(" ")
            print("c12 =", str(c12), 'GPa')
            print(" ")
            print("c44 =", str(c44), 'GPa')
            print(" ")
            print("Warning: If these elastic constants are not the same as in the input elastic tensor file", str(self.args.elas[0]),", then check that the format of the elastic tensor is exactly the same as in the standard output file ELADAT generated by AELAS code (see Example folder)")
            print(" ")
            print(" ")
            print("Isotropic Magnetostrictive coefficient:")
            print(" ")
            print("Using the convention in reference E. D. T. de Lacheisserie, Magnetostriction: Theory and Application of Magnetoelasticity (CRC Press, Boca Raton, FL, 1993).")
            print(" ")
            print("\u03BB\u03B1,2 =", lambdaalpha2*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("Using the convention in reference J.R. Cullen et al., in Materials, Science and Technology (VCH Publishings, 1994), pp.529-565:")
            print(" ")
            print("\u03BB\u03B1 =", (lambdaalpha2/3.0)*1e6,u'x 10\u207B\u2076')
            print(" ")
            print(" ------------------------------ ")
            print(" ")
            print("Spontaneous volume magnetostriction:")
            print(" ")
            print("\u03C9_s =", lambdaalpha2*1e6,u'x 10\u207B\u2076')
            print(" ")


        
        # Hexagonal I, Trigonal I or Tetragonal I:
        
        elif (177 <= sg <= 194) or (149 <= sg <= 167) or (89 <= sg <= 142):


          print("")
          print("Fit of quadractic function f(x)=A*x^2+B*x+C to energy vs strain")
          print("")
          print("-------------------------")
          print("Calculation of b_11:")
          print("-------------------------")
          print(" ")
          print('Deformation x-direction (\u03B5xx) and y-direction (\u03B5yy) , spin1=[0,0,1]')
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

          print("Fitting parameters:")
          print("A =", params[0][0], ", B =", params[0][1], ", C =", params[0][2])

          r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
          print("R-squared =", r_squared)
          print("")

          if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_1_1.png")
            print("")

          b11 = (1.0/2.0)*params[0][1]



        #make figure dE_1.png

          fig = 'fit_ene_1_1.png'
          tit = "Calculation of b_11"

          plt.plot(x, y, 'o', label='data in ene_1_1.dat')
          t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
          plt.plot(t, K(t, params[0][0],params[0][1],params[0][2]), 'r-', label='fit')
          plt.legend()
          ax = plt.gca()
          ax.xaxis.set_major_locator(plt.MaxNLocator(5))

          ylabel ='Energy (eV)'
          plt.ylabel(ylabel)
          label = "Strain \u03B5xx , \u03B5yy"
          plt.xlabel(label)
          plt.title(tit)
          plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
          plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
          plt.savefig(fig)
          plt.close()


          print("")
          print("-------------------------")
          print("Calculation of b_12:")
          print("-------------------------")
          print(" ")
          print('Deformation z-direction (\u03B5zz), spin1=[0,0,1]')
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

          print("Fitting parameters:")
          print("A =", params[0][0], ", B =", params[0][1], ", C =", params[0][2])

          r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
          print("R-squared =", r_squared)
          print("")

          if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_2_1.png")
            print("")

          b12 = params[0][1]



        #make figure dE_2.png

          fig = 'fit_ene_2_1.png'
          tit = "Calculation of b_12"

          plt.plot(x, y, 'o', label='data in ene_2_1.dat')
          t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
          plt.plot(t, K(t, params[0][0],params[0][1],params[0][2]), 'r-', label='fit')
          plt.legend()
          ax = plt.gca()
          ax.xaxis.set_major_locator(plt.MaxNLocator(5))

          ylabel ='Energy (eV)'
          plt.ylabel(ylabel)
          label = "Strain \u03B5zz"
          plt.xlabel(label)
          plt.title(tit)
          plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
          plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
          plt.savefig(fig)
          plt.close()



          print(" ")
          print("----------------------------------------------")
          print("Isotropic magnetoelastic constant:")
          print("----------------------------------------------")
          print(" ")
          print(" ")
          print("Using the convention in reference P. Nieves et al., Comput. Phys. Commun. 264, 107964 (2021):")
          print(" ")
          print("b_11 =", b11,u' eV')
          print(" ")
          print("b_12 =", b12,u' eV')
          print(" ")
          print("b_11 =", b11/nat,u' eV/atom')
          print(" ")
          print("b_12 =", b12/nat,u' eV/atom')
          print(" ")
          print("b_11 =", b11/vol,u' eV/A^3')
          print(" ")
          print("b_12 =", b12/vol,u' eV/A^3')
          print(" ")
          print("b_11 =", (b11/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
          print(" ")
          print("b_12 =", (b12/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
          print(" ")
          print("Using the convention in reference E. D. T. de Lacheisserie, Magnetostriction: Theory and Application of Magnetoelasticity (CRC Press, Boca Raton, FL, 1993):")
          b1a0 = 2.0*b11+b12
          b2a0 = -math.sqrt(2)*(b11-b12)
          print(" ")
          print("b 1\u03B1,0 =", (b1a0/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
          print(" ")
          print("b 2\u03B1,0 =", (b2a0/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
          print(" ")


          if self.args.delas == True:
            print(" ")
            print(" ")
            print("----------------------------------------------")
            print("Calculation of isotropic magnetostrictive coefficient:")
            print("----------------------------------------------")
            print(" ")
            print("Reading the elastic tensor file =", str(self.args.elas[0]))
            print(" ")

            if 177 <= sg <= 194:
            
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


              print("c11 =", str(c11), 'GPa')
              print(" ")
              print("c12 =", str(c12), 'GPa')
              print(" ")
              print("c13 =", str(c13), 'GPa')
              print(" ")
              print("c33 =", str(c33), 'GPa')
              print(" ")
              print("c44 =", str(c44), 'GPa')

            
            elif 149 <= sg <= 167:
            
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

            
            
            elif 89 <= sg <= 142:
            

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





            bb11 = (b11/vol)*1.602176565e-19*1e30*1e-9  #GPa
            bb12 = (b12/vol)*1.602176565e-19*1e30*1e-9  #GPa

            
            lambdaalpha10 = (bb11*c33+bb12*c13)/(c33*(c11+c12)-2.0*c13**2.0)
            
            lambdaalpha20 = (2.0*bb11*c13-bb12*(c11+c12))/(c33*(c11+c12)-2.0*c13**2.0)   


            print(" ")
            print("Warning: If these elastic constants are not the same as in the input elastic tensor file", str(self.args.elas[0]),", then check that the format of the elastic tensor is exactly the same as in the standard output file ELADAT generated by AELAS code (see Example folder)")
            print(" ")
            print(" ")
            print("Isotropic Magnetostrictive coefficients:")
            print(" ")
            print("Using the convention in reference P. Nieves et al., Comput. Phys. Commun. 264, 107964 (2021):")
            print(" ")
            print("\u03BB\u03B1 1,0 =", lambdaalpha10*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB\u03B1 2,0 =", lambdaalpha20*1e6,u'x 10\u207B\u2076')
            print(" ")
            if 177 <= sg <= 194:
              print("Using the convention in reference E. D. T. de Lacheisserie, Magnetostriction: Theory and Application of Magnetoelasticity (CRC Press, Boca Raton, FL, 1993):")
              print(" ")
              print("\u03BB 1\u03B1,0 =", (2.0*lambdaalpha10+lambdaalpha20)*1e6,u'x 10\u207B\u2076')
              print(" ")
              print("\u03BB 2\u03B1,0 =", (-lambdaalpha10+lambdaalpha20)*1e6,u'x 10\u207B\u2076')
              print(" ")
            elif 89 <= sg <= 142:
              print("Using the convention in reference E. D. T. de Lacheisserie, Magnetostriction: Theory and Application of Magnetoelasticity (CRC Press, Boca Raton, FL, 1993):")
              print(" ")
              print("\u03BB 1\u03B1,0 =", (2.0/3.0)*(2.0*lambdaalpha10+lambdaalpha20)*1e6,u'x 10\u207B\u2076')
              print(" ")
              print("\u03BB 2\u03B1,0 =", (2.0/3.0)*(-lambdaalpha10+lambdaalpha20)*1e6,u'x 10\u207B\u2076')
              print(" ")
              
            print(" ------------------------------ ")
            print(" ")
            print("Spontaneous volume magnetostriction:")
            print(" ")
            print("\u03C9_s =", (2*lambdaalpha10+lambdaalpha20)*1e6,u'x 10\u207B\u2076')
            print(" ")



        # Orthorhombic:
        
        elif 16 <= sg <= 74:


          print("")
          print("Fit of quadractic function f(x)=A*x^2+B*x+C to energy vs strain")
          print("")
          print("-------------------------")
          print("Calculation of b_1\u03B1,0:")
          print("-------------------------")
          print(" ")
          print('Deformation \u03B5xx=s,\u03B5yy=s,\u03B5zz=s:')
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

          print("Fitting parameters:")
          print("A =", params[0][0], ", B =", params[0][1], ", C =", params[0][2])

          r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
          print("R-squared =", r_squared)
          print("")

          if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_1_1.png")
            print("")

          b01 = params[0][1]



        #make figure dE_1.png

          fig = 'fit_ene_1_1.png'
          tit = "Calculation of b_1\u03B1,0"

          plt.plot(x, y, 'o', label='data in ene_1_1.dat')
          t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
          plt.plot(t, K(t, params[0][0],params[0][1],params[0][2]), 'r-', label='fit')
          plt.legend()
          ax = plt.gca()
          ax.xaxis.set_major_locator(plt.MaxNLocator(5))

          ylabel ='Energy (eV)'
          plt.ylabel(ylabel)
          label = "Strain \u03B5xx=s,\u03B5yy=s,\u03B5zz=s"
          plt.xlabel(label)
          plt.title(tit)
          plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
          plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
          plt.savefig(fig)
          plt.close()


          print("")
          print("-------------------------")
          print("Calculation of b_2\u03B1,0:")
          print("-------------------------")
          print(" ")
          print('Deformation \u03B5xx=-s/2,\u03B5yy=-s/2,\u03B5zz=s')
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

          print("Fitting parameters:")
          print("A =", params[0][0], ", B =", params[0][1], ", C =", params[0][2])

          r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
          print("R-squared =", r_squared)
          print("")

          if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_2_1.png")
            print("")

          b02 = math.sqrt(2)*params[0][1]



        #make figure dE_2.png

          fig = 'fit_ene_2_1.png'
          tit = "Calculation of b 2\u03B1,0"

          plt.plot(x, y, 'o', label='data in ene_2_1.dat')
          t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
          plt.plot(t, K(t, params[0][0],params[0][1],params[0][2]), 'r-', label='fit')
          plt.legend()
          ax = plt.gca()
          ax.xaxis.set_major_locator(plt.MaxNLocator(5))

          ylabel ='Energy (eV)'
          plt.ylabel(ylabel)
          label = "Strain \u03B5xx=-s/2,\u03B5yy=-s/2,\u03B5zz=s"
          plt.xlabel(label)
          plt.title(tit)
          plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
          plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
          plt.savefig(fig)
          plt.close()


          print("")
          print("-------------------------")
          print("Calculation of b_3\u03B1,0:")
          print("-------------------------")
          print(" ")
          print('Deformation \u03B5xx=s,\u03B5yy=-s')
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

          print("Fitting parameters:")
          print("A =", params[0][0], ", B =", params[0][1], ", C =", params[0][2])

          r_squared = r2_score(y, K(x,params[0][0],params[0][1],params[0][2]))
          print("R-squared =", r_squared)
          print("")

          if r_squared < 0.98:
            print("WARNING!! R-squared is lower than 0.98. Check figure fit_ene_3_1.png")
            print("")

          b03 = -math.sqrt(3.0/2.0)*params[0][1]



        #make figure dE_3.png

          fig = 'fit_ene_3_1.png'
          tit = "Calculation of b 3\u03B1,0"

          plt.plot(x, y, 'o', label='data in ene_3_1.dat')
          t = np.arange(min(x)-0.05*(max(x)-min(x)), max(x)+0.05*(max(x)-min(x)), 0.0001)
          plt.plot(t, K(t, params[0][0],params[0][1],params[0][2]), 'r-', label='fit')
          plt.legend()
          ax = plt.gca()
          ax.xaxis.set_major_locator(plt.MaxNLocator(5))

          ylabel ='Energy (eV)'
          plt.ylabel(ylabel)
          label = "Strain \u03B5xx=s,\u03B5yy=-s"
          plt.xlabel(label)
          plt.title(tit)
          plt.tight_layout(pad=6, h_pad=None, w_pad=None, rect=None)
          plt.ticklabel_format(axis='both', style='plain', useOffset=False, useMathText=True)
          plt.savefig(fig)
          plt.close()




          print(" ")
          print("----------------------------------------------")
          print("Isotropic magnetoelastic constant:")
          print("----------------------------------------------")
          print(" ")
          print(" ")
          print("Using the convention in reference E. D. T. de Lacheisserie, Magnetostriction: Theory and Application of Magnetoelasticity (CRC Press, Boca Raton, FL, 1993):")
          print(" ")
          print("b_1\u03B1,0 =", b01,u' eV')
          print(" ")
          print("b_2\u03B1,0 =", b02,u' eV')
          print(" ")
          print("b_3\u03B1,0 =", b03,u' eV')
          print(" ")
          print("b_1\u03B1,0 =", b01/nat,u' eV/atom')
          print(" ")
          print("b_2\u03B1,0 =", b02/nat,u' eV/atom')
          print(" ")
          print("b_3\u03B1,0 =", b03/nat,u' eV/atom')
          print(" ")
          print("b_1\u03B1,0 =", b01/vol,u' eV/A^3')
          print(" ")
          print("b_2\u03B1,0 =", b02/vol,u' eV/A^3')
          print(" ")
          print("b_3\u03B1,0 =", b03/vol,u' eV/A^3')
          print(" ")
          print("b_1\u03B1,0 =", (b01/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
          print(" ")
          print("b_2\u03B1,0 =", (b02/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
          print(" ")
          print("b_3\u03B1,0 =", (b03/vol)*1.602176565e-19*1e30*1e-6 ,u' MPa')
          print(" ")



          if self.args.delas == True:
            print(" ")
            print(" ")
            print("----------------------------------------------")
            print("Calculation of isotropic magnetostrictive coefficient:")
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


            b1a0 = (b01/vol)*1.602176565e-19*1e30*1e-9  #GPa
            b2a0 = (b02/vol)*1.602176565e-19*1e30*1e-9  #GPa
            b3a0 = (b03/vol)*1.602176565e-19*1e30*1e-9  #GPa


            ca11 = (1.0/3.0)*(c11+2.0*c12+2.0*c13+c22+2.0*c23+c33)
            ca22 = (1.0/6.0)*(c11+2.0*c12-4.0*c13+c22-4.0*c23+4.0*c33)
            ca33 = (1.0/2.0)*(c11-2.0*c12+c22)
            ca12 = (-c11-2.0*c12+c13-c22+c23+2.0*c33)/(3.0*math.sqrt(2))
            ca13 = (c11+c13-c22-c23)/(math.sqrt(6))
            ca23 = (-c11+2.0*c13+c22-2.0*c23)/(2.0*math.sqrt(3))

            lmb1a0 = (-b3a0*ca13*ca22+b3a0*ca12*ca23+b2a0*ca13*ca23-b1a0*ca23**2-b2a0*ca12*ca33+b1a0*ca22*ca33)/(ca13**2*ca22-2.0*ca12*ca13*ca23+ca12**2*ca33+ca11*(ca23**2-ca22*ca33))
            
            lmb2a0 = (b3a0*ca12*ca13-b2a0*ca13**2-b3a0*ca11*ca23+b1a0*ca13*ca23+b2a0*ca11*ca33-b1a0*ca12*ca33)/(math.sqrt(2)*(ca13**2*ca22-2.0*ca12*ca13*ca23+ca12**2*ca33+ca11*(ca23**2-ca22*ca33)))
            
            lmb3a0 = (-b3a0*ca12**2+b2a0*ca12*ca13+b3a0*ca11*ca22-b1a0*ca13*ca22-b2a0*ca11*ca23+b1a0*ca12*ca23)/(math.sqrt(6)*(ca13**2*ca22-2.0*ca12*ca13*ca23+ca12**2*ca33+ca11*(ca23**2-ca22*ca33)))  



            print("Isotropic Magnetostrictive coefficients:")
            print(" ")
            print("Using the convention in reference E. D. T. de Lacheisserie, Magnetostriction: Theory and Application of Magnetoelasticity (CRC Press, Boca Raton, FL, 1993):")
            print(" ")
            print("\u03BB 1\u03B1,0 =", lmb1a0*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB 2\u03B1,0 =", lmb2a0*1e6,u'x 10\u207B\u2076')
            print(" ")
            print("\u03BB 3\u03B1,0 =", lmb3a0*1e6,u'x 10\u207B\u2076')
            print(" ")
            print(" ------------------------------ ")
            print(" ")
            print("Spontaneous volume magnetostriction:")
            print(" ")
            print("\u03C9_s =", (lmb1a0+(1.0/3.0)*lmb2a0)*1e6,u'x 10\u207B\u2076')
            print(" ")



        
        
        