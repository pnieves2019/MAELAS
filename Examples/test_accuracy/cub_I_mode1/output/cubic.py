#!/bin/bash

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import stat



c11=243.0*10**9  # same values as in ELADAT!
c12=138.0*10**9
c44=122.0*10**9
b0=0.0
b1=-4.1*10**6
b2=10.9*10**6
a=2.8293474424254716
v0=(a*10**(-10))**3.0  # volume poscar!

lmb001=-2*b1/(3*(c11-c12))
lmb111=-b2/(3*c44)

ndist = 7
strain = 0.01


for i in range(int(ndist)):


    strain1 = - float(strain)+2.0*(float(strain)/(float(ndist)-1.0))*i

    print("strain", strain1)


        #Generation POSCAR file

        #lambda_001

    fzz = 1.0 + strain1
    fxx = 1.0/math.sqrt(fzz)
    fyy = fxx
    fxy=0.0
    fxz = fxy
    fyx = fxy
    fyz = fxy
    fzx = fxy
    fzy = fxy
    
    exx=fxx-1.0
    exy=0.5*(fxy+fyx)
    exz=0.5*(fxz+fzx)
    eyy=fyy-1.0
    eyz=0.5*(fyz+fzy)
    ezz=fzz-1.0
    
    ax=0.0
    ay=0.0
    az=1.0
    
    elas=0.5*c11*(exx**2+eyy**2+ezz**2)+c12*(exx*eyy+exx*ezz+eyy*ezz)+2*c44*(exy**2+exz**2+eyz**2)
    em=b0*(exx+eyy+ezz)+b1*(ax**2*exx+ay**2*eyy+az**2*ezz)+2*b2*(ax*ay*exy+ax*az*exz+ay*az*eyz)
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_1_" + str(i+1) + "_1"
    
    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write(repr(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    
    
    
    
    ax=1.0
    ay=0.0
    az=0.0
     
    elas=0.5*c11*(exx**2+eyy**2+ezz**2)+c12*(exx*eyy+exx*ezz+eyy*ezz)+2*c44*(exy**2+exz**2+eyz**2)
    em=b0*(exx+eyy+ezz)+b1*(ax**2*exx+ay**2*eyy+az**2*ezz)+2*b2*(ax*ay*exy+ax*az*exz+ay*az*eyz)
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_1_" + str(i+1) + "_2"


    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write(repr(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    

    
############ lmb 111
    
    
    const = (4/(4-3*(strain1**2)+strain1**3))**(1/3)
    
    fzz = const
    fxx = const
    fyy = const
    fxy=const*strain1*0.5
    fxz = fxy
    fyx = fxy
    fyz = fxy
    fzx = fxy
    fzy = fxy
    
    exx=fxx-1.0
    exy=0.5*(fxy+fyx)
    exz=0.5*(fxz+fzx)
    eyy=fyy-1.0
    eyz=0.5*(fyz+fzy)
    ezz=fzz-1.0
    
    ax=1.0/math.sqrt(3)
    ay=1.0/math.sqrt(3)
    az=1.0/math.sqrt(3)
    
    elas=0.5*c11*(exx**2+eyy**2+ezz**2)+c12*(exx*eyy+exx*ezz+eyy*ezz)+2*c44*(exy**2+exz**2+eyz**2)
    em=b0*(exx+eyy+ezz)+b1*(ax**2*exx+ay**2*eyy+az**2*ezz)+2*b2*(ax*ay*exy+ax*az*exz+ay*az*eyz)
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_2_" + str(i+1) + "_1"
    
    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write(repr(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    
    
    
    
    ax=1.0/math.sqrt(2)
    ay=0.0
    az=-1.0/math.sqrt(2)
     
    elas=0.5*c11*(exx**2+eyy**2+ezz**2)+c12*(exx*eyy+exx*ezz+eyy*ezz)+2*c44*(exy**2+exz**2+eyz**2)
    em=b0*(exx+eyy+ezz)+b1*(ax**2*exx+ay**2*eyy+az**2*ezz)+2*b2*(ax*ay*exy+ax*az*exz+ay*az*eyz)
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_2_" + str(i+1) + "_2"


    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write(repr(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()


lmb_name = "output_exact.dat"

b1 = b1*10**-9
b2 = b2*10**-9

    
dat = open(lmb_name,'w')
dat.write("lmb001 = ")
dat.write(repr(lmb001))
dat.write('\n')
dat.write("lmb111 = ")
dat.write(repr(lmb111))
dat.write('\n')
dat.write("b1 = ")
dat.write(repr(b1))
dat.write(" GPa")
dat.write('\n')
dat.write("b2 = ")
dat.write(repr(b2))
dat.write(" GPa")
dat.write('\n')

dat.close()
    

