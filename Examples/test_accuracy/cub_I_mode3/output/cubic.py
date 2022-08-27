#!/bin/bash

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import stat



c11=243.0*10**9
c12=138.0*10**9
c44=122.0*10**9
balpha=-1.0*10**9
b1=0*10**6
b2=0*10**6
a=2.8293474424254716
v0=(a*10**(-10))**3.0

lmbalpha2=-balpha/(c11+2*c12)
lmbalpha=lmbalpha2/3.0
w_s=lmbalpha2

ndist = 7
strain = 0.01


for i in range(int(ndist)):


    strain1 = - float(strain)+2.0*(float(strain)/(float(ndist)-1.0))*i

    print("strain", strain1)


        #Generation POSCAR file

        #b0

    fxx = 1.0 + strain1
    fyy = 1.0 + strain1
    fzz = 1.0 + strain1
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
    
    ax=1.0
    ay=0.0
    az=0.0
    
    elas=0.5*c11*(exx**2+eyy**2+ezz**2)+c12*(exx*eyy+exx*ezz+eyy*ezz)+2*c44*(exy**2+exz**2+eyz**2)
    em=(1/3)*balpha*(exx+eyy+ezz)+b1*((ax**2-(1/3))*exx+(ay**2-(1/3))*eyy+(az**2-(1/3))*ezz)+2*b2*(ax*ay*exy+ax*az*exz+ay*az*eyz)
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_1_" + str(i+1) + "_1"
    
    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write(repr(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    

lmb_name = "output_exact.dat"
    
dat = open(lmb_name,'w')
dat.write("b alpha,2 = ")
dat.write(repr(balpha*1e-9))
dat.write(" GPa")
dat.write('\n')
dat.write("lmbalpha,2 = ")
dat.write(repr(lmbalpha2))
dat.write('\n')
dat.write("lmbalpha = ")
dat.write(repr(lmbalpha))
dat.write('\n')
dat.write("w_s = ")
dat.write(repr(w_s))
dat.write('\n')


dat.close()
    

