#!/bin/bash

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import stat



c11=428.0*10**9
c12=164.0*10**9
c13=133.0*10**9
c14=-27.0*10**9
c33=434.0*10**9
c44=118*10**9
b11=20.0*10**6
b12=-50.0*10**6
b21=43.1*10**6
b22=(-34.2)*10**6
b3=60.7*10**6
b4=(-34.3)*10**6
b14=(-42.4)*10**6
b34=55.4*10**6


a=3.9248957323000999
c=4.8311031919453633
v0=(((math.sqrt(3.0)*a**2.0)/2.0)*c)*((10.0**(-10.0))**3.0)

#lmb_alpha10=(b11*c33+b12*c13)/(c33*(c11+c12)-2.0*c13**2.0)

lmb_alpha10=(-b11*c33+b12*c13)/(c33*(c11+c12)-2.0*c13**2.0)
lmb_alpha20=(2*b11*c13-b12*(c11+c12))/(c33*(c11+c12)-2.0*c13**2.0)
lmba12=(-b21*c33+b22*c13)/(c33*(c11+c12)-2.0*c13**2.0)
lmba22=(2.0*b21*c13-b22*(c11+c12))/(c33*(c11+c12)-2.0*c13**2.0)
lmbg1 = 2.0*(0.25*c14*b14-0.25*c44*b3)/(0.5*c44*(c11-c12)-c14**2.0)
lmbg2 = 2.0*(-0.25*b4*(c11-c12)+0.5*b34*c14)/(0.5*c44*(c11-c12)-c14**2.0)              
lmb12 = 2.0*(0.5*c14*b4-0.5*c44*b34)/(0.5*c44*(c11-c12)-c14**2.0)               
lmb21 = 2.0*(-0.25*b14*(c11-c12)+0.5*b3*c14)/(0.5*c44*(c11-c12)-c14**2.0)


w_s=2*lmb_alpha10+lmb_alpha20

ndist = 7
strain = 0.01


    
def elas0(eexx,eexy,eexz,eeyy,eeyz,eezz,cc11,cc12,cc13,cc14,cc33,cc44):
    
    ene=0.5*cc11*(eexx**2+eeyy**2)+cc12*(eexx*eeyy)+cc13*(eexx+eeyy)*eezz+0.5*cc33*(eezz**2)+2.0*cc44*(eexz**2.0+eeyz**2.0)+(cc11-cc12)*eexy**2.0+cc14*(4.0*eexy*eexz+2.0*eexx*eeyz-2.0*eeyy*eeyz)
    
    return ene

    
def em0(eexx,eexy,eexz,eeyy,eeyz,eezz,bb11,bb12,bb21,bb22,bb3,bb4,bb14,bb34,aax,aay,aaz):
    
    ene0 = bb11*(eexx+eeyy)+bb12*eezz
    ene0 = ene0 + bb21*(aaz**2.0-(1.0/3.0))*(eexx+eeyy) + bb22*(aaz**2.0-(1.0/3.0))*eezz
    ene0 = ene0 + bb3*(0.5*(aax**2.0-aay**2.0)*(eexx-eeyy)+2.0*aax*aay*eexy)
    ene0 = ene0 + 2.0*bb4*(aax*aaz*eexz+aay*aaz*eeyz)
    ene0 = ene0 + bb14*((aax**2-aay**2)*eeyz+2.0*aax*aay*eexz)
    ene0 = ene0 + bb34*(0.5*aay*aaz*(eexx-eeyy)+2.0*aax*aaz*eexy)
    
    return ene0




for i in range(int(ndist)):


    strain1 = - float(strain)+2.0*(float(strain)/(float(ndist)-1.0))*i

    print("strain", strain1)


        #Generation POSCAR file

        # b_11

    fxx = 1.0 + strain1
    fyy = 1.0 + strain1
    fzz = 1.0
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
    
    ax=0
    ay=0
    az=1.0
    
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c14,c33,c44)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b11,b12,b21,b22,b3,b4,b14,b34,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_1_" + str(i+1) + "_1"
    
    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    

    

    
############ b_12
    
    
    fzz = 1.0 + strain1
    fxx = 1.0
    fyy = 1.0
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
    
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c14,c33,c44)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b11,b12,b21,b22,b3,b4,b14,b34,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_2_" + str(i+1) + "_1"
    
    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    



lmb_name = "output_exact.dat"
    
dat = open(lmb_name,'w')
dat.write("lmb_alpha10 = ")
dat.write(repr(lmb_alpha10))
dat.write('\n')
dat.write("lmb_alpha20= ")
dat.write(repr(lmb_alpha20))
dat.write('\n')


dat.write("b11 = ")
dat.write(repr(b11*1e-6))
dat.write(" MPa")
dat.write('\n')
dat.write("b12= ")
dat.write(repr(b12*1e-6))
dat.write(" MPa")
dat.write('\n')
dat.write("w_s= ")
dat.write(repr(w_s))
dat.write('\n')



dat.close()
    
    

