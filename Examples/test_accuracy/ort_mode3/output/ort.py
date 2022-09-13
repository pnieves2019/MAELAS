#!/bin/bash

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import stat



c11=76.0*10**9
c12=45.0*10**9
c13=48.0*10**9
c22=102.0*10**9
c23=55.0*10**9
c33=141.0*10**9
c44=40.0*10**9
c55=27.0*10**9
c66=39.0*10**9


b01=0.3*10**9
b02=-0.2*10**9
b03=0.7*10**9
b1=0*10**6
b2=0*10**6
b3=0*10**6
b4=0*10**6
b5=0*10**6
b6=0*10**6
b7=0*10**6
b8=0*10**6
b9=0*10**6


a=4.0686275541928172 
b=10.3156917156588115
c=3.8956767252497921
v0=(a*b*c)*((10.0**(-10.0))**3.0)

# Tremolet convention of isotropic magentoelastic constants

b1a0 = b01  #GPa
b2a0 = b02  #GPa
b3a0 = b03  #GPa


ca11 = (1.0/3.0)*(c11+2.0*c12+2.0*c13+c22+2.0*c23+c33)
ca22 = (1.0/6.0)*(c11+2.0*c12-4.0*c13+c22-4.0*c23+4.0*c33)
ca33 = (1.0/2.0)*(c11-2.0*c12+c22)
ca12 = (-c11-2.0*c12+c13-c22+c23+2.0*c33)/(3.0*math.sqrt(2))
ca13 = (c11+c13-c22-c23)/(math.sqrt(6))
ca23 = (-c11+2.0*c13+c22-2.0*c23)/(2.0*math.sqrt(3))

lmb1a0 = (-b3a0*ca13*ca22+b3a0*ca12*ca23+b2a0*ca13*ca23-b1a0*ca23**2-b2a0*ca12*ca33+b1a0*ca22*ca33)/(ca13**2*ca22-2.0*ca12*ca13*ca23+ca12**2*ca33+ca11*(ca23**2-ca22*ca33))
            
lmb2a0 = (b3a0*ca12*ca13-b2a0*ca13**2-b3a0*ca11*ca23+b1a0*ca13*ca23+b2a0*ca11*ca33-b1a0*ca12*ca33)/(math.sqrt(2)*(ca13**2*ca22-2.0*ca12*ca13*ca23+ca12**2*ca33+ca11*(ca23**2-ca22*ca33)))
            
lmb3a0 = (-b3a0*ca12**2+b2a0*ca12*ca13+b3a0*ca11*ca22-b1a0*ca13*ca22-b2a0*ca11*ca23+b1a0*ca12*ca23)/(math.sqrt(6)*(ca13**2*ca22-2.0*ca12*ca13*ca23+ca12**2*ca33+ca11*(ca23**2-ca22*ca33)))  



w_s = lmb1a0

ndist = 7
strain = 0.01


def elas0(eexx,eexy,eexz,eeyy,eeyz,eezz,cc11,cc12,cc13,cc22,cc23,cc33,cc44,cc55,cc66):
    
    ene=0.5*cc11*eexx**2.0+0.5*cc22*eeyy**2.0+cc12*eexx*eeyy
    ene=ene+cc13*eexx*eezz+cc23*eeyy*eezz+0.5*cc33*eezz**2.0
    ene=ene+2.0*cc44*eeyz**2.0+2.0*cc55*eexz**2.0+2.0*cc66*eexy**2.0
    
    return ene

# Tremolet convention for the isotropic magnetoelastic energy

def em0(eexx,eexy,eexz,eeyy,eeyz,eezz,bb01,bb02,bb03,bb1,bb2,bb3,bb4,bb5,bb6,bb7,bb8,bb9,aax,aay,aaz):
    
#    ene0 = bb01*eexx + bb02*eeyy + bb03*eezz
#    ene0 = ene0 + bb1*aax**2.0*eexx + bb2*aay**2.0*eexx + bb3*aax**2.0*eeyy + bb4*aay**2.0*eeyy
#    ene0 = ene0 + bb5*aax**2.0*eezz + bb6*aay**2.0*eezz
#    ene0 = ene0 + 2.0*bb7*aax*aay*eexy + 2.0*bb8*aax*aaz*eexz + 2.0*bb9*aay*aaz*eeyz

    ene0 = (1/3)*bb01*(eexx+eeyy+eezz)+(math.sqrt(2)/3)*bb02*(eezz-0.5*(eexx+eeyy))+(1/math.sqrt(6))*bb03*(eexx-eeyy)

    return ene0






for i in range(int(ndist)):


    strain1 = - float(strain)+2.0*(float(strain)/(float(ndist)-1.0))*i

    print("strain", strain1)


        #Generation POSCAR file

        #b01

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
    
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c22,c23,c33,c44,c55,c66)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b01,b02,b03,b1,b2,b3,b4,b5,b6,b7,b8,b9,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_1_" + str(i+1) + "_1"
    
    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    
    
    
    
    

    
############ b02
    
    
    fxx = 1.0 - 0.5*strain1
    fyy = 1.0 - 0.5*strain1
    fzz = 1.0 + strain1
    fxy = 0.0
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
    ay=1.0
    az=0.0
    
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c22,c23,c33,c44,c55,c66)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b01,b02,b03,b1,b2,b3,b4,b5,b6,b7,b8,b9,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_2_" + str(i+1) + "_1"
    
    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    
    
    
############ b03
    
    
    fyy = 1.0 + strain1
    fxx = 1.0 - strain1
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
    
    ax=1.0
    ay=0.0
    az=0.0
    
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c22,c23,c33,c44,c55,c66)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b01,b02,b03,b1,b2,b3,b4,b5,b6,b7,b8,b9,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_3_" + str(i+1) + "_1"
    
    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    
    

lmb_name = "output_exact.dat"
    
dat = open(lmb_name,'w')
dat.write("lmb_1alpha0 = ")
dat.write(repr(lmb1a0))
dat.write('\n')
dat.write("lmb_2alpha0= ")
dat.write(repr(lmb2a0))
dat.write('\n')
dat.write("lmb_3alpha0= ")
dat.write(repr(lmb3a0))
dat.write('\n')


dat.write("b_1alpha0 = ")
dat.write(repr(b01*1e-6))
dat.write(" MPa")
dat.write('\n')
dat.write("b_2alpha0= ")
dat.write(repr(b02*1e-6))
dat.write(" MPa")
dat.write('\n')
dat.write("b_3alpha0= ")
dat.write(repr(b03*1e-6))
dat.write(" MPa")
dat.write('\n')
dat.write("w_s= ")
dat.write(repr(w_s))
dat.write('\n')



dat.close()
    
    

