#!/bin/bash

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import stat



c11=307.0*10**9
c12=165.0*10**9
c13=103.0*10**9
c33=358.0*10**9
c44=75.0*10**9
b21=-31.9*10**6
b22=25.5*10**6
b3=-8.1*10**6
b4=42.9*10**6


#c11=150.0*10**9
#c12=80.0*10**9
#c13=50.0*10**9
#c33=100.0*10**9
#c44=100.0*10**9
#b21=-150.0*10**6
#b22=100.0*10**6
#b3=-150.0*10**6
#b4=100.0*10**6



a=2.4560768153817771
c=3.9821032982002493
v0=(((math.sqrt(3.0)*a**2.0)/2.0)*c)*((10.0**(-10.0))**3.0)

lmba12=(-b21*c33+b22*c13)/(c33*(c11+c12)-2.0*c13**2.0)
lmba22=(2.0*b21*c13-b22*(c11+c12))/(c33*(c11+c12)-2.0*c13**2.0)
lmbg2=-b3/(c11-c12)
lmbe2=-b4/(2.0*c44)

ndist = 7
strain = 0.01


def elas0(eexx,eexy,eexz,eeyy,eeyz,eezz,cc11,cc12,cc13,cc33,cc44):
    
    ene=0.5*cc11*(eexx**2+eeyy**2)+cc12*(eexx*eeyy)+cc13*(eexx+eeyy)*eezz+0.5*cc33*(eezz**2)+2.0*cc44*(eexz**2.0+eeyz**2.0)+(cc11-cc12)*eexy**2.0
 
    return ene

def em0(eexx,eexy,eexz,eeyy,eeyz,eezz,bb21,bb22,bb3,bb4,aax,aay,aaz):
    
    ene0 = bb21*(aaz**2.0-(1.0/3.0))*(eexx+eeyy) + bb22*(aaz**2.0-(1.0/3.0))*eezz + bb3*(0.5*(aax**2.0-aay**2.0)*(eexx-eeyy)+2.0*aax*aay*eexy)
    ene0 = ene0 + 2.0*bb4*(aax*aaz*eexz+aay*aaz*eeyz)
    
    return ene0


#def elas0(eexx,eexy,eexz,eeyy,eeyz,eezz,cc11,cc12,cc13,cc33,cc44):
    
#    ene=0.5*cc11*(eexx**2+eeyy**2)+cc12*(eexx*eeyy)+cc13*(eexx+eeyy)*eezz+0.5*cc33*(eezz**2)+0.5*cc44*(eexz**2+eeyz**2)+0.25*(cc11-cc12)*eexy**2
 
#    return ene

#def em0(eexx,eexy,eexz,eeyy,eeyz,eezz,bb21,bb22,bb3,bb4,aax,aay,aaz):
    
#    ene0 = bb21*(aaz**2-(1/3))*(eexx+eeyy) + bb22*(aaz**2-(1/3))*eezz + bb3*(0.5*(aax**2-aay**2)*(eexx-eeyy)+aax*aay*eexy)
#    ene0 = ene0 + bb4*(aax*aaz*eexz+aay*aaz*eeyz)
    
#    return ene0



for i in range(int(ndist)):


    strain1 = - float(strain)+2.0*(float(strain)/(float(ndist)-1.0))*i

    print("strain", strain1)


        #Generation POSCAR file

        #lambda_alpha12

    fxx = 1.0 + 0.5*strain1
    fyy = 1.0 + 0.5*strain1
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
    
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c33,c44)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b21,b22,b3,b4,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_1_" + str(i+1) + "_1"
    
    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    
    
    
    
    ax=1.0/math.sqrt(2)
    ay=1.0/math.sqrt(2)
    az=0.0
     
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c33,c44)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b21,b22,b3,b4,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_1_" + str(i+1) + "_2"


    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    

    
############ lmb_alpha_22
    
    
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
    
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c33,c44)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b21,b22,b3,b4,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_2_" + str(i+1) + "_1"
    
    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    
    
    
    
    ax=1.0/math.sqrt(2)
    ay=1.0/math.sqrt(2)
    az=0.0
     
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c33,c44)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b21,b22,b3,b4,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_2_" + str(i+1) + "_2"


    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    
    
    
############ lmb_gamma_2
    
    
    fxx = 1.0 + 0.5*strain1
    fyy = 1.0 - 0.5*strain1
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
    
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c33,c44)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b21,b22,b3,b4,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_3_" + str(i+1) + "_1"
    
    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    
    
    ax=0.0
    ay=1.0
    az=0.0
     
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c33,c44)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b21,b22,b3,b4,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_3_" + str(i+1) + "_2"


    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()   
    
############ lmb_epsilon_2
    
    
            
    fxx = 1.0
    fxy = 0.0
    fxz = strain1
    fyx = 0.0
    fyy = 1.0
    fyz = 0.0
    fzx = strain1
    fzy = 0.0
    fzz = 1.0
    
    
    exx=fxx-1.0
    exy=0.5*(fxy+fyx)
    exz=0.5*(fxz+fzx)
    eyy=fyy-1.0
    eyz=0.5*(fyz+fzy)
    ezz=fzz-1.0
    
    ax=1.0/math.sqrt(2)
    ay=0.0
    az=1.0/math.sqrt(2)
    
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c33,c44)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b21,b22,b3,b4,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_4_" + str(i+1) + "_1"
    
    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    
    
    ax=-1.0/math.sqrt(2)
    ay=0.0
    az=1.0/math.sqrt(2)
     
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c33,c44)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b21,b22,b3,b4,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_4_" + str(i+1) + "_2"


    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()  




lmb_name = "output_exact.dat"
    
dat = open(lmb_name,'w')
dat.write("lmba12 = ")
dat.write(repr(lmba12))
dat.write('\n')
dat.write("lmba22= ")
dat.write(repr(lmba22))
dat.write('\n')
dat.write("lmbg2= ")
dat.write(repr(lmbg2))
dat.write('\n')
dat.write("lmbe2= ")
dat.write(repr(lmbe2))
dat.write('\n')

dat.write("b12 = ")
dat.write(repr(b21*1e-6))
dat.write(" MPa")
dat.write('\n')
dat.write("b22= ")
dat.write(repr(b22*1e-6))
dat.write(" MPa")
dat.write('\n')
dat.write("b3= ")
dat.write(repr(b3*1e-6))
dat.write(" MPa")
dat.write('\n')
dat.write("b4= ")
dat.write(repr(b4*1e-6))
dat.write(" MPa")
dat.write('\n')



dat.close()
    
    

