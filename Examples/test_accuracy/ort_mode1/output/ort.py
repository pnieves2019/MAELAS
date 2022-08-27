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



b1=43.1*10**6
b2=-34.2*10**6
b3=60.7*10**6
b4=-34.3*10**6
b5=-42.4*10**6
b6=55.4*10**6
b7=35.4*10**6
b8=-22.6*10**6
b9=38.7*10**6


a=4.0686275541928172 
b=10.3156917156588115
c=3.8956767252497921
v0=(a*b*c)*((10.0**(-10.0))**3.0)


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


ndist = 7
strain = 0.01


def elas0(eexx,eexy,eexz,eeyy,eeyz,eezz,cc11,cc12,cc13,cc22,cc23,cc33,cc44,cc55,cc66):
    
    ene=0.5*cc11*eexx**2.0+0.5*cc22*eeyy**2.0+cc12*eexx*eeyy
    ene=ene+cc13*eexx*eezz+cc23*eeyy*eezz+0.5*cc33*eezz**2.0
    ene=ene+2.0*cc44*eeyz**2.0+2.0*cc55*eexz**2.0+2.0*cc66*eexy**2.0
    
    return ene

def em0(eexx,eexy,eexz,eeyy,eeyz,eezz,bb1,bb2,bb3,bb4,bb5,bb6,bb7,bb8,bb9,aax,aay,aaz):
    
    ene0 = bb1*aax**2.0*eexx + bb2*aay**2.0*eexx + bb3*aax**2.0*eeyy + bb4*aay**2.0*eeyy
    ene0 = ene0 + bb5*aax**2.0*eezz + bb6*aay**2.0*eezz
    ene0 = ene0 + 2.0*bb7*aax*aay*eexy + 2.0*bb8*aax*aaz*eexz + 2.0*bb9*aay*aaz*eeyz
    
    return ene0






for i in range(int(ndist)):


    strain1 = - float(strain)+2.0*(float(strain)/(float(ndist)-1.0))*i

    print("strain", strain1)


        #Generation POSCAR file

        #lambda_1

    fxx = 1.0 + strain1
    fyy = 1.0/math.sqrt(fxx)
    fzz = fyy
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
    em = em0(exx,exy,exz,eyy,eyz,ezz,b1,b2,b3,b4,b5,b6,b7,b8,b9,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_1_" + str(i+1) + "_1"
    
    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    
    
    
    
    ax=0.0
    ay=0.0
    az=1.0
     
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c22,c23,c33,c44,c55,c66)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b1,b2,b3,b4,b5,b6,b7,b8,b9,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_1_" + str(i+1) + "_2"


    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    

    
############ lmb_2
    
    
    fxx = 1.0 + strain1
    fyy = 1.0/math.sqrt(fxx)
    fzz = fyy
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
    ay=1.0
    az=0.0
    
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c22,c23,c33,c44,c55,c66)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b1,b2,b3,b4,b5,b6,b7,b8,b9,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_2_" + str(i+1) + "_1"
    
    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    
    
    
    
    ax=0.0
    ay=0.0
    az=1.0
     
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c22,c23,c33,c44,c55,c66)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b1,b2,b3,b4,b5,b6,b7,b8,b9,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_2_" + str(i+1) + "_2"


    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    
    
    
############ lmb_3
    
    
    fyy = 1.0 + strain1
    fxx = 1.0/math.sqrt(fyy)
    fzz = fxx
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
    em = em0(exx,exy,exz,eyy,eyz,ezz,b1,b2,b3,b4,b5,b6,b7,b8,b9,ax,ay,az)
    
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
    ay=0.0
    az=1.0
     
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c22,c23,c33,c44,c55,c66)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b1,b2,b3,b4,b5,b6,b7,b8,b9,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_3_" + str(i+1) + "_2"


    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()   
    
############ lmb_4
    
    fyy = 1.0 + strain1
    fxx = 1.0/math.sqrt(fyy)
    fzz = fxx
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
    ay=1.0
    az=0.0
    
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c22,c23,c33,c44,c55,c66)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b1,b2,b3,b4,b5,b6,b7,b8,b9,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_4_" + str(i+1) + "_1"
    
    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    
    
    ax=0.0
    ay=0.0
    az=1.0
     
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c22,c23,c33,c44,c55,c66)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b1,b2,b3,b4,b5,b6,b7,b8,b9,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_4_" + str(i+1) + "_2"


    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()  


############ lmb_5
    
    fzz = 1.0 + strain1
    fyy = 1.0/math.sqrt(fzz)
    fxx = fyy
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
    em = em0(exx,exy,exz,eyy,eyz,ezz,b1,b2,b3,b4,b5,b6,b7,b8,b9,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_5_" + str(i+1) + "_1"
    
    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    
    
    ax=0.0
    ay=0.0
    az=1.0
     
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c22,c23,c33,c44,c55,c66)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b1,b2,b3,b4,b5,b6,b7,b8,b9,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_5_" + str(i+1) + "_2"


    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close() 



############ lmb_6
    
    fzz = 1.0 + strain1
    fyy = 1.0/math.sqrt(fzz)
    fxx = fyy
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
    ay=1.0
    az=0.0
    
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c22,c23,c33,c44,c55,c66)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b1,b2,b3,b4,b5,b6,b7,b8,b9,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_6_" + str(i+1) + "_1"
    
    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    
    
    ax=0.0
    ay=0.0
    az=1.0
     
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c22,c23,c33,c44,c55,c66)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b1,b2,b3,b4,b5,b6,b7,b8,b9,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_6_" + str(i+1) + "_2"


    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close() 



############ lmb_7
    
    const = (1.0/(1.0-(strain1*0.5)**2.0))**(1.0/3.0)
            
    fxx = const
    fxy = const*strain1*0.5*(b/a)
    fxz = 0.0
    fyx = const*strain1*0.5*(a/b)
    fyy = const
    fyz = 0.0
    fzx = 0.0
    fzy = 0.0
    fzz = const
    
    
    exx=fxx-1.0
    exy=0.5*(fxy+fyx)
    exz=0.5*(fxz+fzx)
    eyy=fyy-1.0
    eyz=0.5*(fyz+fzy)
    ezz=fzz-1.0
    
    ax=1.0/math.sqrt(2)
    ay=1.0/math.sqrt(2)
    az=0.0
    
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c22,c23,c33,c44,c55,c66)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b1,b2,b3,b4,b5,b6,b7,b8,b9,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_7_" + str(i+1) + "_1"
    
    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    
    
    ax=0.0
    ay=0.0
    az=1.0
     
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c22,c23,c33,c44,c55,c66)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b1,b2,b3,b4,b5,b6,b7,b8,b9,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_7_" + str(i+1) + "_2"


    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close() 

############ lmb_8
    
    const = (1.0/(1.0-(strain1*0.5)**2.0))**(1.0/3.0)
            
    fxx = const
    fxy = 0.0
    fxz = const*strain1*0.5*(c/a)
    fyx = 0.0
    fyy = const
    fyz = 0.0
    fzx = const*strain1*0.5*(a/c)
    fzy = 0.0
    fzz = const
    
    exx=fxx-1.0
    exy=0.5*(fxy+fyx)
    exz=0.5*(fxz+fzx)
    eyy=fyy-1.0
    eyz=0.5*(fyz+fzy)
    ezz=fzz-1.0
    
    ax=1.0/math.sqrt(2)
    ay=0.0
    az=1.0/math.sqrt(2)
    
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c22,c23,c33,c44,c55,c66)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b1,b2,b3,b4,b5,b6,b7,b8,b9,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_8_" + str(i+1) + "_1"
    
    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    
    
    ax=0.0
    ay=0.0
    az=1.0
     
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c22,c23,c33,c44,c55,c66)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b1,b2,b3,b4,b5,b6,b7,b8,b9,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_8_" + str(i+1) + "_2"


    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    
    
############ lmb_9
    
    const = (1.0/(1.0-(strain1*0.5)**2.0))**(1.0/3.0)
            
    fxx = const
    fxy = 0.0
    fxz = 0.0
    fyx = 0.0
    fyy = const
    fyz = const*strain1*0.5*(c/b)
    fzx = 0.0
    fzy = const*strain1*0.5*(b/c)
    fzz = const
    
    
    exx=fxx-1.0
    exy=0.5*(fxy+fyx)
    exz=0.5*(fxz+fzx)
    eyy=fyy-1.0
    eyz=0.5*(fyz+fzy)
    ezz=fzz-1.0
    
    ax=0.0
    ay=1.0/math.sqrt(2)
    az=1.0/math.sqrt(2)
    
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c22,c23,c33,c44,c55,c66)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b1,b2,b3,b4,b5,b6,b7,b8,b9,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_9_" + str(i+1) + "_1"
    
    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close()
    
    
    ax=0.0
    ay=0.0
    az=1.0
     
    elas = elas0(exx,exy,exz,eyy,eyz,ezz,c11,c12,c13,c22,c23,c33,c44,c55,c66)
    em = em0(exx,exy,exz,eyy,eyz,ezz,b1,b2,b3,b4,b5,b6,b7,b8,b9,ax,ay,az)
    
    etot=v0*(elas+em)*6241509000000000000.0
    
    osz_name = "OSZICAR_9_" + str(i+1) + "_2"


    dat = open(osz_name,'w')
    dat.write("           ")
    dat.write("{:10.14f}".format(etot))
    dat.write("                                       ")
    dat.write('\n')
    dat.write('end file')
    dat.close() 


lmb_name = "output_exact.dat"
    
dat = open(lmb_name,'w')
dat.write("lmb1 = ")
dat.write(repr(lambda_ortho1))
dat.write('\n')
dat.write("lmb2= ")
dat.write(repr(lambda_ortho2))
dat.write('\n')
dat.write("lmb3= ")
dat.write(repr(lambda_ortho3))
dat.write('\n')
dat.write("lmb4= ")
dat.write(repr(lambda_ortho4))
dat.write('\n')
dat.write("lmb5= ")
dat.write(repr(lambda_ortho5))
dat.write('\n')
dat.write("lmb6= ")
dat.write(repr(lambda_ortho6))
dat.write('\n')
dat.write("lmb7= ")
dat.write(repr(lambda_ortho7))
dat.write('\n')
dat.write("lmb8= ")
dat.write(repr(lambda_ortho8))
dat.write('\n')
dat.write("lmb9= ")
dat.write(repr(lambda_ortho9))
dat.write('\n')


dat.write("b1 = ")
dat.write(repr(b1*1e-9))
dat.write(" GPa")
dat.write('\n')
dat.write("b2= ")
dat.write(repr(b2*1e-9))
dat.write(" GPa")
dat.write('\n')
dat.write("b3= ")
dat.write(repr(b3*1e-9))
dat.write(" GPa")
dat.write('\n')
dat.write("b4= ")
dat.write(repr(b4*1e-9))
dat.write(" GPa")
dat.write('\n')
dat.write("b5= ")
dat.write(repr(b5*1e-9))
dat.write(" GPa")
dat.write('\n')
dat.write("b6= ")
dat.write(repr(b6*1e-9))
dat.write(" GPa")
dat.write('\n')
dat.write("b7= ")
dat.write(repr(b7*1e-9))
dat.write(" GPa")
dat.write('\n')
dat.write("b8= ")
dat.write(repr(b8*1e-9))
dat.write(" GPa")
dat.write('\n')
dat.write("b9= ")
dat.write(repr(b9*1e-9))
dat.write(" GPa")
dat.write('\n')


dat.close()
    
    

