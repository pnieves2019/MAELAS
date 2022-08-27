#!/usr/bin/env python3

import numpy as np
import fileinput
import sys, string, os
import ctypes
import matplotlib.pyplot as plt

from pymatgen.analysis.eos import EOS,EOSBase




f = open('ene_vs_vol_ferro.dat','r')
l = f.readlines()
f.close

x = []
y = []

for i in l:
    x.append(float(i.split()[0]))
    y.append(float(i.split()[1]))

vol = np.array(x)
ene = np.array(y)

#eos = EOS(eos_name='vinet')
eos = EOS(eos_name='murnaghan')
eos_fit = eos.fit(vol, ene)




print("EOS Ferromagnetic state:")

#eos_fit.plot()
print("B0=",eos_fit.b0_GPa)
print("E0=",eos_fit.e0)
print("v_FM=",eos_fit.v0)
print("dB0=",eos_fit.b1)


f = open('ene_vs_vol_para.dat','r')
l = f.readlines()
f.close

x = []
y = []

for i in l:
    x.append(float(i.split()[0]))
    y.append(float(i.split()[1]))

volp = np.array(x)
enep = np.array(y)

eosp = EOS(eos_name='murnaghan')
eos_fitp = eosp.fit(volp, enep)


print("------------------")

print("EOS Paramagnetic state:")


eos_fit.plot()
print("B0p=",eos_fitp.b0_GPa)
print("E0p=",eos_fitp.e0)
print("v_PM=",eos_fitp.v0)
print("dB0p=",eos_fitp.b1)

print("----------------")

print("Spontaneous volume magnetostriction:")

print("ws = (v_FM-v_PM)/v_PM = ",((eos_fit.v0-eos_fitp.v0)/eos_fitp.v0))

#print("lmb_alpha = ws/3 = ",((eos_fit.v0-eos_fitp.v0)/eos_fitp.v0)/3.0)

