# -*- coding: utf-8 -*-
"""
Created on Sat Aug 30 15:00:49 2025

@author: FTICR_Kool_Kidz_PC
"""

#import figure_pipelines as fp
import numpy as np
import pandas as pd
import tkinter as tk
import pandas as pd
import pipelines.pipeline_1
import pipelines.pipeline_2
import matplotlib.pyplot as plot
from matplotlib.gridspec import GridSpec
from functions import elements_functions
from functions import polymer_funcs
from tkinter import ttk
import sv_ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import scipy.signal
import matplotlib
import time
from matplotlib.figure import Figure
import os
from scipy.optimize import minimize


#universal constants
da=1.66054e-27
echar=1.60217662e-19
emass=0.0005485833

def me_calc(dm,m):
    me=(dm/m)*1E6
    return me
def getfreq(mass,charge,magnet): #get the frequency of a m/z
    mass=mass*da #work out charge in kg
    charge=charge*echar #work out charge in C
    freq=(magnet*charge)/(2*np.pi*mass) #determine frequency using base ICR equation . assumes no electric fields this is close then calibration eq get closer
    return freq # returns frequency

def recalibrate(A,B,freqs):
    Aterm=np.multiply(np.reciprocal(freqs),A)
    Bterm=np.multiply(np.reciprocal(np.power(freqs,2)),B)
    return np.add(Aterm,Bterm)
def recalibrate2(x,coeffs):
    return(np.subtract(np.divide(x,1/coeffs[0]),coeff[1]))

c=['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    
data=pd.read_excel("C:/Users/FTICR_Kool_Kidz_PC/Downloads/Peg_41k_CID.xlsx")#"C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_2006hrs.xlsx"
mz=data["m/z"].values
intens=data["I"].values
refmass=elements_functions.formula_to_mass("C2H4O")#"C8H8"
mmd=polymer_funcs.get_mmd_array(mz,refmass)
x,y=polymer_funcs.get_mmd_projections(mz,refmass)
peaks,peakints=polymer_funcs.get_peaks_in_projection(x,y)
tol=(2/1E6)*np.mean(mz)
pwp=[]
# print(np.mod(mz,refmass))
# print(mmd)
pwp,iwp=(pipelines.pipeline_2.pwp(mz,mmd,peaks,tol,ints=intens))
errors=[]
mzs=[]
fj=[]
fi=[]
y=[]
cs=[]
ints=[]
a=0
md=[]
#print(pwp)
for i in range(len(pwp)):
    #pwpi=np.sort(pwp[i])
    pwpi=pwp[i]
    iwpi=iwp[i]
    sub=np.subtract(pwpi,refmass)
    
    for j in range(len(sub)):
        if j!=0:
            error=sub[j]-pwpi[j-1]
            if abs(error)<1:
                md.append(error)
                mzs.append(pwpi[j])
                errors.append(me_calc(error, sub[j]))
                fj.append(getfreq(sub[j], 1, 15))
                fi.append(getfreq(pwpi[j-1], 1, 15))
                y.append(refmass)
                ints.append(iwpi[j])
                cs.append(c[a])
    a=a+1
    if a>len(c)-1:
        a=0

plot.scatter(mzs,md,s=np.divide(ints,max(ints)/5),c=cs)
coeff=np.polyfit(mzs,md,1)
print(coeff)
bfy=np.poly1d(coeff)(mzs)
plot.plot(mzs,bfy)
plot.show()
print(np.std(errors))
# md2=recalibrate2(md, coeff)
# plot.scatter(mzs,md2,s=np.divide(ints,max(ints)/5),c=cs)
# coeff2=np.polyfit(mzs,md2,1)
# print(coeff2)
# bfy2=np.poly1d(coeff2)(mzs)
# plot.plot(mzs,bfy2)
# plot.show()
mz2=[]
for i in mz:
    newmz=i+(i/(1/coeff[0]))-coeff[1]
    mz2.append(newmz)
#plot.stem(mz2,intens)

mz=mz2
refmass=elements_functions.formula_to_mass("C2H4O")#"C8H8"
mmd=polymer_funcs.get_mmd_array(mz,refmass)
x,y=polymer_funcs.get_mmd_projections(mz,refmass)
peaks,peakints=polymer_funcs.get_peaks_in_projection(x,y)
tol=(2/1E6)*np.mean(mz)
pwp=[]
pwp,iwp=(pipelines.pipeline_2.pwp(mz,mmd,peaks,tol,ints=intens))
errors=[]
mzs=[]
fj=[]
fi=[]
y=[]
cs=[]
ints=[]
a=0
md=[]
#print(pwp)
for i in range(len(pwp)):
    #pwpi=np.sort(pwp[i])
    pwpi=pwp[i]
    iwpi=iwp[i]
    sub=np.subtract(pwpi,refmass)
    
    for j in range(len(sub)):
        if j!=0:
            error=sub[j]-pwpi[j-1]
            if abs(error)<1:
                md.append(error)
                mzs.append(pwpi[j])
                errors.append(me_calc(error, sub[j]))
                fj.append(getfreq(sub[j], 1, 15))
                fi.append(getfreq(pwpi[j-1], 1, 15))
                y.append(refmass)
                ints.append(iwpi[j])
                cs.append(c[a])
    a=a+1
    if a>len(c)-1:
        a=0
#md=np.multiply(md,ints)
plot.scatter(mzs,md,s=np.divide(ints,max(ints)/5),c=cs)

coeff=np.polyfit(mzs,md,1)
print(coeff)
bfy=np.poly1d(coeff)(mzs)
plot.plot(mzs,bfy)
plot.show()
print(np.std(errors))























#X = np.vstack([X1, X2]).T
#minimise()
# X1 = np.subtract(np.reciprocal(fj),np.reciprocal(fi))
# X2 =np.subtract(np.reciprocal(np.power(fj,2)),np.reciprocal(np.power(fi,2)))
# A_B = np.linalg.lstsq(X, y, rcond=None)[0]
# A=A_B[0]
# B=A_B[1]
# print("A:", A)
# print("B:", B)
# mzs=recalibrate(A, B, fj)
# print(mzs)

#plot.scatter(mzs,errors)
# xs,ys=[],[]
# xax = np.linspace(-1, 1, 100001)
# sqrt_2pi = np.sqrt(2 * np.pi)
# yaxis = np.zeros(100001)
# for i in range(len(mzs)):
#     me=me_calc(errors[i], mzs[i])
#     width=abs(me)/1
#     diff = xax - me
#     exponent = -0.5 * (diff / width)**2
#     gaussian = np.abs(np.exp(exponent) / (sqrt_2pi * width))
#     #plot.plot(xax,np.abs(gaussian))
#     gaussian = np.divide(gaussian,np.max(gaussian))
#     yaxis += gaussian
    
#plot.plot(xax,yaxis)
#print(np.std((errors)))

