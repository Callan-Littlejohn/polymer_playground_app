# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 11:01:47 2025

@author: FTICR_Kool_Kidz_PC
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plot

def open_excel(filename,get_charge=False):
    data=pd.read_excel(filename)
    mz=data["m/z"].values
    intens=data["I"].values
    
    for i in range(len(mz)):
        mz[i]=float(mz[i])
        intens[i]=float(intens[i])
    #intens=np.log10(intens)
    #intens=np.divide(intens, max(intens))
    if get_charge==True:
        z=data["z"].values
        for i in range(len(mz)):
            z[i]=int(z[i][:-1])
        return mz, intens, z
    return mz, intens

def open_csv(filename):
    data=pd.read_csv(filename)
    mz=data["m/z"].values
    z=data["z"].values
    intens=data["I"].values
    for i in range(len(mz)):
        mz[i]=float(mz[i])
        intens[i]=float(intens[i])
    return mz, intens

    