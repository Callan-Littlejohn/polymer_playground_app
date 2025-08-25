# -*- coding: utf-8 -*-
"""
----License----

Copyright Ó 2025  Callan Littlejohn,  Peter O'Connor, and The University of Warwick

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or  any later version. You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

For more details, see the GNU General Public License.

 

----Notes----

This software was produced by Callan Littlejohn as a Research fellow at the University of Warwick. It is intended to make the processing and analysis of polymer MS data data simpler and faster.


@author: Callan Littlejohn
"""

import numpy as np
import matplotlib.pyplot as plot
import pandas as pd
import functions.elements_functions as functions
import scipy.signal
import itertools
da=1.66054e-27
echar=1.60217662e-19
emass=0.0005485833

def get_mmd_array(masses,refmass):
    return np.mod(masses, refmass)

def get_mmd_projections(masses,refmass):
    mmd=get_mmd_array(masses,refmass)
    ppm=1E-6
    widths=np.divide(masses,float(1/ppm))
    for i in range(len(widths)):
        if widths[i]<0:
            widths[i]=abs(widths[i])*2
    xaxis=np.linspace(0,refmass,10000)
    yaxis=np.zeros(10000)
    zipped=zip(masses,mmd,widths)
    for (x,y,width) in zipped:
        guassian=np.exp(-0.5*((xaxis-y)/width)**2)/(np.sqrt(2*np.pi)*width)
        guassian=np.divide(guassian,max(guassian))
        yaxis+=guassian
    yaxis[0]=yaxis[0]+yaxis[-1]
    yaxis[-1]=0
    return xaxis,yaxis

def simulate_polydispersity(monomer,Mw,Mn,endcap):
    midmass=round(Mw/monomer)*monomer+endcap
    pdi=Mw/Mn
    mpeak=((2*Mw)-Mn)/(pdi)
    sigpeak=(Mw-Mn)/np.sqrt(2*pdi*pdi)
    peaks=[]
    offset=0
    peaks1,ints1=[],[]
    n=[]
    for i in range(10000):
        peaks1.append(offset+(monomer*i))
    ints=[]
    for i in peaks1:
        ints1.append(100*np.exp(-((i-mpeak)**2)/(2*sigpeak**2)))
    for i in range(len(ints1)):
        if ints1[i]>1:
            peaks.append(peaks1[i])
            ints.append(ints1[i])
            n.append(i)
    return peaks,ints,n

def make_lists_internal(monomer, n0,n1, endgroup,adduct,charge):
    mm=functions.formula_to_mass(monomer)
    am=functions.formula_to_mass(adduct)
    em=functions.formula_to_mass(endgroup)
    masses,labels,charges=[],[],[]
    for i in range(n0,n1):
        m=(i*mm+em+am-(emass*charge))/charge
        masses.append(round(m,6))
        labels.append("("+monomer+")"+str(i)+" "+endgroup+" "+adduct+"+")
        charges.append(charge)
    return masses,labels,charges

def make_lists(monomer, n0,n1, endgroup,adducts,max_charge):
    adducts.append("")
    potential_adduct_combos=[]
    potential_adduct_combo_iter=itertools.combinations_with_replacement(adducts, max_charge)
    charges=[]
    assigned=[]
    potential_adduct_combo_iter=sorted(set(potential_adduct_combo_iter))
    for i in potential_adduct_combo_iter:
        #charge=max_charge-((i.count("")))
        addcount=""
        for j in adducts:
            if j!=""and i.count(j)!=0:
                addcount=addcount+j+str(i.count(j)) 
            elif j=="":
                charge=max_charge-((i.count("")))
            #addcount.append(i.count(j))
        if addcount!="":
            charges.append(charge)
            potential_adduct_combos.append(addcount)
    potential_adduct_combos=potential_adduct_combos[:-1]
    charges=charges[:-1]
    #print(len(potential_adduct_combos),len(charges))
    masslist,labellist,chargesout=[],[],[]
    for i in range(len(potential_adduct_combos)):
        for j in endgroup:
            masses,labels,chargesout=make_lists_internal(monomer, n0, n1, j, potential_adduct_combos[i], charges[i])
            masslist.extend(masses)
            labellist.extend(labels)
    return masslist, labellist,chargesout


# worker funcs

def gaussian(x, amp, cent,width):
    return amp *np.exp(-((x-cent)**2)/(2*width**2))
def fitgaussian(x,y,po=None):
    params,covar=scipy.optimize.curve_fit(gaussian,x,y,po)
    return params

def get_peaks_in_projection(xaxis,yaxis,height=3):
    peaks,_=scipy.signal.find_peaks(yaxis,height=height)
    maxx=max(xaxis)
    width=10000/(maxx*4)
    widtheitherway=int(width/2)
    widtheitherway=30
    moi,intens=[],[]
    #print(peaks)
    for i in peaks:
        if i>width and i<10000-widtheitherway:
            try:
                peak=fitgaussian(xaxis[i-widtheitherway:i+widtheitherway], yaxis[i-widtheitherway:i+widtheitherway],po=[yaxis[i],xaxis[i],0.01])          
                moi.append(peak[1])
                intens.append(peak[0])
            except:
                ma=1
                moi.append(xaxis[i])
                intens.append(yaxis[i])
                #print("didnt work for :",xaxis[i])
                #plot.plot(xaxis[i-widtheitherway:i+widtheitherway], yaxis[i-widtheitherway:i+widtheitherway])
                #plot.show()
                
    return moi,intens

def gen_mmdpro_optimised(mz, refmass):
    xax = np.linspace(0, refmass, 10001)
    yaxis = np.zeros(10001)
    mmds = np.mod(mz, refmass)
    ppmwidths = mz / 1e6
    
    sqrt_2pi = np.sqrt(2 * np.pi)
    
    for i in range(len(mz)):
        diff = xax - mmds[i]
        exponent = -0.5 * (diff / ppmwidths[i])**2
        gaussian = np.exp(exponent) / (sqrt_2pi * ppmwidths[i])
        gaussian /= np.max(gaussian)
        yaxis += gaussian
    
    return xax, yaxis

def getidscore(mz,i,refmass):
    xax,yax=gen_mmdpro_optimised(mz,  refmass)
    #plot.plot(xax,yax)
    #plot.xlim(61.8,62)
    # yax=np.sort(yax)
    # minx=np.argmin(np.abs(np.subtract(yax,3.0)))
    #print(minx).
    
    peaks=scipy.signal.find_peaks(yax,threshold=3)[0]
    if len(peaks)==0:
        return(np.mean(yax))
    xx=[]
    for i in peaks:
            xx.append(yax[i])
        
    idscore=np.sum(xx[:])
    return idscore
