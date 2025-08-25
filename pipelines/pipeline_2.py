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

from pyximport import install; install()
from functions import fileio
from functions import polymer_funcs
from functions import elements_functions
from functions import acuratemass
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plot
import time


expected_adducts=["Na","H","NH4","K","Ag","Ca","Li"]
expected_endgroups=["H2O","","C2H2O","CH3","CH3O","OH","O","C","C2","H2","C2H2","H"]

def pipeline2(infile,  repeat_unit,elements,mins,maxs,maxtol): # targetted pipeline
    mz, intens,z = fileio.open_excel(infile,get_charge=True)
    #mz=np.multiply(mz,z)
    #z=np.ones(len(mz))
    masses=np.array(elements_functions.decharge(mz, z))
    #masses=mz
    maxzarr=np.sort(z)
    mmd=np.array(polymer_funcs.get_mmd_array(masses,elements_functions.formula_to_mass(repeat_unit) ))
    mmdx,mmdy= polymer_funcs.get_mmd_projections(masses,elements_functions.formula_to_mass(repeat_unit) )
    peaks,peakintens=polymer_funcs.get_peaks_in_projection(mmdx, mmdy)
    
    cyt=[]
    pyt=[]
    tol=(maxtol/1E6)*np.mean(mz)
    mip=pwp(mz,mmd,peaks,tol)
    results=[]
    for i in range(len(peaks)):
        #try:
        if 1==1:   
            tol=(maxtol/1E6)*np.mean(mip[i])
            #t0=time.time()
            #form=acuratemass.accuratemass2(i,np.mean(mz))
            #t1=time.time()
            
            form2=findcomp(peaks[i], elements, mins, maxs, tol,min(mip[i]),elements_functions.formula_to_mass(repeat_unit))
            results.append(form2)
            #print(i,form2)
            # t2=time.time()
            # cyt.append(t1-t0)
            # pyt.append(t2-t1)
            #print(form,form2)
        # except:
        #     print("not for",i)
    #print(np.mean(cyt),np.mean(pyt))
    #plot.stem(peaks,peakintens)
    #for i in range(len(peaks)):
    return results

def findcomp(mass,elements,emin, emax,tol,minmass,refmass):
    while mass<minmass:
        results=acmasstest(mass, elements, emin, emax, tol)
        if results!=0:
            return results
        mass=mass+refmass
    return 0

def pwp(masses, mmd, peaks, tol):
    masses=np.array(masses)
    massinpeaks=[]
    for i in peaks:
        #print(i)
        inds=np.where(np.abs(np.subtract(mmd,i))<tol)
        #print(inds)
        #print(masses[inds])
        massinpeaks.append(masses[inds])
    return massinpeaks
def acmasstest(mass,elements,emin, emax,tol):
    melcombos=[]
    widearray=[]
    elnames=[]
    elcombos=[]
    results=[]
    for i in range(len(elements)):
        elnums=np.linspace(emin[i],emax[i],emax[i]-emin[i]+1)
        mel=elements_functions.formula_to_mass(elements[i])
        melcombos.append(np.multiply(elnums,mel))
        temp=[]
        for j in elnums:
            temp.append(elements[i]+str(int(j)))
        elnames.append(temp)
            
    for combo in itertools.product(*elnames): 
        elcombos.append("".join(combo))
    for combo in itertools.product(*melcombos):
        combomass=sum(combo)
        widearray.append(combomass)
    widearray=np.array(widearray)
    widearray[widearray>mass+tol]=0
    widearray[widearray<mass-tol]=0
    if sum(widearray)==0:
        return 0
    #return widearray[widearray!=0], elcombos[widearray!=0]
    for i in range(len(widearray)):
        if widearray[i]!=0:
            results.append([elcombos[i],widearray[i]])
    return results
    

#r=pipeline2("C:/Users/FTICR_Kool_Kidz_PC/Downloads/Peg_41k_CID.xlsx",expected_endgroups,expected_adducts,"C2H4O",4,3)
# #print(acmasstest(58.02,["H","C","O"],[0,2,0],[10,3,1],0,0,0.5) )
#print(r)