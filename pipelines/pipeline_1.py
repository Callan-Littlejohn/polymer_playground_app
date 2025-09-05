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
import itertools

import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plot
import time
import math


def pipeline2(infile,  expected_endgroups,expected_adducts,repeat_unit,maxcharge,maxtol,neg_mode=False): # targetted pipeline
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
    mip,mipd,iip=pwp([mz,intens],masses,mmd,peaks,tol*1.5)
    results=[]
    mmdtheos,theoforms,theocharges=mmdrefmaker(repeat_unit, expected_endgroups , expected_adducts, maxcharge)
    for i in range(len(peaks)):
        tol=(maxtol/1E6)*np.mean(mip[i])
        finderarray=np.abs(np.subtract(mmdtheos,peaks[i]))
        if min(finderarray)<tol:
            closestind=np.argmin(finderarray)
            results.append([peaks[i],mmdtheos[closestind],theoforms[closestind],mip[i],mipd[i],theocharges[closestind],iip[i]])
    assignments=[]
    for i in results:
        addform=i[2]
        maxn=math.ceil(max(i[4])/elements_functions.formula_to_mass(repeat_unit))
        minn=math.floor(min(i[4])/elements_functions.formula_to_mass(repeat_unit))
        masslist, labellist,chargesout=polymer_funcs.make_lists_internal(repeat_unit, minn, maxn, addform,"",i[5],neg_mode)
        #print(masslist,i[5])
        for j in range(len(i[3])):
            closestind=np.argmin(np.abs(np.subtract(masslist,i[3][j])))
            me=elements_functions.masserror(i[3][j], masslist[closestind])           
            if abs(me)<2*maxtol:
                assignments.append([i[3][j],i[6][j],masslist[closestind],labellist[closestind],me,i[5]])
                #([labellist[closestind],masslist[closestind],i[3][j],elements_functions.masserror(i[3][j], masslist[closestind])])
    assigned_peaks=assignments
    massestoplot,istoplot=[],[]
    for i in assigned_peaks:
        massestoplot.append(i[0])
        istoplot.append(i[1])
    # plot.stem(mz,intens,markerfmt="",linefmt="b")
    # plot.stem(massestoplot,istoplot,linefmt="r",markerfmt="")
    residual_spectram,residual_spectrai=[],[]
    for i in range(len(mz)):
        if mz[i] not in massestoplot:
            residual_spectram.append(mz[i])
            residual_spectrai.append(intens[i])
    mmdrx,mmdry=polymer_funcs.get_mmd_projections(residual_spectram,elements_functions.formula_to_mass("C2H4O") )
    respeaks=polymer_funcs.get_peaks_in_projection(mmdrx, mmdry,height=10)
    #assignments=np.transpose(assigned_peaks)
    return [mmdx,mmdy],[mmdrx,mmdry],[mz,intens,mmd],[peaks,respeaks],assigned_peaks
    return assignments

def findcomp(mass,elements,emin, emax,tol,minmass,refmass):
    while mass<minmass:
        results=acmasstest(mass, elements, emin, emax, tol)
        if results!=0:
            return results
        mass=mass+refmass
    return 0

def pwp(masses,m, mmd, peaks, tol):
    ints=np.array(masses[1])
    masses=np.array(masses[0])
    massinpeaks=[]
    intsinpeaks=[]
    m2=[]
    for i in peaks:
        #print(i)
        inds=np.where(np.abs(np.subtract(mmd,i))<tol)
        #print(inds)
        #print(masses[inds])
        massinpeaks.append(masses[inds])
        intsinpeaks.append(ints[inds])
        m2.append(m[inds])
        
    return massinpeaks,m2,intsinpeaks
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

def mmdrefmaker(repeat_unit,endgroups, adducts, maxcharge):
    adducts.append("")
    potential_adduct_combos=[]
    potential_adduct_combo_iter=itertools.combinations_with_replacement(adducts, maxcharge)
    charges=[]
    assigned=[]
    for i in potential_adduct_combo_iter:
        charge=maxcharge-((i.count("")))
        charges.append(charge)
        potential_adduct_combos.append("".join(i))
    potential_adduct_combos=potential_adduct_combos[:-1]
    charges=charges[:-1]
    combos=itertools.product(endgroups,potential_adduct_combos)
    combocharges=itertools.product(endgroups,charges)
    charges=[]
    for i in combocharges:
        charges.append(int(i[1]))
    mmdtheo=[]
    combos2=[]
    for i in combos:
        form=" ".join(i)
        combos2.append(form)
        mmdtheo.append(elements_functions.formula_to_mass(form)%elements_functions.formula_to_mass(repeat_unit))
    return mmdtheo,combos2,charges
        