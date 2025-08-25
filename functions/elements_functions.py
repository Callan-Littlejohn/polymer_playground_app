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

import pandas as pd
import matplotlib.pyplot as plot
import re
import numpy as np
import itertools
import math

elementsfile="isotopes_acurate.xlsx"

c13= 13.0033548
c13diff=c13-12
da=1.66054e-27
echar=1.60217662e-19
emass=0.0005485833

common_elements=["C","H","N","O","S"]

def getelements():
    eldata=pd.read_excel(elementsfile)
    elements={}
    element_info={}
    elsym=eldata["Symbol"].values
    elmass=eldata["Mass"].values
    elabundace=eldata["Abund."].values
    for i in range(len(elmass)):
        if type(elmass[i])==str:
            elmass[i]=float(elmass[i].split("(")[0])
        if type(elabundace[i])==str:
            elabundace[i]=float(elabundace[i].split("(")[0])
    
    
    for i in range(len(elsym)):
        if type(elsym[i])==str:
            name=elsym[i].split("(")[0]
        if name not in elements:
            elements.update({name:elmass[i]})
            element_info.update({name:{"masses":[elmass[i]],"abundance":[elabundace[i]]}})
        else:
            element_info[name]["masses"].append(elmass[i])
            element_info[name]["abundance"].append(elabundace[i])
    return elements,element_info

def formula_to_array(formula):
    # Use regular expressions to split the formula into elements and counts
    elements = re.findall(r'([A-Z][a-z]*)(-?\d*)', formula)
    
    # Initialize an empty array to store the elements and counts
    result = []
    
    for element, count in elements:
        # If the count is empty, assume it's 1 atom
        if count == '':
            count = 1
        else:
            count = int(count)
        
        # Add the element and count to the result array
        result.append([element, count])
    
    return result

def formula_to_mass(formula_ui):
    formula=formula_to_array(formula_ui)
    mass=0
    for i in formula:
        mass=mass+(elements[i[0]]*i[1])
    return mass

def decharge(mz,z):
    m= np.multiply(mz,z)
    es=np.multiply(z,emass)
    m=np.add(m, es)
    # m=np.round(m,decimals=6)
    for i in range(len(m)):
        m[i]=round(m[i],6)
    return m

def masserror(obs, theo):
    return((theo-obs)/theo)*1E6

def masserror_at_range(obs, theo,at):
    return((theo-obs)/at)*1E6

def determine_potential_adducts(masses,adducts,center):
    potential_adducts=[]
    for i in masses:
        potential_adducts.append([])
    for i in masses:
        temp=[]
        for j in adducts:
            for k in adducts:
                if 1==1:
                    mass=i-formula_to_mass(j)
                    altmasses=np.subtract(masses,formula_to_mass(k))
                    for l in range(len(altmasses)):
                        if abs(masserror_at_range(altmasses[l],mass,center))<10 and altmasses[l]!=mass:
                            print(i,mass,masses[l],j,k)
                            potential_adducts[l].append(k)
    return potential_adducts

def filter_common_adducts_and_endgroups(masses, adducts, endgroups, center,maxcharge,refmass):
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
    #masses=np.subtract(masses,emass)
    for i in range(len(potential_adduct_combos)):
        for j in endgroups:
            massend=(formula_to_mass(potential_adduct_combos[i])+formula_to_mass(j))%refmass
            #print(potential_adduct_combos[i],j,massend)
            for k in masses:
                m=k+(emass*charges[i])
                if abs(masserror_at_range(m, massend, center))<10:
                    assigned.append([potential_adduct_combos[i],charges[i],j,k])
    return assigned

def exact_mass_calculator(masses, adducts, additional_el,max_charge,refmass,center,additional):
    common_el=common_elements
    masses_assigned,assigned_formula,mes=[],[],[]
    masses=np.add(masses, additional)
    for i in range(max_charge):
        charge=i+1
        for adduct in adducts:
            m_add_removed=np.subtract(np.multiply(masses,charge), ((formula_to_mass(adduct)-emass)*charge))
            for j in range(len(m_add_removed)):
                if m_add_removed[j]<0:
                    m_add_removed[j]=m_add_removed[j]+formula_to_mass(refmass)
            max_el={}
            for el in common_el:
                max_el.update({el:int((max(m_add_removed)/ formula_to_mass(el)))})
            for O in range(max_el["O"]):
                #mass=formula_to_mass("O"+str(O))
                for N in range(max_el["N"]):
                   #mass=formula_to_mass("N"+str(N)+"O"+str(O)) 
                   for C in range(max_el["C"]):
                       #mass=formula_to_mass("C"+str(C)+"N"+str(N)+"O"+str(O)) 
                       for H in range(max_el["H"]):
                           mass=formula_to_mass("C"+str(C)+"H"+str(H)+"N"+str(N)+"O"+str(O)) 
                           masschecker=np.abs(np.subtract(m_add_removed,mass))
                           if min(masschecker)<0.2:
                               closest=m_add_removed[np.argmin(masschecker)]
                               if abs(masserror_at_range(closest, mass, center))<2:
                                   masses_assigned.append(masses[np.argmin(masschecker)])
                                   assigned_formula.append("C"+str(C)+"H"+str(H)+"N"+str(N)+"O"+str(O))
                                   mes.append(masserror_at_range(closest, mass, center))
    return masses_assigned,assigned_formula,mes

def ret_mass_diffs():
    mds=[]
    for i in elements:
        m=elements[i]
        md=m-math.floor(m)
        mds.append([i,md])
    return mds



mds=[]
abunds=[]
key=[]
md2=[]
elements,elements_info=getelements()
# for i in elements_info:
#     masses=elements_info[i]["masses"]
#     ab=elements_info[i]["abundance"]
#     #print(i,masses,ab)
#     if len(masses)==1:
#         mds.append([0])
#         abunds.append([1])
#     else:
#         temp=[]
#         for j in range(len(masses)):
#             temp.append(masses[j]-masses[0])
#             if masses[j]-masses[0]<1.1 and masses[j]-masses[0]>0.9:
#                 md2.append(masses[j]-masses[0])
#                 key.append(i)
#         mds.append(temp)
#         abunds.append(ab)
# # print(abunds,mds)
# # for i in range(len(mds)):
# #     plot.scatter(mds[i],abunds[i])
# # plot.xlim(0.9,1.1)
# plot.bar(key, md2)
# plot.ylim(0.98,1.02)