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

from functions import fileio
from functions import polymer_funcs
from functions import elements_functions
import numpy as np
import pandas as pd
import matplotlib.pyplot as plot


expected_adducts=["Na","H","NH4","K","Ag","Ca","Li"]
expected_endgroups=["H2O","","C2H2O","CH3","CH3O","OH","O","C","C2","H2","C2H2","H"]

def pipeline1(infile,  expected_endgroups,expected_adducts):
    mz, intens,z = fileio.open_excel("C:/Users/FTICR_Kool_Kidz_PC/Downloads/Peg_41k_CID.xlsx",get_charge=True)
    masses=elements_functions.decharge(mz, z)
    maxzarr=np.sort(z)
    mmd=polymer_funcs.get_mmd_array(masses,elements_functions.formula_to_mass("C2H4O") )
    mmdx,mmdy= polymer_funcs.get_mmd_projections(masses,elements_functions.formula_to_mass("C2H4O") )
    peaks=polymer_funcs.get_peaks_in_projection(mmdx, mmdy)
    assignments=elements_functions.filter_common_adducts_and_endgroups(peaks, expected_adducts,expected_endgroups,1000,4,elements_functions.formula_to_mass("C2H4O"))
    assigned_peaks=[]
    assigned_residual_formulas=[]
    assigned_charges=[]
    refmasses=[]
    reflabels=[]
    for i in assignments:
        assigned_residual_formulas.append(i[0]+i[2])
        assigned_charges.append(i[1])
    for i in range(len(assigned_residual_formulas)):
        lists,labels,charges=polymer_funcs.make_lists_internal("C2H4O", 0, 100, assigned_residual_formulas[i], "",assigned_charges[i])
        refmasses.extend(lists)
        reflabels.extend(labels)
    #refmasses,reflabels,charges=polymer_funcs.make_lists("C2H4O",0,100, expected_endgroups,expected_adducts,5)
    mchecker=[]
    for i in range(len(mz)):
        mass=mz[i]
        intensity=intens[i]
        checker=np.abs(np.subtract(refmasses,mass))
        closest=np.argmin(checker)
        if abs(elements_functions.masserror(mass, refmasses[closest]))<2:
            assigned_peaks.append([mass,intensity,refmasses[closest],reflabels[closest],elements_functions.masserror(mass, refmasses[closest])])
            mchecker.append(mass)
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
    

mmd,rmmd, spectra, mmdpeaks,assignments=pipeline1("C:/Users/FTICR_Kool_Kidz_PC/Downloads/Peg_41k_CID.xlsx",expected_endgroups,expected_adducts)


    