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

    