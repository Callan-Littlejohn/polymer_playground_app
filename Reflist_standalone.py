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

import figure_pipelines as fp
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
import math


class gui():
    def __init__(self):
        
        self.make_mainscreen()
        
    def make_mainscreen(self):
        self.mainscreen=tk.Tk()
        self.mainscreen.title("Polymer ref standalone")
        labelmono=ttk.Label(self.mainscreen,text="Repeat Unit").grid(row=0,column=0,sticky="NSEW")
        labeladd=ttk.Label(self.mainscreen,text="Adduct").grid(row=1,column=0,sticky="NSEW")
        labelend=ttk.Label(self.mainscreen,text="End Group").grid(row=2,column=0,sticky="NSEW")
        labelminmax=ttk.Label(self.mainscreen,text="Min/Max m/z").grid(row=3,column=0,sticky="NSEW")
        labelz=ttk.Label(self.mainscreen,text="Charge").grid(row=4,column=0,sticky="NSEW")
        
        self.Entrymono=ttk.Entry(self.mainscreen)
        self.Entrymono.grid(row=0,column=1,columnspan=2,sticky="NSEW")
        self.Entryadd=ttk.Entry(self.mainscreen)
        self.Entryadd.grid(row=1,column=1,columnspan=2,sticky="NSEW")
        self.Entryend=ttk.Entry(self.mainscreen)
        self.Entryend.grid(row=2,column=1,columnspan=2,sticky="NSEW")
        self.Entrymin=ttk.Entry(self.mainscreen)
        self.Entrymin.grid(row=3,column=1,sticky="NSEW")
        self.Entrymax=ttk.Entry(self.mainscreen)
        self.Entrymax.grid(row=3,column=2,sticky="NSEW")
        self.Entryz=ttk.Entry(self.mainscreen)
        self.Entryz.grid(row=4,column=1,columnspan=2,sticky="NSEW")
        
        self.mrbutton=ttk.Button(self.mainscreen,text="Make .ref", command=self.makeref).grid(row=5,column=0,columnspan=3,sticky="NSWE")
        
        sv_ttk.set_theme("dark")
        self.mainscreen.mainloop()
        
    def makeref(self):
        monomer,adduct,endgroup, maxcharge,minmz,maxmz=self.Entrymono.get(),self.Entryadd.get(),self.Entryend.get(),self.Entryz.get(),self.Entrymin.get(),self.Entrymax.get()
        mmass=elements_functions.formula_to_mass(monomer)
        maxn=math.ceil(float(maxmz)/(mmass/int(maxcharge)))
        minn=math.ceil(float(minmz)/(mmass/int(maxcharge)))
        masses,labels,charges=polymer_funcs.make_lists_internal(monomer, minn, maxn, endgroup,adduct, int(maxcharge))
        #print(labels,masses)
        outfile=tk.filedialog.asksaveasfilename()
        with open(outfile+".ref","w") as f:
            for i in range(len(masses)):
                labels[i]=labels[i].replace(" ","")
                f.write(f"{labels[i]}\t{masses[i]}\t{charges[i]}\n")
            f.close()

r=gui()