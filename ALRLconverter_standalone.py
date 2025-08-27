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
        self.infile,self.outfile=0,0
        self.make_mainscreen()
        
    def make_mainscreen(self):
        self.mainscreen=tk.Tk()
        self.mainscreen.title("AL->RL converter")
        labelin=ttk.Label(self.mainscreen,text="infile").grid(row=0,column=0)
        labelout=ttk.Label(self.mainscreen,text="outfile").grid(row=1,column=0)
        self.infilepathlabel=ttk.Label(self.mainscreen,text="  choose peak list  ")
        self.infilepathlabel.grid(row=0,column=1,sticky="NSEW")
        self.infilebrowsebutton=ttk.Button(self.mainscreen,text="browse",command=self.infilebrowse).grid(row=0,column=2, sticky="NSEW")
        self.outfilepathlabel=ttk.Label(self.mainscreen,text="  choose peak list  ")
        self.outfilepathlabel.grid(row=1,column=1,sticky="NSEW")
        self.outfilebrowsebutton=ttk.Button(self.mainscreen,text="browse",command=self.outfilebrowse).grid(row=1,column=2, sticky="NSEW")
        butt=ttk.Button(self.mainscreen,text="convert",command=self.convert).grid(row=2,column=0,columnspan=3,sticky="NSEw")
        sv_ttk.set_theme("dark")
        self.mainscreen.mainloop()
    
    def infilebrowse(self):
        self.infile=tk.filedialog.askopenfilename()
        self.infilepathlabel.configure(text=self.infile)
    def outfilebrowse(self):
        self.outfile=tk.filedialog.asksaveasfilename()
        self.outfilepathlabel.configure(text=self.outfile)
        
        # if self.infile[-4:]=="xlsx":
        #     try:
        #         data=pd.read_excel(self.infile)
        #         self.mz=data["m/z"].values
        #         self.intens=data["I"].values
        #         self.infilepathlabel.configure(background="olive drab")
        #     except:
        #         self.infilepathlabel.configure(background="red4")
        # else:
        #     self.infilepathlabel.configure(background="red2")
    def convert(self):
        if self.infile==0:
            print("no infile selected")
            return
        if self.outfile==0:
            print("no outfile selected")
            return
        try:
            data=pd.read_excel(self.infile)
            masses=data["theoretical mass"].values
            labels=data["label"].values
            charges=data["charge"].values
        except:
            print("error opening file")
            return
        try:
            with open(self.outfile+".ref","w") as f:
                for i in range(len(labels)):
                    labels[i]=labels[i].replace(" ","")
                    f.write(f"{labels[i]}\t{masses[i]}\t{charges[i]}\n")
                f.close
        except:
            print("error in saving")
            
r=gui()