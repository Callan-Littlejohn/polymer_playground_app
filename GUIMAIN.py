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
import pipelines.pipline_1_testing
import matplotlib.pyplot as plot
from matplotlib.gridspec import GridSpec
from functions import elements_functions
from tkinter import ttk
import sv_ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import scipy.signal
import matplotlib
import time
from matplotlib.figure import Figure

class gui():
    def __init__(self):
        self.figsize=[5,5]
        self.make_mainscreen()
        
    def make_mainscreen(self):
        self.mainscreen=tk.Tk()
        self.mainscreen.title("Polymer Playground")
        self.mainscreen.iconbitmap("logo.ico")
        self.mainscreen.resizable(True,True)
        icon=tk.PhotoImage(file="logo.png")
        self.mainscreen.iconphoto(False,icon )
        
        tabcontrol=ttk.Notebook(self.mainscreen)
        self.spectratab=ttk.Frame(tabcontrol)
        tabcontrol.add(self.spectratab,text="Overall_spectra") 
        self.mmdtab=ttk.Frame(tabcontrol)
        tabcontrol.add(self.mmdtab,text="MMD")
        self.mmdprotab=ttk.Frame(tabcontrol)
        tabcontrol.add(self.mmdprotab,text="MMD projection")
        tabcontrol.grid(row=2,column=0,rowspan=15,columnspan=15,sticky="NESW")  
        
        tabcontrol2=ttk.Notebook(self.mainscreen)
        self.fileiotab=ttk.Frame(tabcontrol2)
        tabcontrol2.add(self.fileiotab,text="File IO") 
        self.targettab=ttk.Frame(tabcontrol2)
        tabcontrol2.add(self.targettab,text="Targetting")
        self.savetab=ttk.Frame(tabcontrol2)
        tabcontrol2.add(self.savetab,text="Saving")
        tabcontrol2.grid(row=2,column=16,rowspan=15,columnspan=15,sticky="NESW")  
        
        
        #setup of canvas tabs
        self.spectracanvas=tk.Canvas(self.spectratab,width=500,height=500)
        self.spectracanvas.pack()
        self.spectracanvas.create_image(0,0,anchor="nw",image=icon)
        
        self.mmdcanvas=tk.Canvas(self.mmdtab,width=500,height=500)
        self.mmdcanvas.pack()
        self.mmdcanvas.create_image(0,0,anchor="nw",image=icon)
        
        self.mmdprocanvas=tk.Canvas(self.mmdprotab,width=500,height=500)
        self.mmdprocanvas.pack()
        self.mmdprocanvas.create_image(0,0,anchor="nw",image=icon)
        
        
        # setup of file tab
        infilelabel=ttk.Label(self.fileiotab,text="  Peak List  ").grid(row=0,column=0,columnspan=2,sticky="NSEW")
        self.infilepathlabel=ttk.Label(self.fileiotab,text="  choose peak list  ")
        self.infilepathlabel.grid(row=0,column=3,sticky="NSEW")
        self.infilebrowsebutton=ttk.Button(self.fileiotab,text="browse",command=self.infilebrowse).grid(row=0,column=4, sticky="NSEW")
        
        # outfilelabel=ttk.Label(self.fileiotab,text="  assignment file  ").grid(row=0,column=0,columnspan=2,sticky="NSEW")
        # self.outfilepathlabel=ttk.Label(self.fileiotab,text="  choose peak list  ")
        # self.outfilepathlabel.grid(row=0,column=3,sticky="NSEW")
        # self.outfilebrowsebutton=ttk.Button(self.fileiotab,text="browse",command=self.infilebrowse).grid(row=0,column=4, sticky="NSEW")
        
        
        #setup targetting tab
        monomerlabel=ttk.Label(self.targettab,text="repeat unit").grid(row=0,column=0)
        self.monomerentry=ttk.Entry(self.targettab)
        self.monomerentry.grid(row=0,column=1)
        
        eglabel=ttk.Label(self.targettab,text="expected endgroups").grid(row=1,column=0)
        self.egentry=ttk.Entry(self.targettab)
        self.egentry.grid(row=1,column=1)
        
        adlabel=ttk.Label(self.targettab,text="expected adducts").grid(row=2,column=0)
        self.adentry=ttk.Entry(self.targettab)
        self.adentry.grid(row=2,column=1)
        
        zlabel=ttk.Label(self.targettab,text="max charge").grid(row=3,column=0)
        self.zentry=ttk.Entry(self.targettab)
        self.zentry.grid(row=3,column=1)
        
        
        analysebutt=ttk.Button(self.targettab,text="analyse",command=self.analyse).grid(row=4,column=0,columnspan=2,sticky="NSEW")
        
        
        sv_ttk.set_theme("dark")
        self.mainscreen.mainloop()
        


    def infilebrowse(self):
        self.infile=tk.filedialog.askopenfilename()
        self.infilepathlabel.configure(text=self.infile)
        if self.infile[-4:]=="xlsx":
            try:
                data=pd.read_excel(self.infile)
                mz=data["m/z"].values
                intens=data["I"].values
                self.infilepathlabel.configure(background="olive drab")
            except:
                self.infilepathlabel.configure(background="red4")
        else:
            self.infilepathlabel.configure(background="red2")
        
        
    def analyse(self):
        expected_adducts=self.adentry.get().split(",")
        expected_endgroups=self.egentry.get().split(",")
        repeat_unit=self.monomerentry.get()
        charge=int(self.zentry.get())
        mmd,rmmd, spectra, mmdpeaks,assignments=pipelines.pipline_1_testing.pipeline1(self.infile,expected_endgroups,expected_adducts,repeat_unit,charge)
        slabels, smasses,sints=sort_by_labels(assignments, expected_adducts, expected_endgroups)#print(assignments)
        mmdpeaksx,mmdpeaksy=line_up_mmd_peaks(mmd[0],mmd[1],smasses, repeat_unit)
        size=1
        ss1=np.subtract(spectra[1],min(spectra[1]))
        ss2=np.divide(ss1,max(ss1))
        for widget in self.spectratab.winfo_children():
            widget.destroy()
        for widget in self.mmdtab.winfo_children():
            widget.destroy()
        for widget in self.mmdprotab.winfo_children():
            widget.destroy()
        
        
        self.specfig=Figure(figsize=self.figsize)
        self.axspec=self.specfig.add_subplot(111)
        self.mmdfig=Figure(figsize=self.figsize)
        self.axmmd=self.mmdfig.add_subplot(111)
        self.mmdprofig=Figure(figsize=self.figsize)
        self.axmmdpro=self.mmdprofig.add_subplot(111)
        
        markerline, stemlines, baseline =self.axspec.stem(spectra[0],spectra[1],markerfmt=" ",linefmt="k",basefmt="k")#np.argmax(yaxis)
        plot.setp(stemlines, 'linewidth', 0.3)
        
        self.axspec.scatter(0,1000,s=0.001)
        
        self.axmmd.scatter(spectra[0],spectra[2],s=ss2)
        
        self.axmmdpro.plot(mmd[0],mmd[1],linewidth=0.5,c="k")
        self.axmmdpro.scatter(0,0,s=0.001)
        
        for i in range(len(slabels)):
            for j in range(len(sints[i])):
                sints[i][j]=float(sints[i][j])
                #'print(sints)
                smasses[i][j]=float(smasses[i][j])
            #print(slabels)
            try:
                label="+"+slabels[i][0].split(" ")[1]
                self.axspec.scatter(smasses[i],np.add(sints[i],(max(spectra[1])*1.005)-max(spectra[1])),s=size*3,label=label)
                #ss=np.divide(np.log10(sints[i]),9)
                ss=np.divide(np.subtract(sints[i],min(spectra[1])),max(ss1))
                #print(ss)
                self.axmmd.scatter(smasses[i],np.mod(smasses[i],elements_functions.formula_to_mass(repeat_unit)),s=np.multiply(ss,20))
                #print(mmdpeaksx[i])
                self.axmmdpro.scatter(mmdpeaksx[i],np.add(mmdpeaksy[i],0.8),label=label,s=size*50)#
                
            except:
                x=1
        
        
        self.axspec.legend()
        self.axspec.spines['left'].set_visible(True)
        self.axspec.spines['bottom'].set_visible(True)
        self.axspec.spines['top'].set_visible(False)
        self.axspec.spines['right'].set_visible(False)
        
        self.axmmd.spines['left'].set_visible(True)
        self.axmmd.spines['bottom'].set_visible(True)
        self.axmmd.spines['top'].set_visible(False)
        self.axmmd.spines['right'].set_visible(False)
        
        self.axmmdpro.spines['left'].set_visible(True)
        self.axmmdpro.spines['bottom'].set_visible(True)
        self.axmmdpro.spines['top'].set_visible(False)
        self.axmmdpro.spines['right'].set_visible(False)
        # fig1, ax1 = plot.subplots(figsize=(5, 5), dpi=100)
        # markerline, stemlines, baseline = ax1.stem(spectra[0], spectra[1], markerfmt=" ", linefmt="k", basefmt="k")
        # plot.setp(stemlines, 'linewidth', 0.3)
        # ax1.scatter(1000, 0, s=0.001)
        # for i in range(len(slabels)):
        #     label = "+" + slabels[i][0].split(" ")[1]
        #     ax1.scatter(smasses[i], [s + (max(spectra[1]) * 1.005) - max(spectra[1]) for s in sints[i]], s=size*3, label=label)
        # ax1.legend()
        #ax1.set_title("Overall Spectra")
        self.spectracanvas=FigureCanvasTkAgg(self.specfig,master=self.spectratab)
        self.spectracanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        self.spectratoolbar=NavigationToolbar2Tk(self.spectracanvas,self.spectratab)
        self.spectracanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        
        self.mmdcanvas=FigureCanvasTkAgg(self.mmdfig,master=self.mmdtab)
        self.mmdcanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        self.mmdtoolbar=NavigationToolbar2Tk(self.mmdcanvas,self.mmdtab)
        self.mmdcanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        
        self.mmdprocanvas=FigureCanvasTkAgg(self.mmdprofig,master=self.mmdprotab)
        self.mmdprocanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        self.mmdprotoolbar=NavigationToolbar2Tk(self.mmdprocanvas,self.mmdprotab)
        self.mmdprocanvas.get_tk_widget().pack(fill=tk.BOTH,expand=True)
        
def sort_by_labels(assignments,expected_adducts,expected_endgroups):
    assignmentst=np.transpose(assignments)
    obsmasses=assignmentst[0]
    intensitys=assignmentst[1]
    theomasses=assignmentst[2]
    labels=assignmentst[3]
    mes=assignmentst[4]
    sorted_labels=[]
    sorted_masses=[]
    sorted_intensitys=[]
    for j in expected_endgroups:
        sorted_labels.append([])
        sorted_masses.append([])
        sorted_intensitys.append([])
    for i in range(len(labels)):
        for j in range(len(expected_endgroups)):
            endgroup_ID=" "+expected_endgroups[j]+" "
            if endgroup_ID in labels[i]:
                sorted_labels[j].append(labels[i])
                sorted_masses[j].append(obsmasses[i])
                sorted_intensitys[j].append(intensitys[i])
    return sorted_labels,sorted_masses,sorted_intensitys

def line_up_mmd_peaks(x,y,sorted_peaks,repeat_unit):
    peaks,_=scipy.signal.find_peaks(y,height=3)
    peak_center=[]
    peak_height=[]
    mmdlocs=[]
    mmdheights=[]
    refmass=elements_functions.formula_to_mass(repeat_unit)
    for i in peaks:
        peak_center.append(x[i])
        peak_height.append(y[i])
    for i in sorted_peaks:
        tempx,tempy=[],[]
        for j in i:
            mmd=float(j)%refmass
            checker=np.abs(np.subtract(peak_center,mmd))
            closest=np.argmin(checker)
            #print(peak_center[closest],mmd)
            tempx.append(peak_center[closest])
            tempy.append(peak_height[closest])
        mmdlocs.append(tempx)
        mmdheights.append(tempy)
    return mmdlocs,mmdheights
a=gui()