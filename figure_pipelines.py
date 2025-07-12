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
import pipelines.pipline_1_testing
import matplotlib.pyplot as plot
from matplotlib.gridspec import GridSpec
from functions import elements_functions
import tkinter as tk
import scipy.signal
import matplotlib
import time

font = {'family' : 'normal','weight' : 'normal'}  
plot.rc("font")

def makefigs(infile, expected_adducts,expected_endgroups,repeat_unit,charge,p2=False):
    if p2==False:
        mmd,rmmd, spectra, mmdpeaks,assignments=pipelines.pipline_1_testing.pipeline1(infile,expected_endgroups,expected_adducts,repeat_unit,charge)
    else:
       mmd,rmmd, spectra, mmdpeaks,assignments=pipelines.pipline_1_testing.pipeline2(infile,expected_endgroups,expected_adducts,repeat_unit,charge) 
    #print(mmdpeaks[0])
    # m,f,e=elements_functions.exact_mass_calculator(mmdpeaks[0], expected_adducts, [], charge, repeat_unit, 1000,0)
    # print(np.transpose([m,f,e]))
    # print(len(m),len(mmdpeaks[0]))
    # return
    slabels, smasses,sints=sort_by_labels(assignments, expected_adducts, expected_endgroups)#print(assignments)
    mmdpeaksx,mmdpeaksy=line_up_mmd_peaks(mmd[0],mmd[1],smasses, repeat_unit)
    #print(mmdpeaksy)
    print(len(assignments))
    #print(assignments)
    mainfig=plot.figure(constrained_layout=True,dpi=300,figsize=[10,10])
    gs=GridSpec(2,2,figure=mainfig)
    
    spectrax=mainfig.add_subplot(gs[0,:])
    mmdax=mainfig.add_subplot(gs[1,0])
    mmdproax=mainfig.add_subplot(gs[1,1])
    line_up_mmd_peaks(mmd[0],mmd[1],smasses,repeat_unit)
    size=1
    markerline, stemlines, baseline =spectrax.stem(spectra[0],spectra[1],markerfmt=" ",linefmt="k",basefmt="k")#np.argmax(yaxis)
    plot.setp(stemlines, 'linewidth', 0.3)
    spectrax.scatter(1000,0,s=0.001)
    ss1=np.subtract(spectra[1],min(spectra[1]))
    ss2=np.divide(ss1,max(ss1))
    #ss1=np.divide(np.log10(spectra[1]),9)
    mmdax.scatter(spectra[0],spectra[2],s=ss2)
    #print(ss1)
    mmdproax.plot(mmd[0],mmd[1],linewidth=0.5,c="k")
    mmdproax.scatter(0,0,s=0.001)
    #mmdproax.scatter(mmdpeaks[0],mmdpeaks[2])
    #mmdproax.plot(rmmd[0],rmmd[1],linewidth=0.5,c="r")
    for i in range(len(slabels)):
        for j in range(len(sints[i])):
            sints[i][j]=float(sints[i][j])
            #'print(sints)
            smasses[i][j]=float(smasses[i][j])
        #print(slabels)
        try:
            label="+"+slabels[i][0].split(" ")[1]
            spectrax.scatter(smasses[i],np.add(sints[i],(max(spectra[1])*1.005)-max(spectra[1])),s=size*3,label=label)
            #ss=np.divide(np.log10(sints[i]),9)
            ss=np.divide(np.subtract(sints[i],min(spectra[1])),max(ss1))
            #print(ss)
            mmdax.scatter(smasses[i],np.mod(smasses[i],elements_functions.formula_to_mass(repeat_unit)),s=np.multiply(ss,20))
            #print(mmdpeaksx[i])
            mmdproax.scatter(mmdpeaksx[i],np.add(mmdpeaksy[i],0.8),label=label,s=size*50)#
            
        except:
            x=1
    spectrax.legend()
    spectrax.spines['left'].set_visible(True)
    spectrax.spines['bottom'].set_visible(True)
    spectrax.spines['top'].set_visible(False)
    spectrax.spines['right'].set_visible(False)
    
    mmdax.spines['left'].set_visible(True)
    mmdax.spines['bottom'].set_visible(True)
    mmdax.spines['top'].set_visible(False)
    mmdax.spines['right'].set_visible(False)
    
    mmdproax.spines['left'].set_visible(True)
    mmdproax.spines['bottom'].set_visible(True)
    mmdproax.spines['top'].set_visible(False)
    mmdproax.spines['right'].set_visible(False)
    return assignments

def makefigs_topone(infile, expected_adducts,expected_endgroups,repeat_unit,charge,p2=False):
    if p2==False:
        mmd,rmmd, spectra, mmdpeaks,assignments=pipelines.pipline_1_testing.pipeline1(infile,expected_endgroups,expected_adducts,repeat_unit,charge)
    else:
       mmd,rmmd, spectra, mmdpeaks,assignments=pipelines.pipline_1_testing.pipeline2(infile,expected_endgroups,expected_adducts,repeat_unit,charge) 
    #print(mmdpeaks[0])
    # m,f,e=elements_functions.exact_mass_calculator(mmdpeaks[0], expected_adducts, [], charge, repeat_unit, 1000,0)
    # print(np.transpose([m,f,e]))
    # print(len(m),len(mmdpeaks[0]))
    # return
    
    slabels, smasses,sints=sort_by_labels(assignments, expected_adducts, expected_endgroups)#print(assignments)
    mmdpeaksx,mmdpeaksy=line_up_mmd_peaks(mmd[0],mmd[1],smasses, repeat_unit)
    #print(mmdpeaksy)
    print(len(assignments))
    #print(assignments)
    #mainfig=plot.figure(constrained_layout=True,dpi=300,figsize=[10,10])
    # gs=GridSpec(2,2,figure=mainfig)
    
    # spectrax=mainfig.add_subplot(gs[0,:])
    # mmdax=mainfig.add_subplot(gs[1,0])
    # mmdproax=mainfig.add_subplot(gs[1,1])
    line_up_mmd_peaks(mmd[0],mmd[1],smasses,repeat_unit)
    size=1
    plot.figure(figsize=[5,5],dpi=300)
    markerline, stemlines, baseline =plot.stem(spectra[0],spectra[1],markerfmt=" ",linefmt="k",basefmt="k")#np.argmax(yaxis)
    plot.setp(stemlines, 'linewidth', 0.3)
    plot.scatter(1000,0,s=0.001)
    ss1=np.subtract(spectra[1],min(spectra[1]))
    ss2=np.divide(ss1,max(ss1))
    # #ss1=np.divide(np.log10(spectra[1]),9)
    # mmdax.scatter(spectra[0],spectra[2],s=ss2)
    # #print(ss1)
    # mmdproax.plot(mmd[0],mmd[1],linewidth=0.5,c="k")
    # mmdproax.scatter(0,0,s=0.001)
    #mmdproax.scatter(mmdpeaks[0],mmdpeaks[2])
    #mmdproax.plot(rmmd[0],rmmd[1],linewidth=0.5,c="r")
    for i in range(len(slabels)):
        for j in range(len(sints[i])):
            sints[i][j]=float(sints[i][j])
            #'print(sints)
            smasses[i][j]=float(smasses[i][j])
        #print(slabels)
        try:
            label="+"+slabels[i][0].split(" ")[1]
            plot.scatter(smasses[i],np.add(sints[i],(max(spectra[1])*1.005)-max(spectra[1])),s=size*10,label=label)
            #ss=np.divide(np.log10(sints[i]),9)
            # ss=np.divide(np.subtract(sints[i],min(spectra[1])),max(ss1))
            # #print(ss)
            # mmdax.scatter(smasses[i],np.mod(smasses[i],elements_functions.formula_to_mass(repeat_unit)),s=np.multiply(ss,20))
            # #print(mmdpeaksx[i])
            # mmdproax.scatter(mmdpeaksx[i],np.add(mmdpeaksy[i],0.8),label=label,s=size*50)#
            
        except:
            x=1
    #plot.legend()
    plot.gca().spines['top'].set_visible(False)
    plot.gca().spines['right'].set_visible(False)
    plot.gca().spines['bottom'].set_visible(True)
    plot.gca().spines['left'].set_visible(True)
    # spectrax.legend()
    # spectrax.spines['left'].set_visible(True)
    # spectrax.spines['bottom'].set_visible(True)
    # spectrax.spines['top'].set_visible(False)
    # spectrax.spines['right'].set_visible(False)
    
    # mmdax.spines['left'].set_visible(True)
    # mmdax.spines['bottom'].set_visible(True)
    # mmdax.spines['top'].set_visible(False)
    # mmdax.spines['right'].set_visible(False)
    
    # mmdproax.spines['left'].set_visible(True)
    # mmdproax.spines['bottom'].set_visible(True)
    # mmdproax.spines['top'].set_visible(False)
    # mmdproax.spines['right'].set_visible(False)
    return assignments

def makefigs_mmdpro(infile, expected_adducts,expected_endgroups,repeat_unit,charge,p2=False):
    if p2==False:
        mmd,rmmd, spectra, mmdpeaks,assignments=pipelines.pipline_1_testing.pipeline1(infile,expected_endgroups,expected_adducts,repeat_unit,charge)
    else:
       mmd,rmmd, spectra, mmdpeaks,assignments=pipelines.pipline_1_testing.pipeline2(infile,expected_endgroups,expected_adducts,repeat_unit,charge) 
    #print(mmdpeaks[0])
    # m,f,e=elements_functions.exact_mass_calculator(mmdpeaks[0], expected_adducts, [], charge, repeat_unit, 1000,0)
    # print(np.transpose([m,f,e]))
    # print(len(m),len(mmdpeaks[0]))
    # return
    
    slabels, smasses,sints=sort_by_labels(assignments, expected_adducts, expected_endgroups)#print(assignments)
    mmdpeaksx,mmdpeaksy=line_up_mmd_peaks(mmd[0],mmd[1],smasses, repeat_unit)
    #print(mmdpeaksy)
    #print(assignments)
    #mainfig=plot.figure(constrained_layout=True,dpi=300,figsize=[10,10])
    # gs=GridSpec(2,2,figure=mainfig)
    
    # spectrax=mainfig.add_subplot(gs[0,:])
    # mmdax=mainfig.add_subplot(gs[1,0])
    # mmdproax=mainfig.add_subplot(gs[1,1])
    line_up_mmd_peaks(mmd[0],mmd[1],smasses,repeat_unit)
    size=1
    plot.figure(figsize=[5,5],dpi=300)
    #markerline, stemlines, baseline =plot.stem(spectra[0],spectra[1],markerfmt=" ",linefmt="k",basefmt="k")#np.argmax(yaxis)
   # plot.setp(stemlines, 'linewidth', 0.3)
    #plot.scatter(1000,0,s=0.001)
    ss1=np.subtract(spectra[1],min(spectra[1]))
    ss2=np.divide(ss1,max(ss1))
    # #ss1=np.divide(np.log10(spectra[1]),9)
    # mmdax.scatter(spectra[0],spectra[2],s=ss2)
    # #print(ss1)
    plot.plot(mmd[0],mmd[1],linewidth=0.5,c="k")
    plot.scatter(0,0,s=0.001)
    #scatter(mmdpeaks[0],mmdpeaks[2])
    #plot.plot(rmmd[0],rmmd[1],linewidth=0.5,c="r")
    for i in range(len(slabels)):
        for j in range(len(sints[i])):
            sints[i][j]=float(sints[i][j])
            #'print(sints)
            smasses[i][j]=float(smasses[i][j])
        #print(slabels)
        try:
            label="+"+slabels[i][0].split(" ")[1]
            #plot.scatter(smasses[i],np.add(sints[i],(max(spectra[1])*1.005)-max(spectra[1])),s=size*3,label=label)
            #ss=np.divide(np.log10(sints[i]),9)
            # ss=np.divide(np.subtract(sints[i],min(spectra[1])),max(ss1))
            # #print(ss)
            # mmdax.scatter(smasses[i],np.mod(smasses[i],elements_functions.formula_to_mass(repeat_unit)),s=np.multiply(ss,20))
            # #print(mmdpeaksx[i])
            plot.scatter(mmdpeaksx[i],np.add(mmdpeaksy[i],0.8),label=label,s=size*50)#
            
        except:
            x=1
    #plot.legend(loc="upper right")
    plot.gca().spines['top'].set_visible(False)
    plot.gca().spines['right'].set_visible(False)
    plot.gca().spines['bottom'].set_visible(True)
    plot.gca().spines['left'].set_visible(True)
    # spectrax.legend()
    # spectrax.spines['left'].set_visible(True)
    # spectrax.spines['bottom'].set_visible(True)
    # spectrax.spines['top'].set_visible(False)
    # spectrax.spines['right'].set_visible(False)
    
    # mmdax.spines['left'].set_visible(True)
    # mmdax.spines['bottom'].set_visible(True)
    # mmdax.spines['top'].set_visible(False)
    # mmdax.spines['right'].set_visible(False)
    
    # mmdproax.spines['left'].set_visible(True)
    # mmdproax.spines['bottom'].set_visible(True)
    # mmdproax.spines['top'].set_visible(False)
    # mmdproax.spines['right'].set_visible(False)
    return assignments
    
def makefigs_withnums(infile, expected_adducts,expected_endgroups,repeat_unit,charge,p2=False):
    if p2==False:
        mmd,rmmd, spectra, mmdpeaks,assignments=pipelines.pipline_1_testing.pipeline1(infile,expected_endgroups,expected_adducts,repeat_unit,charge)
    else:
       mmd,rmmd, spectra, mmdpeaks,assignments=pipelines.pipline_1_testing.pipeline2(infile,expected_endgroups,expected_adducts,repeat_unit,charge) 
    #print(mmdpeaks[0])
    # m,f,e=elements_functions.exact_mass_calculator(mmdpeaks[0], expected_adducts, [], charge, repeat_unit, 1000,0)
    # print(np.transpose([m,f,e]))
    # print(len(m),len(mmdpeaks[0]))
    # return
    slabels, smasses,sints=sort_by_labels(assignments, expected_adducts, expected_endgroups)#print(assignments)
    mmdpeaksx,mmdpeaksy=line_up_mmd_peaks(mmd[0],mmd[1],smasses, repeat_unit)
    #print(mmdpeaksy)
    print(len(assignments))
    #print(assignments)
    mainfig=plot.figure(constrained_layout=True,dpi=300,figsize=[10,10])
    gs=GridSpec(2,2,figure=mainfig)
    
    spectrax=mainfig.add_subplot(gs[0,:])
    mmdax=mainfig.add_subplot(gs[1,0])
    mmdproax=mainfig.add_subplot(gs[1,1])
    line_up_mmd_peaks(mmd[0],mmd[1],smasses,repeat_unit)
    size=1
    markerline, stemlines, baseline =spectrax.stem(spectra[0],spectra[1],markerfmt=" ",linefmt="k",basefmt="k")#np.argmax(yaxis)
    plot.setp(stemlines, 'linewidth', 0.3)
    spectrax.scatter(1000,0,s=0.001)
    ss1=np.subtract(spectra[1],min(spectra[1]))
    ss2=np.divide(ss1,max(ss1))
    #ss1=np.divide(np.log10(spectra[1]),9)
    mmdax.scatter(spectra[0],spectra[2],s=ss2,c="k")
    #print(ss1)
    mmdproax.plot(mmd[0],mmd[1],linewidth=0.5,c="k")
    mmdproax.scatter(0,0,s=0.001)
    #mmdproax.scatter(mmdpeaks[0],mmdpeaks[2])
    #mmdproax.plot(rmmd[0],rmmd[1],linewidth=0.5,c="r")
    count=1
    for i in range(len(slabels)):
        for j in range(len(sints[i])):
            sints[i][j]=float(sints[i][j])
            #'print(sints)
            smasses[i][j]=float(smasses[i][j])
        #print(slabels)
        if len(slabels[i])>1:#try:
            label="+"+slabels[i][0].split(" ")[1]
            spectrax.scatter(smasses[i],np.add(sints[i],(max(spectra[1])*1.05)-max(spectra[1])),marker="${}$".format(str(count)),s=40,label=label,color="k")
            #spectrax.text(smasses[i],np.add(sints[i],(max(spectra[1])*1.005)-max(spectra[1])),str(count),color="k", fontsize=12,label=label)
            #ss=np.divide(np.log10(sints[i]),9)
            ss=np.divide(np.subtract(sints[i],min(spectra[1])),max(ss1))
            #print(ss)
            mmdax.scatter(smasses[i],np.mod(smasses[i],elements_functions.formula_to_mass(repeat_unit)),s=np.multiply(ss,20),c="k")
            mmdax.scatter(np.multiply(np.divide(smasses[i],smasses[i]),max(spectra[0])),np.round(np.mod(smasses[i],elements_functions.formula_to_mass(repeat_unit)),2),marker="${}$".format(str(count)),s=40,label=label,color="k")
            mmdax.hlines(np.round(np.mod(smasses[i],elements_functions.formula_to_mass(repeat_unit)),2),0,max(spectra[0]),colors="r",linestyles="dotted")
            #print(mmdpeaksx[i])
            mmdproax.scatter(mmdpeaksx[i],np.add(mmdpeaksy[i],0.8),label=label,marker="${}$".format(str(count)),s=40,color="k")#
            count+=1
            
        # except:
        #     x=1
    spectrax.legend()
    spectrax.set_ylim(0,(max(spectra[1])*1.8))
    spectrax.set_xlim(0,3000)
    spectrax.spines['left'].set_visible(True)
    spectrax.spines['bottom'].set_visible(True)
    spectrax.spines['top'].set_visible(False)
    spectrax.spines['right'].set_visible(False)
    
    mmdax.spines['left'].set_visible(True)
    mmdax.spines['bottom'].set_visible(True)
    mmdax.spines['top'].set_visible(False)
    mmdax.spines['right'].set_visible(False)
    
    mmdproax.spines['left'].set_visible(True)
    mmdproax.spines['bottom'].set_visible(True)
    mmdproax.spines['top'].set_visible(False)
    mmdproax.spines['right'].set_visible(False)
    return assignments
    
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
        
def make_assignment_table(outfile, assignments):
    ass=np.transpose(assignments)
    data={"label":ass[3],"theoretical mass": ass[2],"charge":ass[5],"observed mass":ass[0],"intensity":ass[1],"mass error":ass[4]}
    df=pd.DataFrame(data)
    df.to_excel(outfile)
    
  


# expected_adducts=["Na","H","NH4","K","Ag","Ca","Li"]
# expected_endgroups=["H2O","","C2H2O","CH3","CH3O","OH","O","C","C2","H2","C2H2","H"]
# ass=makefigs("C:/Users/FTICR_Kool_Kidz_PC/Downloads/Peg_41k_CID.xlsx",expected_adducts,expected_endgroups,"C2H4O",4,p2=False)

# expected_adducts=["H"]
# expected_endgroups=["H3O","H2O","CH5O","CH3O","OH","CH3O2"]
# ass=makefigs_topone("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PEG_Dart.xlsx",expected_adducts,expected_endgroups,"C2H4O",1,p2=False)
# make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PEG_41K(Assigned).xlsx", ass)
# expected_adducts=list(set(["Na","H","NH4","K"]))#"Ag","Ca","Li"]
# expected_endgroups=["","C4H2O","H","C8H6O"]
# makefigs_topone("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PMAA312k.xlsx",expected_adducts,expected_endgroups,"C2H4ONa",10,p2=False)
# plot.xlim(0,2000)
# expected_adducts=["r"]
# expected_endgroups=["C4H10","C7H8","C7H5O","C7H7O","C7H5O2","C5H11O","C4H9O","C8H9O","C8H5","C15H10O","C7H3O2","C12H22N","C8H15","C6H11","C12H16O"]
# makefigs("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_2006hrs.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=True)

# expected_adducts=["r"]
# expected_endgroups=["C4H10","C8H15","C7H8","C14H12"]
# makefigs_topone("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_noDegrade.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=False)#"C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_iso994_CID.xlsx""C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_noDegrade.xlsx"
# plot.show()

# expected_adducts=["r"]
# expected_endgroups=["C4H10","C8H15","C7H8","C14H12"]
# ass=makefigs_mmdpro("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_noDegrade.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=False)#"C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_iso994_CID.xlsx""C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_noDegrade.xlsx"
# #plot.legend()
# make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_noDegrade(Assigned).xlsx", ass)

# expected_adducts=["r"]
# expected_endgroups=["C4H10","C2H3","C7H8"]
# ass=makefigs_mmdpro("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_iso994_CID.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=False)#"C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_iso994_CID.xlsx""C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_noDegrade.xlsx"
# make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_iso994_CID(Assigned).xlsx", ass)
# # expected_adducts=["r"]
# expected_endgroups=["C4H10","C4H10O","C7H8","C9H9","C8H15","C15H13","C11H21","C14H17","C11H13"]
# makefigs("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_40uv_degraded2.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=True)
# expected_adducts=["r"]
# expected_endgroups=["C4H10","C8H15","C7H8","C14H12"]
# makefigs("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_iso994_CID.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=True)#"C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_iso994_CID.xlsx""C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_noDegrade.xlsx"

# expected_adducts=["r"]
# expected_endgroups=list(set(["C4H10","C7H8","C7H5O","C7H7O","C7H5O2","C5H11O","C4H9O","C8H9O","C8H5","C15H10O","C7H3O2","C8H15","C6H11","C12H16O","C4H10O","C7H8","C9H9","C8H15","C15H13","C11H21","C14H17","C11H13"]))
# ass=makefigs("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_noDegrade.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=False)
# make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PSnodegrade(Assigned).xlsx", ass)
                      
# expected_adducts=["r"]
# expected_endgroups=list(set(["C4H10","C7H8","C7H5O","C7H7O","C7H5O2","C5H11O","C4H9O","C8H9O","C8H5","C15H10O","C7H3O2","C12H22N","C8H15","C6H11","C12H16O","C4H10O","C7H8","C9H9","C8H15","C15H13","C11H21","C14H17","C11H13"]))
# ass=makefigs("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS100c2hrs.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=False)
# make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS100_2hrs(Assigned).xlsx", ass)

# expected_adducts=["r"]
# expected_endgroups=list(set(["C4H10","C7H8","C7H5O","C7H7O","C7H5O2","C5H11O","C4H9O","C8H9O","C8H5","C15H10O","C7H3O2","C12H22N","C8H15","C6H11","C12H16O","C4H10O","C7H8","C9H9","C8H15","C15H13","C11H21","C14H17","C11H13"]))
# ass=makefigs("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS1004hr.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=False)
# make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS100_4hrs(Assigned).xlsx", ass)

# expected_adducts=["r"]
# expected_endgroups=list(set(["C4H10","C7H8","C7H5O","C7H7O","C7H5O2","C5H11O","C4H9O","C8H9O","C8H5","C15H10O","C7H3O2","C12H22N","C8H15","C6H11","C12H16O","C4H10O","C7H8","C9H9","C8H15","C15H13","C11H21","C14H17","C11H13"]))
# ass=makefigs("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_1504hr.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=False)
# make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS150_4hrs(Assigned).xlsx", ass)


# expected_adducts=[["Na","H","NH4","K","Ag","Ca","Li"]]
# expected_endgroups=list(set("H2O",""))#list(set(["C4H10","C7H8","C7H5O","C7H7O","C7H5O2","C5H11O","C4H9O","C8H9O","C7H3O2","C12H22N","C8H15","C6H11","C12H16O","C4H10O","C7H8","C9H9","C8H15","C15H13","C11H21","C14H17","C11H13"]))
# ass=makefigs("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_2006hrs.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=False)
# #make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS200_6hrs(Assigned).xlsx", ass)

# expected_adducts=["r"]
# expected_endgroups=list(set(["C4H10","C7H8","C7H5O","C7H7O","C7H5O2","C5H11O","C4H9O","C8H9O","C7H3O2","C12H22N","C8H15","C6H11","C12H16O","C4H10O","C7H8","C9H9","C8H15","C15H13","C11H21","C14H17","C11H13"]))
# ass=makefigs_withnums("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_2006hrs.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=False)
# #make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS200_6hrs(Assigned).xlsx", ass)

# expected_adducts=["r"]
# expected_endgroups=list(set(["C4H10","C7H8","C7H5O","C7H7O","C7H5O2","C5H11O","C4H9O","C8H9O","C8H5","C15H10O","C7H3O2","C12H22N","C8H15","C6H11","C12H16O","C4H10O","C7H8","C9H9","C8H15","C15H13","C11H21","C14H17","C11H13"]))
# ass=makefigs("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_2p5_degraded.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=False)
# make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_2p5(Assigned).xlsx", ass)

# expected_adducts=["r"]
# expected_endgroups=list(set(["C4H10","C7H8","C7H5O","C7H7O","C7H5O2","C5H11O","C4H9O","C8H9O","C8H5","C15H10O","C7H3O2","C12H22N","C8H15","C6H11","C12H16O","C4H10O","C7H8","C9H9","C8H15","C15H13","C11H21","C14H17","C11H13"]))
# ass=makefigs("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_24_degraded.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=False)
# make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_24(Assigned).xlsx", ass)

# expected_adducts=["r"]
# expected_endgroups=list(set(["C4H10","C7H8","C7H5O","C7H7O","C7H5O2","C5H11O","C4H9O","C8H9O","C8H5","C15H10O","C7H3O2","C12H22N","C8H15","C6H11","C12H16O","C4H10O","C7H8","C9H9","C8H15","C15H13","C11H21","C14H17","C11H13"]))
# ass=makefigs("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_40uv_degraded2.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=False)
# make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_40(Assigned).xlsx", ass)

# expected_adducts=["r"]
# expected_endgroups=list(set(["C4H10","C7H8","C7H5O","C7H7O","C7H5O2","C5H11O","C4H9O","C8H9O","C8H5","C15H10O","C7H3O2","C12H22N","C8H15","C6H11","C12H16O","C4H10O","C7H8","C9H9","C8H15","C15H13","C11H21","C14H17","C11H13"]))
# ass=makefigs("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_TGA_50_degraded.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=False)
# make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_501hrs(Assigned).xlsx", ass)

# expected_adducts=["Na","H","NH4","K","Ag"]#,"Ca","Li"]
# expected_endgroups=list(set(["C4H10","C7H8","C7H5O","C7H7O","C7H5O2","C5H11O","C4H9O","C8H9O","C8H5","C15H10O","C7H3O2","C12H22N","C8H15","C6H11","C12H16O","C4H10O","C7H8","C9H9","C8H15","C15H13","C11H21","C14H17","C11H13"]))
# ass=makefigs("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_maldi_sodiated.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=False)
# make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_maldisodiated(Assigned).xlsx", ass)

# expected_adducts=["Na","H","NH4","K","Ag"]#,"Ca","Li"]
# expected_endgroups=list(set(["C4H10","C7H8","C7H5O","C7H7O","C7H5O2","C5H11O","C4H9O","C8H9O","C8H5","C15H10O","C7H3O2","C12H22N","C8H15","C6H11","C12H16O","C4H10O","C7H8","C9H9","C8H15","C15H13","C11H21","C14H17","C11H13"]))
# ass=makefigs("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_MALDI_Ag.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=False)
# make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_maldag(Assigned).xlsx", ass)

# expected_adducts=["Na","H","NH4","K","Ag"]#,"Ca","Li"]
# expected_endgroups=["H2O","",]
# ass=makefigs_topone("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PEG_1K.xlsx",expected_adducts,expected_endgroups,"C2H4O",4,p2=False)
# plot.xlim(0,2200)
# plot.show()
# expected_adducts=["Na","H","NH4","K","Ag"]#,"Ca","Li"]
# expected_endgroups=["H2O","",]
# ass=makefigs_mmdpro("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PEG_1K.xlsx",expected_adducts,expected_endgroups,"C2H4O",4,p2=False)
# make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PEG_1K(Assigned).xlsx", ass)

# # # #print(ass)
# make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PEG_1K(Assigned).xlsx", ass)
# expected_adducts=["Na","H","NH4","K","Ag"]#,"Ca","Li"]
# expected_endgroups=["CH4O","H2O","","C2H2O","CH3","CH3O","OH","O","C","C2","H2","C2H2","H"]
# ass=makefigs_topone("C:/Users/FTICR_Kool_Kidz_PC/Downloads/peg2kmethylether.xlsx",expected_adducts,expected_endgroups,"C2H4O",4,p2=False)
#print(ass)
#make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PEG_methylether(Assigned).xlsx", ass)

# expected_adducts=["Na","H","NH4","K","Ag"]#,"Ca","Li"]
# expected_endgroups=["H2O","","CH4O"]
# ass=makefigs_topone("C:/Users/FTICR_Kool_Kidz_PC/Downloads/peg2kmethylether.xlsx",expected_adducts,expected_endgroups,"C2H4O",4,p2=False)
# plot.show()
# expected_adducts=["Na","H","NH4","K","Ag"]#,"Ca","Li"]
# expected_endgroups=["H2O","","CH4O"]
# ass=makefigs_mmdpro("C:/Users/FTICR_Kool_Kidz_PC/Downloads/peg2kmethylether.xlsx",expected_adducts,expected_endgroups,"C2H4O",4,p2=False)
# make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/peg2kmethylether(Assigned).xlsx", ass)


# expected_adducts=["Na","H","NH4","K","Ag"]#,"Ca","Li"]
# expected_endgroups=["H2O","","CH4O"]
# ass=makefigs_topone("C:/Users/FTICR_Kool_Kidz_PC/Downloads/Pegmethylether20k.xlsx",expected_adducts,expected_endgroups,"C2H4O",4,p2=False)

# expected_adducts=["Na","H","NH4","K","Ag"]#,"Ca","Li"]
# expected_endgroups=["H2O",""]#,"CH4O","C2H2O","CH3","CH3O","OH","O","C","C2","H2","C2H2","H"]
# ass=makefigs_topone("C:/Users/FTICR_Kool_Kidz_PC/Downloads/Peg_41k_CID.xlsx",expected_adducts,expected_endgroups,"C2H4O",4,p2=False)

# expected_adducts=["Na","H","NH4","K","Ag"]#,"Ca","Li"]
# expected_endgroups=["H2O",""]#,"CH4O","C2H2O","CH3","CH3O","OH","O","C","C2","H2","C2H2","H"]
# ass=makefigs_topone("C:/Users/FTICR_Kool_Kidz_PC/Downloads/Unknown_PEG.xlsx",expected_adducts,expected_endgroups,"C2H4O",4,p2=False)

# expected_adducts=["Na","H","NH4","K","Ag"]#,"Ca","Li"]
# expected_endgroups=["H2O",""]#,"CH4O","C2H2O","CH3","CH3O","OH","O","C","C2","H2","C2H2","H"]
# ass=makefigs_mmdpro("C:/Users/FTICR_Kool_Kidz_PC/Downloads/Unknown_PEG.xlsx",expected_adducts,expected_endgroups,"C2H4O",4,p2=False)


# expected_adducts=["Na","H","NH4","K","Ag"]#,"Ca","Li"]
# expected_endgroups=["H2O",""]
# ass=makefigs_topone("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PEG_10kCID.xlsx",expected_adducts,expected_endgroups,"C2H4O",4,p2=False)
# #print(ass)

# #make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PEG_21kCID(Assigned).xlsx", ass)

# expected_adducts=["Na","H","NH4","K"]#,"Ca","Li"]
# expected_endgroups=["H2O","","N2H4","CH4O","CH2","C2H2","C","CH6N","C6H7NO","C5H7N",""]
# ass=makefigs_mmdpro("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PVP_21K_new.xlsx",expected_adducts,expected_endgroups,"C6H9NO",4,p2=False)
#expected_adducts=["Na","H","NH4","K"]#,"Ca","Li"]
# expected_endgroups=["H2O",""]
# ass=makefigs_mmdpro("C:/Users/FTICR_Kool_Kidz_PC/Downloads/SB_YN36.xlsx",expected_adducts,expected_endgroups,"C10H20O",4,p2=False)

# expected_adducts=["H"]#,"Ca","Li"]
# expected_endgroups=["H2O","","C10H18"]

# ass=makefigs_topone("C:/Users/FTICR_Kool_Kidz_PC/Documents/SB_YN_Collaboration/mass lists/YN_36_APPI_937CID.xlsx",expected_adducts,expected_endgroups,"C10H20O",1,p2=False)
# expected_adducts=["H","Na","Ca","Li"]
# expected_endgroups=["H2O",""]

# ass=makefigs_mmdpro("C:/Users/FTICR_Kool_Kidz_PC/Downloads/c9_test.xlsx",expected_adducts,expected_endgroups,"C9H18O",2,p2=False)
# make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/C9(Assigned).xlsx", ass)


# expected_adducts=["Na","H","NH4","K","Ag"]#,"Ca","Li"]
# expected_endgroups=["H2O","",]#,"CH4O","C2H2O","CH3","CH3O","OH","O","C","C2","H2","C2H2","H"]
# ass=makefigs_topone("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PLA_Prusa.xlsx",expected_adducts,expected_endgroups,"C3H4O2",4,p2=False)
# #plot.ylim(0,0.5E9)
# print(len(ass))

# expected_adducts=["Na","H","NH4","K","Ag"]#,"Ca","Li"]
# expected_endgroups=["H2O","",]#,"CH4O","C2H2O","CH3","CH3O","OH","O","C","C2","H2","C2H2","H"]
# ass=makefigs_topone("C:/Users/FTICR_Kool_Kidz_PC/Downloads/pla_generic_white.xlsx",expected_adducts,expected_endgroups,"C3H4O2",4,p2=False)
# #plot.ylim(0,0.5E9)
# print(len(ass))

# # expected_adducts=["Na","H","NH4"]#,"Ca","Li"]
# expected_endgroups=["H2","H2O","","CH4",]#,"CH4O","C2H2O","CH3","CH3O","OH","O","C","C2","H2","C2H2","H"]
# # ass=makefigs_topone("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PMAA_1k.xlsx",expected_adducts,expected_endgroups,"C4H6O2",1,p2=False)
# # plot.show()
# expected_adducts=["Na","H","NH4"]#,"K","Ag"]
# t0=time.time()
# ass=makefigs_topone("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PMAA_1k.xlsx",expected_adducts,expected_endgroups,"C4H6O2",3,p2=False)
# t1=time.time()
# print(t1-t0)
# make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PMAA_1k(assigned.xlsx", ass)
#plot.xlim(2,4)
#plot.ylim(0,0.5E9)
# expected_adducts=["Ag","Na"]#,"Ca","Li"]
# expected_endgroups=["C4H10","",]
# ass=makefigs_mmdpro("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_MALDI_Ag.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=False)
# plot.show()
# expected_adducts=["Ag","Na"]#,"Ca","Li"]
# expected_endgroups=["C4H10","",]
# t0=time.time()
# ass=makefigs_topone("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_MALDI_Ag.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=False)
# t1=time.time()
# print(t1-t0)
# # make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_MALDI_Ag(assigned).xlsx", ass)
# expected_adducts=["Na"]#,"Ca","Li"]
# expected_endgroups=["C4H10","",]
# ass=makefigs_mmdpro("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_maldi_sodiated.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=False)
# plot.show()
# expected_adducts=["Ag","Na"]#,"Ca","Li"]
# expected_endgroups=["C4H10","",]
# t0=time.time()
# ass=makefigs_topone("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_maldi_sodiated.xlsx",expected_adducts,expected_endgroups,"C8H8",1,p2=False)
# t1=time.time()
# plot.ylim(0,3E8)
# print(t1-t0)
# make_assignment_table("C:/Users/FTICR_Kool_Kidz_PC/Downloads/PS_MALDI_Na(assigned).xlsx", ass)

# #plot.xlim(60,65)
# print(len(ass))

# expected_adducts=["H","Na","NH4","Ca","Li"]
# expected_endgroups=["H2O","","C10H20O",""]
# t0=time.time()
# ass=makefigs_mmdpro("C:/Users/FTICR_Kool_Kidz_PC/Downloads/jingyi_c11.xlsx",expected_adducts,expected_endgroups,"C11H22O",2,p2=False)
# print(len(ass))
# plot.legend()