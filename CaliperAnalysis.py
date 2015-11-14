from __future__ import division, absolute_import, \
                                    print_function, unicode_literals

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as stats
import xlsxwriter
import Tkinter
import tkFileDialog
import os

pd.options.mode.chained_assignment = None
# -*- coding: utf-8 -*-

def pullstdcurvedata(df, stdcurverow, stdcurveplate, peakname):
    ''' Finds the herceptin curve data and returns it. Assumes that there
    is a duplicate on the next row and it takes up 12 wells/row'''
    
    startwell = stdcurverow.upper() + '1'
    # Find start of herceptin curve
    rowpos = df['Sample Name'] == startwell
    startrowpos = df[rowpos]['Type']=='LM' # every well starts with a LM peak
    startingstdindex= df[rowpos][startrowpos].index.tolist()
    ilocstartindex= startingstdindex[stdcurveplate-1]
    df_afterstart = df.iloc[ilocstartindex:]
    
    # Find end of herceptin curve
    lastwellofstd = df_afterstart['Sample Name'] ==chr(ord(stdcurverow.upper())+1) +'12'
    lastwellindices =  df_afterstart[lastwellofstd].index.tolist()
    ilocendindex = lastwellindices[-1]
    herceptinstdwells = df.iloc[ilocstartindex:(ilocendindex+1)]  
    notherceptinstdwells1 = df.iloc[0:ilocstartindex]
    notherceptinstdwells2 = df.iloc[ilocendindex:] # right index?
    notherceptinstdwells = pd.concat([notherceptinstdwells1, notherceptinstdwells2])
    
    # pull out desired peaks 
    identifiedpeaks = herceptinstdwells['Type'] == peakname
    forcedpeaks = herceptinstdwells['Type'] == peakname + '*'
    allpeaks = identifiedpeaks | forcedpeaks
    df_stdcurve = herceptinstdwells[allpeaks]
    
    return df_stdcurve, notherceptinstdwells
    
def findmissingwells(df, wholeplate):
    """ Will return an array of strings with which wells are missing in the 
    dataframe passed using the Sample Name column. Assuming full rows!!"""
    wells = df['Sample Name']
    missingwells = []
    startingwell = wells.iloc[0]
    startingrow = startingwell[0]
    startingcolumn = startingwell[1:3]
    if startingcolumn != '1':
        for i in np.linspace(1, int(startingcolumn)-1, int(startingcolumn)-1):
            missingwells += [startingrow + str(int(i))]
        
    if wholeplate:
        lastwell = 'H12'
    else:
        lastwell = wells.iloc[-1]
    theoreticalnumwells = 12* (ord(lastwell[0]) - ord(startingrow)+1)
    for j in np.linspace(0, theoreticalnumwells-1, theoreticalnumwells):
        currwellnum = (int(j) % 12)+1
        currwellrow = chr(int(ord(startingrow) + j/12))
        currwell = currwellrow + str(int(currwellnum))
        if not currwell in wells.values:
           missingwells += [currwell] 
    return list(set(missingwells))
      
def findblankwells(df, df_scaff):
    """Will return an array of strings with which wells are blank or missing the desired peak
    in the dataframe compared to the scaffold file."""
    wells = df_scaff['Well Label']
    missingwells = []
    for i in np.linspace(0, len(wells)-1, len(wells)):
        currwell = wells.iloc[i]
        if not currwell in df['Well Label'].values:
           missingwells += [currwell] 
    return list(set(missingwells))    
    
    
def calcavgstd(df_stdcurve):
    """ Calculate the average of the two herceptin standard curve. Return an array
    with the average values. """
    avgcorrareas = []
    stdcorrareas = df_stdcurve['Corr. Area']
    if len(stdcorrareas) > 24:
        print("You have extra peaks in your herceptin curve, wtf")
    else:
        missingwells = findmissingwells(df_stdcurve, False)
        startingwell = df_stdcurve['Sample Name'].iloc[0]   
        row1 = startingwell[0]
        row2 = chr(int(ord(startingwell[0]) + 1))
        row1val = []
        row2val = []
        for i in np.linspace(0, 11, 12):
            currwellnum = int((i % 12))+1
            currwell1 = row1 + str(currwellnum)
            currwell2 = row2 + str(currwellnum)
            
            if currwell1 in missingwells:
                row1val += [np.nan]
            else:
                temp = df_stdcurve[df_stdcurve['Sample Name'] == currwell1]
                row1val += [temp['Corr. Area'].iloc[0]]
                
            if currwell2 in missingwells:
                row2val += [np.nan]
            else:
                temp = df_stdcurve[df_stdcurve['Sample Name'] == currwell2]
                row2val += [temp['Corr. Area'].iloc[0]]    
        corrareas = np.array([row1val, row2val])
        avgcorrareas = stats.nanmean(corrareas)
    return avgcorrareas        

def stdcurvequadfit(df_stdcurve, dispgraph):
    ''' Fits a quadratic curve to the average of the herceptin standard data''' 
    
    x = np.array([1000, 800, 600, 400, 200, 150, 100, 80, 60 ,40, 20, 10])
    avgstdcurve = calcavgstd(df_stdcurve)
    avgstdcurvetrim = avgstdcurve
    counter =0
    for j in np.linspace(0, len(avgstdcurve)-1, len(avgstdcurve)):
        if np.isnan(avgstdcurve[j]):
            x= np.delete(x, j-counter)
            avgstdcurvetrim = np.delete(avgstdcurvetrim, j-counter)
            counter= counter + 1               

    # do polyfit...
    fitparams = np.polyfit(avgstdcurvetrim, x, 2)
    
    if dispgraph:
        fit = avgstdcurvetrim**2 * fitparams[0] + avgstdcurvetrim * fitparams[1] + fitparams[2]

        plt.plot(avgstdcurvetrim, x, 'k.')
        plt.plot(avgstdcurvetrim, fit)
        plt.xlabel('Corr. Area')
        plt.ylabel('Std Protein Concentration')
        plt.show()
        
    return fitparams
    
def pulldata(notstddata, scaff):
    # pull out desired peaks 
    scaffolds = scaff['Scaffold'].unique()
    forcedscaffolds = [scaff + "*" for scaff in scaffolds]
    identifiedpeaks = notstddata['Type'].isin(scaffolds)
    forcedpeaks = notstddata['Type'].isin(forcedscaffolds)
    allpeaks = identifiedpeaks | forcedpeaks
    sampledata = notstddata[allpeaks]
    
    return sampledata

def calculateconc(sampledata, curvefitparam, df_scaff):
    sampcorrarea = sampledata['Corr. Area']
    concentrations = curvefitparam[0]*sampcorrarea**2 + curvefitparam[1] *sampcorrarea + curvefitparam[2]
    conversions = {'IgG': 150000, 'GFP': 28000, 'scFvFc': 100000, 'scfv': 23000, 'stump':75000}
    combineddf_full = pd.merge(sampledata, df_scaff, on=['Well Label'])
    conc =  []
    nMconc = []
    for i in np.linspace(0, len(combineddf_full)-1, len(combineddf_full)):
        scaffnorm = combineddf_full['Scaffold'].iloc[i]
        nMconc += [(concentrations.iloc[i] / conversions[scaffnorm]) * 10**6]
        conc += [concentrations.iloc[i]]
    combineddf_full['[IgG], ug/mL'] = conc
    combineddf_full['nM'] = nMconc

    return combineddf_full

def notdetectedaddition(addwells, df_scaff, finaldf):
    for i in range(len(addwells)):
        sample = addwells[i]
        if len(sample) == 2:
            welllabel= sample[0] + '0' + sample[1]
        else:
            welllabel = sample
        a=df_scaff[df_scaff['Well Label'] == welllabel]
        varID = a['Variant ID'].iloc[0]
        scaff = a['Scaffold'].iloc[0]
        finaldf.loc[i+len(finaldf)] = [welllabel, sample, varID, scaff, scaff, 'Not Detected', 0, 'N/A']
    finaldf=finaldf.sort(['Well Label'])
    return finaldf
    
def loadscaffname(): 
    df_scaff_name.set(tkFileDialog.askopenfilename(filetypes =( ("csv files", ".csv"), ("all files", "*.*"))))
    return df_scaff_name
    
def loadpeakname(): 
    peak_table_name.set(tkFileDialog.askopenfilename(filetypes =( ("csv files", ".csv"), ("all files", "*.*"))))
    return peak_table_name

def choosesavelocation(): 
    savelocation.set(tkFileDialog.askdirectory())
    return savelocation   

def close_window():
    tk.destroy()
    
    
tk = Tkinter.Tk()
frame = Tkinter.Frame(tk, padx=3, pady=4)
frame.grid(column=0, row=0)
tk.title("Caliper Concentration Calculator")
df_scaff_name = Tkinter.StringVar()
peak_table_name = Tkinter.StringVar()
savelocation = Tkinter.StringVar()

Tkinter.Label(frame, text='Scaffold File').grid(row=0)
Tkinter.Button(frame, text = "Browse", command = loadscaffname).grid(row=0, column=1, sticky ='W', padx=15)
Tkinter.Label(frame, text='Peak Table File').grid(row=1)
Tkinter.Button(frame, text = "Browse", command = loadpeakname).grid(row=1, column=1 , sticky ='W', padx=15)
Tkinter.Button(frame, text = "Save Location", command = choosesavelocation).grid(row=7, column=0, padx=15)
Tkinter.Label(frame, textvariable= df_scaff_name).grid(column = 2, row=0, sticky= 'W', columnspan=2)
Tkinter.Label(frame, textvariable= peak_table_name).grid(column = 2, row=1, sticky= 'W', columnspan=2)
Tkinter.Label(frame, textvariable= savelocation).grid(column = 1, row=7, sticky= 'W', columnspan=2)

stdrow_tk = Tkinter.StringVar()
stdcurveplate_tk = Tkinter.StringVar()
stdpeakname_tk = Tkinter.StringVar()
filename_tk = Tkinter.StringVar()
Tkinter.Label(frame, text='Standard Curve Info:').grid(row=2, columnspan=2, pady=6, sticky='W')
displaygraph = Tkinter.BooleanVar()
check = Tkinter.Checkbutton(tk, text='Display Curve Fit', variable=displaygraph, onvalue=True, offvalue=False).grid(row=7, column = 0, sticky='W')

Tkinter.Label(frame, text='Starting row').grid(row=3, sticky='W')
stdrow_entry = Tkinter.Entry(frame, width=7, textvariable=stdrow_tk).grid(row=3, column=1, sticky='W')
Tkinter.Label(frame, text='Plate').grid(row=4, sticky='W')
stdcurveplate_entry = Tkinter.Entry(frame, width=7, textvariable=stdcurveplate_tk).grid(row=4, column=1, sticky='W')
Tkinter.Label(frame, text='Peak name').grid(row=5, sticky='W')
stdpeakname_entry = Tkinter.Entry(frame, width=7, textvariable=stdpeakname_tk).grid(row=5, column=1, sticky='W', columnspan=2)

Tkinter.Label(frame, text='File name').grid(row=6, sticky='W')
filename_entry = Tkinter.Entry(frame, width=20, textvariable=filename_tk).grid(row=6, column=1, pady = 12, sticky='W', columnspan=2)


Tkinter.Button(frame, text="Calculate", command= close_window).grid(column=1, row=8,  pady = 12, sticky='W')

tk.mainloop()
os.chdir(savelocation.get())

fname_raw = peak_table_name.get()
df_raw = pd.read_csv(fname_raw, comment='#')
fname_scaff = df_scaff_name.get()
scaff = pd.read_csv(fname_scaff, comment='#')
df_scaff = scaff.dropna()

stdcurverow = stdrow_tk.get()
stdcurveplate = int(stdcurveplate_tk.get())
stdpeakname = stdpeakname_tk.get()
stddata, notstddata = pullstdcurvedata(df_raw, stdcurverow, stdcurveplate, stdpeakname)
equationparams = stdcurvequadfit(stddata, displaygraph.get())
sampledata = pulldata(notstddata, df_scaff)
fulldf = calculateconc(sampledata, equationparams,df_scaff)
trimmeddf = pd.concat([fulldf['Well Label'], fulldf['Sample Name'], fulldf['Variant ID'], fulldf['Scaffold'], fulldf['Type'].str.rstrip('*'), fulldf['[IgG], ug/mL'], fulldf['nM'], fulldf['% Purity']], axis=1)

addwells = findblankwells(trimmeddf, df_scaff)

trimmeddf = notdetectedaddition(addwells, df_scaff, trimmeddf)
desiredscaffind = trimmeddf['Scaffold'] == trimmeddf['Type']
exportdf_trim = trimmeddf[desiredscaffind]
addwells2 = findblankwells(exportdf_trim, df_scaff)
exportdf= notdetectedaddition(addwells2, df_scaff, exportdf_trim)

del exportdf['Type']


# Getting dumped: C:\Users\ltong
filename = filename_tk.get() + '.xlsx'
writer = pd.ExcelWriter(filename, engine='xlsxwriter')
exportdf.to_excel(writer, sheet_name='Sheet1')
trimmeddf.to_excel(writer, sheet_name='Sheet2')
workbook = writer.book
worksheet1 = writer.sheets['Sheet1']
worksheet2 = writer.sheets['Sheet2']
worksheet1.conditional_format('G2:G98', {'type': '3_color_scale', 'min_color': '#008000', 'max_color' : '#FF0000'})
worksheet2.conditional_format('H2:H98', {'type': '3_color_scale', 'min_color': '#008000', 'max_color' : '#FF0000'})
writer.save()


