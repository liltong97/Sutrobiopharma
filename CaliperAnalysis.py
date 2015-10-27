from __future__ import division, absolute_import, \
                                    print_function, unicode_literals

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as stats
import xlsxwriter
import Tkinter

def pullstdcurvedata(df, stdcurverow, stdcurveplate, peakname):
    ''' Finds the herceptin curve data and returns it. Assumes that there
    is a duplicate on the next row and it takes up 12 wells/row'''
    
    startwell = stdcurverow.upper() + '1'
    # Find start of herceptin curve
    rowpos = df['Sample Name'] == startwell
    startrowpos = df[rowpos]['Type']=='?' # every well starts with a ? peak
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
    print(startingcolumn)
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
    conversions = {'IgG': 150000, 'GFP': 28000, 'scFvFc': 100000, 'scfv': 23000}
    combineddf_full = pd.merge(sampledata, df_scaff, on=['Well Label'])
    #combineddf_full = combineddf_full.replace({"Scaffold":conversions})
    conc =  []
    nMconc = []
    for i in np.linspace(0, len(combineddf_full)-1, len(combineddf_full)):
        scaffnorm = combineddf_full['Scaffold'].iloc[i]
        nMconc += [(concentrations.iloc[i] / conversions[scaffnorm]) * 10**6]
        conc += [concentrations.iloc[i]]
    combineddf_full['[IgG], ug/mL'] = conc
    combineddf_full['nM'] = nMconc

    return combineddf_full

top = Tkinter.Tk()
top.mainloop()

fname_raw = '..\ltong\Documents\Caliper Data\Raw Data\pHexpsamples_2015-10-22_03-27-53_PeakTable.csv'
df_raw = pd.read_csv(fname_raw, comment='#')
fname_scaff = '..\ltong\Documents\Caliper Data\Raw Data\scaffold_1988.csv'
df_scaff = pd.read_csv(fname_scaff, comment='#')

stdcurverow = 'G'
stdcurveplate = 2
stdpeakname = 'IgG'
stddata, notstddata = pullstdcurvedata(df_raw, stdcurverow, stdcurveplate, stdpeakname)
equationparams = stdcurvequadfit(stddata, True)
sampledata = pulldata(notstddata, df_scaff)
fulldf = calculateconc(sampledata, equationparams,df_scaff)
trimmeddf = pd.concat([fulldf['Well Label'], fulldf['Sample Name'], fulldf['Variant ID'], fulldf['Scaffold'], fulldf['Type'].str.rstrip('*'), fulldf['[IgG], ug/mL'], fulldf['nM']], axis=1)

addwells = findmissingwells(trimmeddf, True)

for i in range(len(addwells)):
    sample = addwells[i]
    if len(sample) == 2:
        welllabel= sample[0] + '0' + sample[1]
    else:
        welllabel = sample
    a=df_scaff[df_scaff['Well Label'] == welllabel]
    varID = a['Variant ID'].iloc[0]
    scaff = a['Scaffold'].iloc[0]
    trimmeddf.loc[i+len(trimmeddf)] = [welllabel, sample, varID, scaff, scaff, 'Not Detected', 0]
desiredscaffind = trimmeddf['Scaffold'] == trimmeddf['Type']
exportdf_trim = trimmeddf[desiredscaffind]
addwells2 = findmissingwells(exportdf_trim, True)
for i in range(len(addwells2)):
    sample = addwells2[i]
    if len(sample) == 2:
        welllabel= sample[0] + '0' + sample[1]
    else:
        welllabel = sample
    a=df_scaff[df_scaff['Well Label'] == welllabel]
    varID = a['Variant ID'].iloc[0]
    scaff = a['Scaffold'].iloc[0]
    exportdf_trim.loc[i+len(exportdf_trim)] = [welllabel, sample, varID, scaff, scaff, 'Not Detected', 0]
exportdf = exportdf_trim.sort(['Well Label'])
del exportdf['Type']
exportdf_full = trimmeddf.sort(['Well Label'])
    
# Getting dumped: C:\Users\ltong
writer = pd.ExcelWriter('simple.xlsx', engine='xlsxwriter')
exportdf.to_excel(writer, sheet_name='Sheet1')
exportdf_full.to_excel(writer, sheet_name='Sheet2')
writer.save()


