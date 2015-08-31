from __future__ import division, absolute_import, \
                                    print_function, unicode_literals

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

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
    herceptinstdwells = df.iloc[ilocstartindex:ilocendindex]  
    
    # pull out desired peaks 
    identifiedpeaks = herceptinstdwells['Type'] == peakname
    forcedpeaks = herceptinstdwells['Type'] == peakname + '*'
    allpeaks = identifiedpeaks | forcedpeaks
    df_stdcurve = herceptinstdwells[allpeaks]
    
    return df_stdcurve
    
def findmissingwells(df):
    """ Will return an array of strings with which wells are missing in the 
    dataframe passed using the Sample Name column."""
    wells = df['Sample Name']
    missingwells = []
    startingwell = wells.iloc[0]
    startingrow = startingwell[0]
    startingcolumn = startingwell[1:3]
    if startingcolumn != '1':
        for i in np.linspace(1, int(startingcolumn), int(startingcolumn)+1):
            missingwells += [startingrow + str(i)]
        
    lastwell = wells.iloc[-1]
    theoreticalnumwells = int(lastwell[1:3]) - int(startingcolumn) + 12* (ord(lastwell[0]) - ord(startingrow))
    print(theoreticalnumwells)
    
    for i in np.linspace(0, theoreticalnumwells-1, theoreticalnumwells):
        currwellnum = (i % 12)+1
        currwellrow = chr(int(ord(startingrow) + i/12))
        currwell = currwellrow + str(int(currwellnum))
        if not currwell in wells.values:
           missingwells += [currwell] 
    return missingwells
      
      
def calcavgstd(df_stdcurve):
    """ Calculate the average of the two herceptin standard curve. Return an array
    with the average values. """
    avgcorrareas = []
    stdcorrareas = df_stdcurve['Corr. Area']
    if len(stdcorrareas) != 24:
        if len(stdcorrareas) > 24:
            print("You have extra peaks in your herceptin curve, wtf")
        else:
            misswells = findmissingwells(df_stdcurve)
        
    else:
        for i in np.linspace(0, 11, 12):
            avgcorrareas += [(stdcorrareas.iloc[i] + stdcorrareas.iloc[i+12]) / 2]
    

    return avgcorrareas        

def stdcurvequadfit(df_stdcurve, dispgraph):
    ''' Fits a quadratic curve to the average of the herceptin standard data''' 
    fitdata = df_stdcurve['Corr. Area']
    x= [10, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000] #doublecheck
    avgstdcurve = calcavgstd(df_stdcurve)
    # do polyfit...
    np.polyfit(x, avgstdcurve , 2)
    
    if dispgraph:
        plot(x, fitdata)
        
        # plot fitted data ontop...
        
#def calculateconc(df, curvefitparam):
    
    


fname = '..\ltong\Documents\Python test runs\CaliperPeakTable.csv'
df = pd.read_csv(fname, comment='#')

stdcurverow = 'A'
stdcurveplate = 2
stdpeakname = 'IgG'
stddata = pullstdcurvedata(df, stdcurverow, stdcurveplate, stdpeakname)

    
    
    