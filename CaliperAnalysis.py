from __future__ import division, absolute_import, \
                                    print_function, unicode_literals

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def pullstdcurvedata(stdcurverow, stdcurveplate, peakname):
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

def stdcurvequadfit(df_stdcurve, dispgraph):
    ''' Fits a quadratic curve to the herceptin standard data''' 
    fitdata = df_stdcurve['Corr. Area']
    x= [10, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000] #doublecheck
    # do polyfit...
    
    if dispgraph:
        plot(x, fitdata)
        
        # plot fitted data ontop...
        
def calculateconc(df, curvefitparam):
    
    


fname = '..\ltong\Documents\Python test runs\CaliperPeakTable.csv'
df = pd.read_csv(fname, comment='#')

stdcurverow = 'A'
stdcurveplate = 2
stdpeakname = 'IgG'
stddata = pullstdcurvedata(stdcurverow, stdcurveplate, stdpeakname)

    
    
    