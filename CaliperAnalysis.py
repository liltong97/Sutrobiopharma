from __future__ import division, absolute_import, \
                                    print_function, unicode_literals

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def pullstdcurvedata(stdcurverow, stdcurveplate):
    ''' Finds the herceptin curve data and returns it. Assumes that there
    is a duplicate on the next row and it takes up 12 wells/row'''
    rowApos = df['Sample Name'] == 'A1'
    startrowApos = df[rowApos]['Type']=='?'
    

    
    

    


fname = '..\ltong\Documents\Python test runs\CaliperPeakTable.csv'
df = pd.read_csv(fname, comment='#')

stdcurverow = 'A'
stdcurveplate = 2

# stddata = pullsstdcurvedata(stdcurverow, stdcurveplate)

    
    
    