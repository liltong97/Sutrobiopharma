from __future__ import division, absolute_import, \
                                    print_function, unicode_literals

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


fname = '..\ltong\Documents\Python test runs\CaliperPeakTable.csv'
df = pd.read_csv(fname, comment='#')

stdcurverow = 'A'
stdcurveplate = 2

def pullstdcurvedata(stdcurverow, stdcurveplate):
    