import os
import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.signal import savgol_filter
import pingouin as pg
plt.rcParams['svg.fonttype'] = 'none'

# TODO Load in data 
# TODO plot the pump bias
# First five minutes are dead
# Left Right bias 

lw = 1
PROPS = {
    'boxprops':{'edgecolor':'black','linewidth': 0},
    'medianprops':{'color':'black','linewidth': lw},
    'whiskerprops':{'color':'black','linewidth': lw},
    'flierprops':{'markeredgewidth': 0,'linewidth': 0,'ms':2.5,'color': 'black'},
    'capprops':{'color':'black','linewidth': lw,}}

datafile = '/home/gummi/banditExperiment/dataset.data'
def nan_helper(y):
    return np.isnan(y), lambda z: z.nonzero()[0]

def interpolate_nans(y):
    nans, x= nan_helper(y)
    y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    return y

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

cohort = 'WT'
mode = 'Det'
state = 'Fasted'

df = pd.DataFrame(columns=['Mouse','Sex','Cohort','State','Day','Trials'])
i = 0
with h5py.File(datafile,'r') as file:
    for mouse in file[cohort].keys():
        
        # Determine sex of mouse
        sex = file[cohort][mouse].attrs['Sex']
        
        # Sort experiment days
        exp_days = list(file[cohort][mouse][mode][state].keys())
        exp_days.sort()
        
        for ei, exp in enumerate(exp_days):
            expgrp = file[cohort][mouse][mode][state][exp] # Shortcut to exp
            print(expgrp["Action"]["Bandit"][()].astype(str))
            