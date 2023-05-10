# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 15:00:57 2023

@author: atesmer
"""

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
             
#%%============================================================================
#                       Re-Generate dataset
# =============================================================================
#generate_h5py(datafile)

#%%============================================================================
# 
# =============================================================================

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
            
            n_trials = len(expgrp['Action']['Bandit'])/expgrp['Action'].attrs['Duration']*60


            df.loc[i] = [mouse,sex,cohort,state,ei+1,n_trials]
            i += 1
            
f, ax = plt.subplots(1,1,figsize=(4,4),dpi=400)
sns.lineplot(df,x='Day',y='Trials',hue='Sex',units='Mouse',estimator=None,palette=['teal','crimson'],alpha=0.75)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylabel('Trials Per Minute')
ax.set_xlabel('Training Day')
ax.set_xticks(df['Day'].unique())
f.savefig(r'P:\Alexander\Bandit\Figures\Training_Trials.svg')
# ax.set_ylim([0,5])
plt.show()
#%%============================================================================
#                       Check Slopes
# =============================================================================
cohort = 'WT'
mode = 'Prob'
state = 'FastedVeh'
day = 1

df = pd.DataFrame(columns=['Mouse','Sex','Cohort','State','Day','Time','Trials'])
i = 0
with h5py.File(datafile,'r') as file:
    for mouse in file[cohort].keys():
        if mode not in file[cohort][mouse].keys():
            continue
        
        
        # Determine sex of mouse
        sex = file[cohort][mouse].attrs['Sex']
        
        # Sort experiment days
        exp_days = list(file[cohort][mouse][mode][state].keys())
        exp_days.sort()
        
        for ei, exp in enumerate(exp_days):
            expgrp = file[cohort][mouse][mode][state][exp] # Shortcut to exp
            
            time = expgrp['Action']['Time'][()]
            trials = np.cumsum(np.ones(len(time)))

            df.loc[i] = [mouse,sex,cohort,state,ei,time,trials]
            i += 1
            
f, ax = plt.subplots(1,1,figsize=(2,2),dpi=400)


for mouse in df[df['Day']==day]['Mouse'].unique():
    dfm = df[(df['Mouse']==mouse)&(df['Day']==day)]
    x = dfm['Time'].to_numpy().item()
    y = dfm['Trials'].to_numpy().item()
    sex = dfm['Sex'].item()

    ax.plot(x,y,color='teal' if sex == 'm' else 'orchid',lw=1,alpha=1)
    
ax.set_ylim(bottom=0)
ax.set_xlim([0,30*60])
ax.set_xticks(np.arange(0,30*60+1,60*10))
ax.set_xticklabels(np.arange(0,30*60+1,60*10)//60)
ax.set_xlabel('Time (Minutes)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylabel('Cummulative Trials')
ax.set_title('[{}] {}: Day {}'.format(mode,state,day+1))
f.savefig(r'P:\Alexander\Bandit\Figures\Slope.svg')
plt.show()

#%%============================================================================
#                           Deterministic Biases
# =============================================================================

cohort = 'WT'

state = 'Fasted'

df = pd.DataFrame(columns=['Mouse','Sex','Cohort','State','Day','Pump','Prob'])
i = 0
with h5py.File(datafile,'r') as file:
    for mouse in file[cohort].keys():
        if mode not in file[cohort][mouse].keys():
            continue
        
        # Determine sex of mouse
        sex = file[cohort][mouse].attrs['Sex']
        
        # Sort experiment days
        exp_days = list(file[cohort][mouse]['Det'][state].keys())
        exp_days.sort()
        if len(exp_days)!=4:
            continue
        
        for ei, exp in enumerate(exp_days[-1:]):
            expgrp = file[cohort][mouse]['Det'][state][exp] # Shortcut to exp
                        
            choice = expgrp['Action']['Choice'][()].astype(str)
            
            pumps, counts = np.unique(choice,return_counts=True)
            counts = counts.astype(float)/sum(counts)
            
            for pump, count in zip(pumps,counts):
                df.loc[i] = [mouse,sex,cohort,state,ei,pump,count]
                i += 1

f, ax = plt.subplots(1,1,figsize=(0.75,1.5),dpi=400)
hue = None
g=sns.barplot(ax=ax,data=df,x='Pump',y='Prob',color='grey',hue=hue,palette=['grey','grey'],errorbar=('se',1),errwidth=1,capsize=0.1)
sns.lineplot(data=df,x='Pump',y='Prob',color='black',units='Mouse',estimator=None,lw=0.25,alpha=0.25)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.axhline(1/3,color='crimson',lw=0.5,ls='--',zorder=0)
ax.set_ylim([0,0.5])
ax.tick_params(axis='y',labelsize=8)
ax.set_title('Deterministic',fontsize=8)

dfttest = df.pivot(index='Mouse',columns='Pump',values='Prob').reset_index()
for pump in ['A','B','C']:
    t,p = stats.ttest_1samp(dfttest[pump],1/3)
    print('[{}] p={:.4f}'.format(pump,p))
    
f.savefig(r'P:\Alexander\Bandit\Figures\Pump-Bias.svg')
#==============================================================================
#                           Bandit Side-Biases
# =============================================================================

cohort = 'WT'
mode = 'Det'
state = 'Fasted'

df = pd.DataFrame(columns=['Mouse','Sex','Cohort','State','Day','Bandit','Prob'])
i = 0
with h5py.File(datafile,'r') as file:
    for mouse in file[cohort].keys():
        if mode not in file[cohort][mouse].keys():
            continue
        
        # Determine sex of mouse
        sex = file[cohort][mouse].attrs['Sex']
        
        # Sort experiment days
        exp_days = list(file[cohort][mouse][mode][state].keys())
        exp_days.sort()
        
        for ei, exp in enumerate(exp_days[-1:]):
            expgrp = file[cohort][mouse][mode][state][exp] # Shortcut to exp
            
            probs = [expgrp['Action'].attrs['{}_Prob'.format(_)] for _ in ['A','B','C']]
            
            bandit = expgrp['Action']['Bandit'][()].astype(str)
            choice = expgrp['Action']['Choice'][()].astype(str)

            for b in ['AB','BC','AC']:
                pump, count = np.unique(choice[bandit==b],return_counts=True)
                
                if b == 'AB':
                    left = 'A'
                elif b == 'BC':
                    left = 'B'
                elif b == 'AC':
                    left = 'C'
                    
                count = float(count[pump==left])/sum(count)
                
                count -= 0.5
                count = abs(count)

                df.loc[i] = [mouse,sex,cohort,state,ei,b,count]
                i += 1


df = df.groupby(['Mouse','Sex','Cohort','State','Day']).mean(numeric_only=True).reset_index()

hue = 'Sex'

f, ax = plt.subplots(1,1,figsize=(1,2),dpi=400)
# g=sns.boxplot(data=df,x='Day',y='Prob',hue='Sex',color='grey',palette=['teal','orchid'],**PROPS,whis=1.5e30)
g=sns.barplot(ax=ax,data=df,x='State',y='Prob',color='grey',hue=hue,palette=['teal','orchid'],errorbar=('se',1),errwidth=1,capsize=0.1)
g.legend_.remove()
sns.stripplot(data=df,x='State',y='Prob',hue='Sex',color='black',dodge=True,s=1,legend=False)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# ax.axhline(0.5,color='crimson',lw=0.5,ls='--',zorder=0)
ax.set_ylim([0,0.5])
ax.tick_params(axis='both',labelsize=8)
ax.set_title('Deterministic',fontsize=8)
ax.set_ylabel('Left or Right Bias')
ax.set_xlabel('')
plt.show()
f.savefig(r'P:\Alexander\Bandit\Figures\Side-Bias.svg')


#%%============================================================================
#                           Bandit Prob-Biases over Time
# =============================================================================
from scipy import signal
cohort = 'WT'
mode = 'Prob'
state = 'FastedVeh'

df = pd.DataFrame(columns=['Mouse','Sex','Cohort','State','Day','Bandit','T','Prob'])
i = 0
with h5py.File(datafile,'r') as file:
    for mouse in file[cohort].keys():
        if mode not in file[cohort][mouse].keys():
            continue
        # if mouse != '1384':
        #     continue

        # Determine sex of mouse
        sex = file[cohort][mouse].attrs['Sex']
        
        # Sort experiment days
        try:
            exp_days = list(file[cohort][mouse][mode][state].keys())
            exp_days.sort()
        except KeyError:
            continue
        
        for ei, exp in enumerate(exp_days[:]):
            expgrp = file[cohort][mouse][mode][state][exp] # Shortcut to exp
                        
            bandit = expgrp['Action']['Bandit'][()].astype(str)
            reward = expgrp['Action']['Reward'][()].astype(str)
            choice = expgrp['Action']['Choice'][()].astype(str).astype(object)
            time = expgrp['Action']['Time'][()].astype(float)
            for pump in ['A','B','C']: 
                choice[choice==pump] = expgrp['Action'].attrs[pump+'_Prob']
                
            print([expgrp['Action'].attrs[pump+'_Prob'] for pump in ['A','B','C']])

            
            dur = expgrp['Action'].attrs['Duration']
            t_bins = np.arange(0,dur+1,5*60)
            # t_bins = np.array([0,*np.linspace(300,30*60,1+1)])
            # t_bins = [0,300,600,900]
            

            for b in ['AB','BC','AC']:
                pump = np.sort([expgrp['Action'].attrs[_+'_Prob'] for _ in b])
                
                choice_b = choice[bandit==b]
                reward_b = reward[bandit==b]
                time_b = time[bandit==b]
                
                y = choice_b==pump[np.argmax(pump)]
                # y = signal.savgol_filter(y,8,1)
                
                data = stats.binned_statistic(time_b,y,bins=t_bins).statistic
                # data -= data[0]
                data -= 0.5
                if np.isnan(data).any():
                    print(data)
                                    
                choices = '-'.join(list(pump.astype(str)))

                for t, d in zip(t_bins,data):
                    df.loc[i] = [mouse,sex,cohort,state,ei,choices,t//60,d]
                    i += 1

df = df.groupby(['Mouse','Sex','Cohort','State','Bandit','T']).mean(numeric_only=True).reset_index()

f, ax = plt.subplots(1,3,figsize=(5,2.5),dpi=400)

for i, bandit in enumerate(['0.25-1.0','0.5-1.0','0.25-0.5']):
    dfb = df[df['Bandit']==bandit]
    hue = None
    
    # g=sns.boxplot(ax=ax[i],data=dfb,x='T',y='Prob',color='grey',hue=hue,**PROPS,whis=1.5e10,palette=['teal','crimson'])
    g=sns.barplot(ax=ax[i],data=dfb,x='T',y='Prob',color='grey',hue=hue,palette=['teal','orchid'],errorbar=('se',1),errwidth=1,capsize=0.1)

    try:
        g.legend_.remove()
    except:
        pass
    
    for mouse in dfb['Mouse'].unique():
        
        dfm = dfb[dfb['Mouse']==mouse]
        dfm = dfm.pivot(index=['Mouse','Sex'],columns='T',values='Prob').reset_index()
        y = dfm[np.array(t_bins[:-1])//60].to_numpy().flatten()
        x = np.arange(len(y)) + (-0.2 if dfm['Sex'].item()=='m'else 0.2)+0.1
        ax[i].scatter(x,y,color='black',alpha=0.75,marker='o',s=1,lw=0)
        # ax[i].plot(np.arange(len(y)),y,color='black',lw=0.25,alpha=0.5)
    
    
    # sns.lineplot(ax=ax[i],data=dfb,x='T',y='Prob',color='dimgrey',units='Mouse',estimator=None,lw=0.25,alpha=0.5,legend=False,)
    ax[i].axvspan(-0.5,0.5,color='lightgrey',zorder=0)


    ax[i].spines['top'].set_visible(False)
    ax[i].spines['right'].set_visible(False)
    # ax[i].set_ylim([-0.7,0.7])
    ax[i].set_ylim([-0.51,0.51])
    ax[i].axhline(0,lw=0.5,ls='--',color='black',zorder=0)
    ax[i].set_yticks(np.linspace(-0.5,0.5,6))
    ax[i].set_yticklabels(['0.0','0.3','0.4','0.6','0.8','1.0'])

    
    ax[i].set_xticks(np.arange(len(t_bins)-1))
    ax[i].set_xticklabels(['-']+list(t_bins[1:-1]//60))
    gs = (np.array(bandit.split('-')).astype(float)*100).astype(int)
    ax[i].set_title('G$_{{{}}}^{{{}}}$'.format(*gs))
    
    ax[i].tick_params(axis='both',labelsize=8)
    ax[i].set_xlabel('')
    # ax[i].set_ylabel('Δ Prob Higher')
    ax[i].set_ylabel('p(Higher)-p(Bias)')
    ax[i].set_xlim(left=-0.5)
    
# # g.legend_.remove()

f.savefig(r'P:\Alexander\Bandit\Figures\Gamble_Sucess.svg')
#%%============================================================================
#                           Per-Bandit Success Rate
# =============================================================================
cohort = 'WT'
mode = 'Prob'
states = ['FastedVeh','FastedAlm']
# states = ['Fed','Fasted','FastedVeh']
# states = ['Fed','Fasted']


begin, end = [int(5*60),30*60]
# begin, end = [int(5*60),int(17.5*60)]
# begin, end = [0,5*60]


hue= 'State'
# palette=['teal','orchid']
palette=['grey','forestgreen']
# palette=['black','grey']


# =============================================================================
# 
# =============================================================================
if hue is None:
    contrast = 'Cohort'
else:
    contrast = hue
    if hue == 'State':
        hue_order = states
    elif hue == 'Sex':
        hue_order = ['m','f']
    else:
        hue_order = None
        
df = pd.DataFrame(columns=['Mouse','Sex','Cohort','State','Day','Bandit','Prob'])
i = 0
with h5py.File(datafile,'r') as file:
    for mouse in file[cohort].keys():
        for state in states:

            # Determine sex of mouse
            sex = file[cohort][mouse].attrs['Sex']
            # if sex != 'f':
            #     continue
            
            # Sort experiment days
            try:
                exp_days = list(file[cohort][mouse][mode][state].keys())
                exp_days.sort()
            except KeyError:
                continue
            
            for ei, exp in enumerate(exp_days[:]):
                expgrp = file[cohort][mouse][mode][state][exp] # Shortcut to exp
                            
                bandit = expgrp['Action']['Bandit'][()].astype(str)
                reward = expgrp['Action']['Reward'][()].astype(bool)
                choice = expgrp['Action']['Choice'][()].astype(str).astype(object)
                time = expgrp['Action']['Time'][()].astype(float)
                for pump in ['A','B','C']: 
                    choice[choice==pump] = expgrp['Action'].attrs[pump+'_Prob']
                    
                # print([expgrp['Action'].attrs[pump+'_Prob'] for pump in ['A','B','C']])
    
                dur = expgrp['Action'].attrs['Duration']
    
                for b in ['AB','BC','AC']:
                    pump = np.sort([expgrp['Action'].attrs[_+'_Prob'] for _ in b])
                    
                    choice_b = choice[bandit==b]
                    reward_b = reward[bandit==b]
                    time_b = time[bandit==b]
                    
                    y = choice_b==pump[np.argmax(pump)]                    
                   
                    data = y[(time_b>=begin)&(time_b<=end)].mean()
                    # baseline = y[(time_b>=0)&(time_b<5*60)].mean()
                    # data -= (baseline-0.5)
                    # print(baseline)
                                        
                    choices = '-'.join(list(pump.astype(str)))
    
                    df.loc[i] = [mouse,sex,cohort,state,ei,choices,data]
                    i += 1

if len(states) ==3:
    df['State'] = [_ if _ != 'FastedVeh' else 'Fasted' for _ in df['State'] ]
    hue_order = [_ for _ in states if _ != 'FastedVeh']

df = df.groupby(['Mouse',contrast,'Bandit'])['Prob'].mean(numeric_only=True).reset_index()
order = ['0.25-0.5','0.5-1.0','0.25-1.0',]
df['Bandit'] = pd.Categorical(df['Bandit'],categories=order,ordered=True)
# df = df[df['Sex']=='m']

f, ax = plt.subplots(1,1,figsize=(1,2),dpi=400)

# g=sns.barplot(ax=ax,data=df,x='Bandit',y='Prob',color='grey',hue=hue,palette=palette,errorbar=('se',1),errwidth=1,capsize=0.1)
g=sns.lineplot(ax=ax,data=df,x='Bandit',y='Prob',color='grey',hue=hue,palette=palette,errorbar=('se',1),marker='s',err_style='bars',alpha=0.8,ms=4,mew=0,hue_order=hue_order)

try:
    g.legend_.remove()
except:
    pass

# Lines
for b, bandit in enumerate(order):
    

    dfb = df[df['Bandit']==bandit].pivot(index=['Mouse'],columns=hue,values='Prob').reset_index()
    xx = np.array([-0.2,0.2])+b
    
    for m, mouse in enumerate(dfb['Mouse'].unique()):
        yy = dfb[dfb['Mouse']==mouse][hue_order].to_numpy().flatten()
        ax.plot(xx,yy,color='black',lw=0.25,alpha=0.25,zorder=0)
        
    # Do ANOVAs and T-Tests
    ttests = np.vstack(dfb[hue_order].to_numpy())
    if hue == 'Sex':
        if b == 0:
            anova = pg.mixed_anova(data=df,within='Bandit',between=hue,dv='Prob',subject='Mouse')
            print(anova.round(4))
            print('\n')
            posthoc = pg.pairwise_tests(data=df,within='Bandit',between=hue,dv='Prob',subject='Mouse')
            print(posthoc)
        t,p = stats.ttest_ind(ttests[:,0],ttests[:,1],nan_policy='omit')
        print('[{}] t={:.3f} p={:.3f} ind'.format(bandit,t,p))
    else:
        if b == 0:
            anova = pg.rm_anova(data=df,within=['Bandit',hue],dv='Prob',subject='Mouse')
            print(anova.round(4))
            print('\n')
            # for h in hue_order:
            #     anova = pg.rm_anova(data=df[df[hue]==h],within=['Bandit'],dv='Prob',subject='Mouse')
            #     print(anova.round(4))
            #     print('\n')
            #     posthoc = pg.pairwise_tests(data=df[df[hue]==h],within='Bandit',dv='Prob',subject='Mouse')
            #     print(posthoc)
            # posthoc = pg.pairwise_tests(data=df,within='Bandit',between=hue,dv='Prob',subject='Mouse')
            # print(posthoc)
        # t,p = stats.ttest_rel(ttests[:,0],ttests[:,1])
        # print('[{}] t={:.3f} p={:.3f} paired'.format(bandit,t,p))


sns.stripplot(ax=ax,data=df,x='Bandit',y='Prob',palette=palette,hue=hue,jitter=0,dodge=True,size=1,legend=False,order=order,alpha=0.5,hue_order=hue_order)


ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylim([0.2,0.9])
ax.set_ylim([0.3,0.9])
ax.set_yticks(np.arange(0.3,0.9+0.1,0.1))
ax.axhline(0.5,lw=0.5,ls='--',color='black',zorder=0)


ax.set_xticks([0,1,2])

xticklabels = ['G$_{{{}}}^{{{}}}$'.format(*(np.array(bandit.split('-')).astype(float)*100).astype(int)) for bandit in order]
ax.set_xticklabels(xticklabels)

ax.tick_params(axis='both',labelsize=8)
ax.set_xlabel('')
ax.set_ylabel('p(higher)')
ax.set_xlabel('Bandit')

f.savefig(r'P:\Alexander\Bandit\Figures\Gamble_Sucess.svg')

#%%============================================================================
#                Reward Identification Over Time       
# =============================================================================
from scipy.ndimage import uniform_filter1d
cohort = 'WT'
mode = 'Prob'
states = ['Fed','Fasted']
states = ['Fed','Fasted','FastedVeh','FastedAlm']

hue = 'State'
# =============================================================================
# 
# =============================================================================
if hue is None:
    contrast = 'Cohort'
    oneplot = True
else:
    contrast = hue
    if hue == 'State':
        hue_order = states
    elif hue == 'Sex':
        hue_order = ['m','f']
    else:
        hue_order = None
    oneplot = False

df = pd.DataFrame(columns=['Mouse','Sex','Cohort','State','Day','Data25','Data50','Data100'])
i = 0
with h5py.File(datafile,'r') as file:
    for mouse in file[cohort].keys():
        for state in states:
        
            # Determine sex of mouse
            sex = file[cohort][mouse].attrs['Sex']
            
            # Sort experiment days
            try:
                exp_days = list(file[cohort][mouse][mode][state].keys())
                exp_days.sort()
            except KeyError:
                continue
            
            for ei, exp in enumerate(exp_days[:]):
                expgrp = file[cohort][mouse][mode][state][exp] # Shortcut to exp
                
                dur = expgrp['Action'].attrs['Duration']
                
                mb = 5
                t_bins = np.arange(0,dur+1,mb*60)
                
                time = expgrp['Action']['Time'][()]
                bandit = expgrp['Action']['Bandit'][()].astype(str)
                choice = expgrp['Action']['Choice'][()].astype(str).astype(object)
                for pump in ['A','B','C']:
                    choice[choice==pump] = expgrp['Action'].attrs[pump+'_Prob']
                    
                data = []
                for pump in [0.25,0.5,1]:
                    out = stats.binned_statistic(time,choice==pump,bins=t_bins).statistic
                    out = interpolate_nans(out)
                    # out = savgol_filter(out,5,1)
                    # out = uniform_filter1d(out,5,mode='nearest')
                    # out -= (out[0]-1/3)
                    data.append(out)
    
                df.loc[i] = [mouse,sex,cohort,state,ei,*data]
                i += 1
            
df = df.groupby(['Mouse',hue])[['Data25','Data50','Data100']].mean(numeric_only=False).reset_index()


if oneplot:
    f, ax = plt.subplots(1,1,figsize=(1,2),dpi=400)
    for i, dset in enumerate(['Data25','Data50','Data100']):
        color = ['#D81B60','#004D40','#FFC107',][i]
        
            
        y = np.vstack(df[dset].to_numpy())
        ym = np.nanmean(y,axis=0)
        ye = np.nanstd(y,axis=0)/np.sqrt(y.shape[0])
        x = np.arange(len(ym))+0.5
        ax.plot(x,ym,color=color,lw=1,alpha=0.75,marker='s',ms=2,mew=0)
        ax.fill_between(x,ym+ye,ym-ye,color=color,lw=0,alpha=0.25)
            
            
        for mouse in df['Mouse'].unique():
            dfm = df[df['Mouse']==mouse]            
            y = np.vstack(dfm[dset].to_numpy()).flatten()
            x = np.arange(len(y))+0.5
            ax.plot(x,y,color=color,lw=0.25,alpha=0.1)

    ax.axhline(1/3,color='black',lw=1,ls='--',zorder=1)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim([1/6,3/6])
    # ax.set_title('')
    ax.set_xticks(np.arange(0,len(t_bins),1))
    xticklabels = (t_bins//60).astype(object)
    xticklabels[1::2] = ''
    # xticklabels = np.array(['']*len(t_bins)).astype(object)
    # xticklabels[::5] = (t_bins//60).astype(object)[::5]
    ax.set_xticklabels(xticklabels)
    ax.set_xlabel('Time Bin (min)')
    ax.axvspan(0,5/mb,color='lightgrey',zorder=0)
    ax.set_xlim(0,len(t_bins)-1)
    ax.tick_params(axis='both',labelsize=8)
    ax.set_ylabel('p(chosen)')
    f.savefig(r'P:\Alexander\Bandit\Figures\Probability_Preference.svg')
    plt.show()

else:
    f, axes = plt.subplots(1,len(hue_order),figsize=(len(hue_order),1.5),dpi=400)
    for ax, hueval in zip(axes,hue_order):
        dfhue = df[df[hue]==hueval]
    
        for i, dset in enumerate(['Data25','Data50','Data100']):
            color = ['#D81B60','#004D40','#FFC107',][i]
            
                
            y = np.vstack(dfhue[dset].to_numpy())
            ym = np.mean(y,axis=0)
            ye = y.std(0)/np.sqrt(y.shape[0])
            x = np.arange(len(ym))+0.5
            ax.plot(x,ym,color=color,lw=2,alpha=0.75,marker='s',ms=4,mew=0)
            ax.fill_between(x,ym+ye,ym-ye,color=color,lw=0,alpha=0.25)
                
                
            for mouse in dfhue['Mouse'].unique():
                dfm = dfhue[dfhue['Mouse']==mouse]            
                y = np.vstack(dfm[dset].to_numpy()).flatten()
                x = np.arange(len(y))+0.5
                ax.plot(x,y,color=color,lw=0.5,alpha=0.1)
    
        ax.axhline(1/3,color='black',lw=1,ls='--',zorder=1)
    
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylim([1/6,3/6])
        # ax.set_title('')
        ax.set_xticks(np.arange(0,len(t_bins),1))
        xticklabels = (t_bins//60).astype(object)
        xticklabels[1::2] = ''
        # xticklabels = np.array(['']*len(t_bins)).astype(object)
        # xticklabels[::5] = (t_bins//60).astype(object)[::5]
        ax.set_xticklabels(xticklabels)
        ax.set_xlabel('Time Bin (min)')
        ax.axvspan(0,5/mb,color='lightgrey',zorder=0)
        ax.set_xlim(0,len(t_bins)-1)
        ax.tick_params(axis='both',labelsize=8)
        ax.set_ylabel('p(chosen)')
        ax.set_title(hueval)
    f.savefig(r'P:\Alexander\Bandit\Figures\Probability_Preference.svg')
    plt.show()

#%%============================================================================
#                           Circle-Back Markov
# =============================================================================
cohort = 'WT'
mode = 'Prob'
states = ['FastedVeh','FastedAlm']
# states = ['Fed','Fasted']
# states = ['FastedAlm']
# states = ['FastedVeh','Fasted']

begin, end = [int(5*60),30*60]

hue= 'State'
# palette=['teal','orchid']
palette=['dimgrey','lightgrey','lightgrey','forestgreen']
# palette=['grey','black']


# =============================================================================
# 
# =============================================================================
if hue is None:
    contrast = 'Cohort'
else:
    contrast = hue
    if hue == 'State':
        hue_order = states
    elif hue == 'Sex':
        hue_order = ['m','f']
    else:
        hue_order = None
        
df = pd.DataFrame(columns=['Mouse','Sex','Cohort','State','Day','PumpProb','CircleBackProb'])
i = 0
with h5py.File(datafile,'r') as file:
    for mouse in file[cohort].keys():
        for state in states:

            # Determine sex of mouse
            sex = file[cohort][mouse].attrs['Sex']
            # if sex != 'f':
            #     continue
            
            # Sort experiment days
            try:
                exp_days = list(file[cohort][mouse][mode][state].keys())
                exp_days.sort()
            except KeyError:
                continue
            
            for ei, exp in enumerate(exp_days[:]):
                expgrp = file[cohort][mouse][mode][state][exp] # Shortcut to exp
                            
                reward = expgrp['Action']['Reward'][()].astype(bool)
                choice = expgrp['Action']['Choice'][()].astype(str).astype(object)
                time = expgrp['Action']['Time'][()].astype(float)
                for pump in ['A','B','C']: 
                    choice[choice==pump] = expgrp['Action'].attrs[pump+'_Prob']
                    
                reward = reward[(time>=begin)&(time<=end)]
                choice = choice[(time>=begin)&(time<=end)]
                dur = expgrp['Action'].attrs['Duration']
                
                # circle_back = {0.25:[],0.5:[],1.0:[]}
                # for ci,c in enumerate(choice[:-2]):
                #     # if not reward[ci+1]:
                #         circle_back[c].append(choice[ci+2] == c)
                # for _ in [0.25,0.5,1.0]:
                #     circle_back_mean = np.mean(circle_back[_])
                #     df.loc[i] = [mouse,sex,cohort,state,ei,_,circle_back_mean]
                #     i += 1
                
                circle_back = []
                for ci,c in enumerate(choice[:-2]):
                    # if not reward[ci]:
                        circle_back.append(choice[ci+2] == c)
                        
                circle_back_mean = np.mean(circle_back)
                df.loc[i] = [mouse,sex,cohort,state,ei,'All',circle_back_mean]
                i += 1

# if len(states) ==3:
#     df['State'] = [_ if _ != 'FastedVeh' else 'Fasted' for _ in df['State']]
#     hue_order = ['Fed','Fasted']

df = df.groupby(['Mouse',contrast,'PumpProb'])['CircleBackProb'].mean(numeric_only=True).reset_index()

f, ax = plt.subplots(1,1,figsize=(0.5,2),dpi=400)

g = sns.barplot(data=df,x=hue,y='CircleBackProb',errorbar=('se',1),capsize=0.1,errwidth=1,order=hue_order,palette=palette)
sns.stripplot(data=df,x=hue,y='CircleBackProb',s=1,order=hue_order,palette='dark:black',jitter=0)

# g.legend_.remove()

dfb = df.pivot(index=['Mouse'],columns=hue,values='CircleBackProb').reset_index()
xx = np.arange(len(hue_order))

for m, mouse in enumerate(dfb['Mouse'].unique()):
    yy = dfb[dfb['Mouse']==mouse][hue_order].to_numpy().flatten()
    ax.plot(xx,yy,color='black',lw=0.25,alpha=0.25,zorder=3)

anova = pg.rm_anova(data=df,dv='CircleBackProb',subject='Mouse',within='State')
print(anova)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylim([0.0,1.0])
ax.set_yticks(np.arange(0,1.1,0.2))
ax.axhline(0.5,lw=0.5,ls='--',color='black',zorder=0)

ax.tick_params(axis='both',labelsize=8)
ax.set_xlabel('')
ax.set_ylabel('p(return)')
ax.set_xlabel('Bandit')
plt.show()

#%%============================================================================
#                               Performance Analysis
# =============================================================================
cohort = 'WT'
mode = 'Prob'
states = ['FastedVeh','FastedAlm']
# states = ['Fed','Fasted']

df = pd.DataFrame(columns=['Mouse','Sex','Cohort','State','Day','R','Re','NTrials'])

i = 0
with h5py.File(datafile,'r') as file:
    for mouse in file[cohort].keys():
        for state in states:
        
            # Determine sex of mouse
            sex = file[cohort][mouse].attrs['Sex']
            
            # Sort experiment days
            try:
                exp_days = list(file[cohort][mouse][mode][state].keys())
                exp_days.sort()
            except KeyError:
                continue
        
            
            for ei, exp in enumerate(exp_days[:]):
                expgrp = file[cohort][mouse][mode][state][exp] # Shortcut to exp
                
                probs = {}
                for _ in ['A','B','C']:
                    probs[_]= expgrp['Action'].attrs[_+'_Prob']
                reward = expgrp['Action']['Reward'][()].astype(bool)
                choice = expgrp['Action']['Choice'][()].astype(str)
                time = expgrp['Action']['Time'][()].astype(float)
                bandit = expgrp['Action']['Bandit'][()].astype(str)

                begin, end = [5,30]
                begin, end = [17.5,30]
                begin, end = [5,17.5]
                begin, end = [5,30]
                reward = reward[(time>=begin*60)&(time<=end*60)]
                choice = choice[(time>=begin*60)&(time<=end*60)]
                bandit = bandit[(time>=begin*60)&(time<=end*60)]
                                
                reward_e = []
                # Calcualte Average Reward if Mouse Chooses Randomly
                for b in bandit:
                    reward_e.append(np.mean([probs[_] for _ in b]))

                assert len(reward_e) == len(reward), 'mismatch'

                re = np.mean(reward_e)
                r = np.mean(reward)
                r -= re
                    
                df.loc[i] = [mouse,sex,cohort,state,ei+1,r,re,len(reward)]
                i += 1
            
df = df.groupby(['Mouse','Sex','Cohort','State']).mean(numeric_only=False).reset_index()
df['Day'] = [1 for _ in df['Day']]

x = 'State'
order=states
hue = None
palette = ['teal','orchid']
palette = ['dimgrey','lightgrey']

f, ax = plt.subplots(1,2,figsize=(1.5,2),dpi=400)
g = sns.barplot(df,ax=ax[0],x=x,y='R',hue=hue,errorbar=('se',1),errwidth=1,capsize=0.1,palette=palette,order=order)
sns.stripplot(df,ax=ax[0],x=x,y='R',hue=hue,palette=['black','black'],jitter=0,dodge=False,size=1,order=order,)
try:
    g.legend_.remove()
except AttributeError:
    pass
ax[0].set_ylabel('p(reward)-p(chance)')
ax[0].set_ylim([-0.06,0.16])
ax[0].axhline(0,color='black',ls='--',lw=0.5,zorder=0)
g = sns.barplot(df,ax=ax[1],x=x,y='NTrials',hue=hue,errorbar=('se',1),errwidth=1,capsize=0.1,palette=palette,order=order)
sns.stripplot(df,ax=ax[1],x=x,y='NTrials',hue=hue,palette=['black','black'],jitter=0,dodge=False,size=1,order=order)

if hue != 'Sex':
    for i, val in enumerate(['R','NTrials']):
        dfx = df.pivot(index='Mouse',columns=x,values=val).reset_index()
        for mouse in dfx['Mouse'].unique():
            dfm = dfx[dfx['Mouse']==mouse]            
            yy = np.vstack(dfm[order].to_numpy()).flatten()
            xx = np.arange(len(yy))
            ax[i].plot(xx,yy,color='black',lw=0.25,alpha=0.25)


try:
    g.legend_.remove()
except AttributeError:
    pass
ax[1].set_ylabel('Valid Trials')
ax[1].set_ylim([0,200])
sns.despine()
for a in ax.flatten():
    a.tick_params(axis='both',labelsize=8)
    a.set_xlabel('')
    a.tick_params(axis='x',rotation=90)
# f.tight_layout()
f.savefig(r'P:\Alexander\Bandit\Figures\Performance_Metrics.svg')
plt.show()