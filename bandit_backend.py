# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 15:25:19 2023

@author: atesmer
"""

import h5py
import os

def generate_h5py(newfilename,data_dir=r'P:/Alexander/Bandit/Data'):
    males = ['1381','1382','1383','1384','1385']
    with h5py.File(newfilename,'w') as file:
        for path in os.scandir(data_dir):
            if not os.path.isdir(path):
                continue
            targ = os.path.basename(path)
    
            file.create_group(targ)
            
            for datapath in os.scandir(path):
                if not datapath.path.endswith('.data'):
                    continue
                
                filename = os.path.basename(datapath)
                
                # Find indexes of underscores in the file name
                u = [i for i, ltr in enumerate(filename) if ltr == '_']
                
                # Extract information from the title
                mouse = filename[:u[0]]
                state = filename[u[0]+1:u[1]]
                mode = filename[u[1]+1:u[2]]
                date = filename[u[2]+1:u[3]]
                
                if mouse not in file[targ].keys():
                    file[targ].create_group(mouse)
                    file[targ][mouse].attrs['Sex'] = 'm' if mouse in males else 'f'
                if mode not in file[targ][mouse].keys():  
                    file[targ][mouse].create_group(mode)
                if state not in file[targ][mouse][mode].keys():  
                    file[targ][mouse][mode].create_group(state)
    
                # Create a new experiment-group
                experiment = file[targ][mouse][mode][state].create_group(date)
                
                # Copy data from individual experiment file
                with h5py.File(datapath.path,'r') as expfile:
                    expfile.copy('Action',experiment)
                    expfile.copy('Lick',experiment)
                