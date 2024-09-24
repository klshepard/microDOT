# -*- coding: utf-8 -*-
"""
@author: Kevin Renehan
"""

import numpy as np
import matplotlib.pyplot as plt

    
#%% Load block average data

npz = np.load('block_average_data.npz')
time = npz['time']
E_hbo = npz['E_hbo']
E_hbr = npz['E_hbr']
M_hbo = npz['M_hbo']
M_hbr = npz['M_hbr']
L_hbo = npz['L_hbo']
L_hbr = npz['L_hbr']

# Task duration
task_duration = 6


#%% Plot settings

plt.style.use('default')
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['lines.markersize'] = 6
legendsize=14


#%% Block average plots

# Detectors to plot
dtp = np.arange(16)

# Create a block average plot for each detector
for d in dtp:
    
    # Figure
    plt.figure()
    
    # Plot 
    plt.plot(time, E_hbo[d], linestyle='-', label='E', color='tab:blue')
    plt.plot(time, M_hbo[d], linestyle='-', label='M', color='tab:orange')
    plt.plot(time, L_hbo[d], linestyle='-', label='L', color='tab:green')
    plt.plot(time, E_hbr[d], linestyle=':', color='tab:blue')
    plt.plot(time, M_hbr[d], linestyle=':', color='tab:orange')
    plt.plot(time, L_hbr[d], linestyle=':', color='tab:green')
    
    plt.axvspan(0, task_duration, facecolor='grey', alpha=0.25)
    plt.xlabel("Time Since Task Onset (s)")
    plt.ylabel(r'$\Delta [Hb(O_2|R)]\ (\mu M)$')
    plt.xticks(ticks=np.arange(-2, 14+2, 2))
    plt.grid(b=True)
    plt.xlim((-2, 14))
    plt.legend(fontsize=legendsize)
    
    