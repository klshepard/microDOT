# -*- coding: utf-8 -*-
"""
@author: Kevin Renehan
"""

import numpy as np
import matplotlib.pyplot as plt

#%% Load block averages data 

npz = np.load('block_averages_data.npz')

time_ba = npz['time_ba']
cw_hbo_ba = npz['cw_hbo_ba']
mt_hbo_ba = npz['mt_hbo_ba']
l_hbo_ba = npz['l_hbo_ba']
cw_hbr_ba = npz['cw_hbr_ba']
mt_hbr_ba = npz['mt_hbr_ba']
l_hbr_ba = npz['l_hbr_ba']

        
#%% Make some plots

plt.style.use('default')
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['lines.markersize'] = 6
legendsize=14


s = 0
dtp = [2,6,9,14]
idcs = np.where(time_ba <= 80)
flist = []

plt.close("all")
for d in dtp:
    f = plt.figure()
    flist.append(f)
    
    # Plot
    plt.plot(time_ba[idcs], cw_hbo_ba[s][d][idcs]*1e6, linestyle='-', label='Intensity', color='tab:blue')
    plt.plot(time_ba[idcs], mt_hbo_ba[s][d][idcs]*1e12, linestyle='-', label='Mean Time', color='tab:orange')
    plt.plot(time_ba[idcs], l_hbo_ba[s][d][idcs]*1e12, linestyle='-', label='Laplace', color='tab:green')
    plt.plot(time_ba[idcs], cw_hbr_ba[s][d][idcs]*1e6, linestyle=':', color='tab:blue')
    plt.plot(time_ba[idcs], mt_hbr_ba[s][d][idcs]*1e12, linestyle=':', color='tab:orange')
    plt.plot(time_ba[idcs], l_hbr_ba[s][d][idcs]*1e12, linestyle=':', color='tab:green')
    
    plt.axvspan(0, 60, facecolor='grey', alpha=0.25)
    plt.xlabel("Time Since Onset (s)")
    plt.ylabel(r'$\Delta [HbO_2]\ or\ \Delta [HbR]\ (\mu M)$')
    plt.xticks(ticks=np.arange(0, 80, 10))
    plt.grid(b=True)
    plt.xlim((-2, 80))
    # plt.title("SDS " + str(distance(s,d)) + "mm")
    plt.legend(loc='upper left',fontsize=legendsize)