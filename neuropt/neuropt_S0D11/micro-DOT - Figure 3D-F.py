# -*- coding: utf-8 -*-
"""
@author: Kevin Renehan
"""

#%% Imports
import numpy as np
import matplotlib.pyplot as plt


#%% Load from npy file

npz = np.load('contrast_data.npz')
depth = npz['depth']
nir_contrast = npz['nir_contrast']
nir_cnr = npz['nir_cnr']
nir_bkg_captures = npz['nir_bkg_captures']



#%% Plot settings

plt.style.use('default')
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['lines.markersize'] = 6
legendsize=14

#%% Figure 3E and 3F


# Source and detector to plot
s_to_plot = [0]
d_to_plot = [11]

# Set full_windows
full_windows = True
contrast_defined = True
cnr_defined = True


# Choose windows to plot (unused if full_windows==False)
w_to_plot = [0,6,20]
legend = ['Early Window', 'Late Window', 'All Windows']

# Which things to plot
plot_contrast = True
plot_cnr = True

# Override windows for early, total, and late
w_to_plot = [0, 1, 2] if not full_windows else w_to_plot

if contrast_defined and plot_contrast:
    
    # Create a plot
    for s in s_to_plot:
        for d in d_to_plot:
            plt.figure()
            title_header = "Contrast"
            for i in w_to_plot:
                plt.plot(depth, nir_contrast[i][s][d], marker='o', markeredgecolor='black', alpha=0.8)
            plt.xlabel("Depth of ROI (mm)")
            plt.xlim((4, 32))
            plt.xticks(ticks=depth)
            plt.ylabel(title_header)
            plt.legend(legend, fontsize=legendsize)
            plt.grid()


if cnr_defined and plot_cnr:
       
    # Create a plot
    for s in s_to_plot:
        for d in d_to_plot:
            plt.figure()
            title_header = "Contrast"
            for i in w_to_plot:
                plt.plot(depth, nir_cnr[i][s][d], marker='o', markeredgecolor='black', alpha=0.8)
            plt.xlabel("Depth of ROI (mm)")
            plt.xlim((4, 32))
            plt.xticks(ticks=depth)
            plt.ylabel(title_header + " to Noise Ratio")
            plt.legend(legend, fontsize=legendsize)
            plt.grid()
    
    
#%% Distance formula

import math

# Find distance between source and detector on 4x4
def distance(source, detector):
    
    # Get rows
    row_source = math.floor(source/4)
    row_detector = math.floor(detector/4)
    
    # Get columns
    col_source = source % 4
    col_detector = detector % 4
    
    # Calculate the distance
    vertical_distance = (row_detector-row_source)*8
    horizontal_distance = (col_detector-col_source)*8
    distance = (horizontal_distance**2 + vertical_distance**2)**0.5
    return round(distance,2)
    
        
#%% Figure 3D

# Sources and detectors to plot
s_to_plot = [0]
d_to_plot = [1, 2, 3, 5, 6, 7,  11, 15]
w = 20
idcs = np.where(depth == 16)[0][0]
hist_t = np.arange(170)*65e-3

import matplotlib.ticker as mticker
def log_tick_formatter(val, pos=None):
    return f"$10^{{{int(val)}}}$"

# Plot
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_xlabel('Arrival Time (ns)')
ax.set_ylabel('SDS (mm)')
ax.set_zlabel('Counts')
ax.zaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
ax.zaxis.set_major_locator(mticker.MaxNLocator(integer=True))
ax.set_zlim((0,5))
ax.invert_yaxis()
ax.set_yticks(np.arange(-5, 40, 5))
ax.set_xticks(np.arange(-5, 40, 5))
ax.view_init(elev=25, azim=60)
pad = np.zeros((10,)) + .1

for s in s_to_plot:
    for d in d_to_plot:
        sds_val = distance(s, d)
        sds_axis = np.zeros_like(hist_t) + sds_val
        data = np.log10(nir_bkg_captures[s][d][idcs])
        data = np.hstack([pad, data])
        ax.plot3D(hist_t, sds_axis, data, alpha=0.8, linewidth=2)
