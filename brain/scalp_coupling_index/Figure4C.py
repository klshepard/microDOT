# -*- coding: utf-8 -*-
"""
@author: Kevin Renehan
"""


import numpy as np
import matplotlib.pyplot as plt

#%% Load data

npz = np.load('scalp_coupling_index_data.npz')
sci = npz['sci']
position = npz['position']

number_of_rows = 4
number_of_cols = 4
number_of_detectors = 16


#%% Plot settings

plt.style.use('default')
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['lines.markersize'] = 6
legendsize=14
    

#%% Create plot

# Create plot data structure for reshaping
sci_rshp = np.zeros((number_of_detectors))

# Fill 
for d in range(number_of_detectors):
    
    # Detector number that will fill this index in sci_rshp
    detector_number = np.where(position.flat == d)
    
    # Check that we found the detector number
    if np.shape(detector_number)[1] == 1:
        detector_number = detector_number[0][0]
    else:
        print("Found " + str(np.shape(detector_number)[1]) + " instances of detector " + str(d) + " in position vector")
        # return -1
    
    # Fill index in sci_rshp with the data for the correct detector number
    sci_rshp[d] = sci[detector_number]
    
# Reshape to 4x4
sci_rshp = np.reshape(sci_rshp, newshape=(number_of_rows, number_of_cols))


# Create colormap - the flip here makes it such that the image matches the printed array
plt.figure()
plt.pcolormesh(np.flip(sci_rshp, axis=0), vmin=-1, vmax=1, cmap='RdYlGn', edgecolors='black', linewidth=1)
cb = plt.colorbar()
cb.ax.tick_params(labelsize=12)
   
# Hide x and y axes
plt.gca().axes.xaxis.set_ticks([])
plt.gca().axes.yaxis.set_ticks([])

# Add chip names
for r in range(number_of_rows):
    for c in range(number_of_cols):
        
        # Get the position
        y_pos = (number_of_rows - r) - 0.5
        x_pos = c + 0.5
        
        # Get the chip number
        chip_number = position[r][c]
        
        # Color of text
        color = 'tab:red' if (chip_number==0) else 'black'
        
        # Text
        if (chip_number==0):
            text = 'S' + str(chip_number+1)
        else:
            text = 'D' + str(chip_number+1) + '\n' + str(round(sci_rshp[r][c], 2))
        
        # Add text labels
        plt.text(x_pos, y_pos, text, horizontalalignment='center', verticalalignment='center', \
                 color=color, fontfamily='sans-serif', fontweight='normal', fontproperties='Arial', \
                     fontsize='x-large', fontstretch='extra-expanded')

        plt.grid()