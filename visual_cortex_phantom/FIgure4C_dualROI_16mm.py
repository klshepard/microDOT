# -*- coding: utf-8 -*-
"""
@author: Kevin Renehan
"""

from scipy.io import loadmat
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np


# Load mat file
mat_file = loadmat("single-and-dual-roi-16mm-recon-outputs.mat")

# # Source locations
sx = (1.2, 0.4, -0.4, -1.2, 1.2, 0.4, -0.4, -1.2, 1.2, 0.4, -0.4, -1.2, 1.2, 0.4, -0.4, -1.2)
sy = (1.2, 1.2, 1.2, 1.2, 0.4, 0.4, 0.4, 0.4, -0.4, -0.4, -0.4, -0.4, -1.2, -1.2, -1.2, -1.2)

# X and Y coordinates
xloc = mat_file['xloc'][0]
yloc = mat_file['yloc'][0]

# Make a mesh
xmesh, ymesh = np.meshgrid(xloc, yloc)

# Map locations numbers to ROI positions
location_to_position = {
    1: (-0.8, 0.8), \
    2: (0.0, 0.8), \
    3: (0.8, 0.8), \
    4: (-0.8, 0.0), \
    5: (0.0, 0.0), \
    6: (0.8, 0.0), \
    7: (-0.8, -0.8), \
    8: (0.0, -0.8), \
    9: (0.8, -0.8), \
    }
    
# Relevant data
spatial_E = ['Integrated Intensity', mat_file['spatial_E']]
spatial_M = ['Mean Time', mat_file['spatial_M']]
spatial_L_2params = ['Laplace', mat_file['spatial_L_2params']]


#%% Plot for visualization

plt.close('all')

plot_sources = False
flist = []

# General settings
from matplotlib import colors
w = colors.to_rgba('white', alpha=0)
db = colors.to_rgba("midnightblue", alpha=1)
cmap = colors.ListedColormap([w, db])

# Plot settings
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['font.serif'] = 'Times New Roman'
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['lines.markersize'] = 6
plt.rcParams['savefig.transparent'] = True
legendsize=14

# Plot for visualization
flist.append(plt.figure())
plt.pcolormesh(xmesh, ymesh, np.rot90(np.flip(spatial_E[1], axis=0), k=-1), shading='auto', cmap=cmap)
plt.title(spatial_E[0], fontsize=16)
flist.append(plt.figure())
plt.pcolormesh(xmesh, ymesh, np.rot90(np.flip(spatial_M[1], axis=0), k=-1), shading='auto', cmap=cmap)
plt.title(spatial_M[0], fontsize=16)
flist.append(plt.figure())
plt.pcolormesh(xmesh, ymesh, np.rot90(np.flip(spatial_L_2params[1], axis=0), k=-1), shading='auto', cmap=cmap)
plt.title(spatial_L_2params[0], fontsize=16)

for f in flist:
    if plot_sources:
        f.gca().scatter(sx, sy, label="Source Location")
    f.gca().set_xlabel("X Position (cm)")
    f.gca().set_ylabel("Y Position (cm)")
    f.gca().set_aspect('equal', adjustable='box')
    f.gca().set_xlim([-2, 2])
    f.gca().set_ylim([-2, 2])
    f.gca().set_xticks([-2, -1, 0, 1, 2])
    f.gca().set_yticks([-2, -1, 0, 1, 2])
    f.gca().grid(color='black', linestyle=':')
    f.gca().add_patch(Circle(np.array((0.1, 0.28)), radius = .18, color='red', fill=False, linewidth=2))
    f.gca().add_patch(Circle(np.array((0.1, -0.55)), radius = .18, color='red', fill=False, linewidth=2))

