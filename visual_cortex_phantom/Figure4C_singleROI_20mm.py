# -*- coding: utf-8 -*-
"""
@author: Kevin Renehan
"""

from scipy.io import loadmat
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np


# Load mat file
mat_file_previous = loadmat('single-roi-16mm-recon-outputs.mat')
mat_file = loadmat("single-roi-20mm-recon-outputs.mat")

# # Source locations
sx = (1.2, 0.4, -0.4, -1.2, 1.2, 0.4, -0.4, -1.2, 1.2, 0.4, -0.4, -1.2, 1.2, 0.4, -0.4, -1.2)
sy = (1.2, 1.2, 1.2, 1.2, 0.4, 0.4, 0.4, 0.4, -0.4, -0.4, -0.4, -0.4, -1.2, -1.2, -1.2, -1.2)

# X and Y coordinates
xloc = mat_file_previous['xloc'][0]
yloc = mat_file_previous['yloc'][0]

# Make a mesh
xmesh, ymesh = np.meshgrid(xloc, yloc)
mesh_size = len(xmesh)

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
number_of_locations = len(location_to_position)
    
# Create numpy array
parameter_idx = {'E': 0, 'M': 1, 'L': 2}
parameter_list = ['E', 'M', 'L']
number_of_parameters = len(parameter_idx)

# Fill data structure
recon_data = np.zeros((number_of_parameters, number_of_locations, mesh_size, mesh_size))
for k in parameter_idx.keys():
    recon_data[parameter_idx[k]] = mat_file[k]


#%% Combined figure

# Close previous figures
plt.close('all')

from matplotlib import colors
w = colors.to_rgba('white', alpha=0)
db = colors.to_rgba("midnightblue", alpha=1)
cmap = colors.ListedColormap([w, db])

# Plot settings
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['lines.markersize'] = 6
plt.rcParams['savefig.transparent'] = True
legendsize=14

# Location to plot
location_list_to_plot = [1, 3, 5, 7, 9]

# Parameter to plot
parameter_to_plot = None

# Create lists
if location_list_to_plot == None:
    location_list_to_plot = np.arange(1, number_of_locations+1)
    
# Create lists
if parameter_to_plot == None:
    parameter_to_plot = parameter_list
else:
    parameter_to_plot = [parameter_to_plot]
    
# Create holders for averaged localization error and spatial resolution
le = np.zeros((len(parameter_to_plot), len(location_list_to_plot)))
sr = np.zeros((len(parameter_to_plot), len(location_list_to_plot)))
    
for i, p in enumerate(parameter_to_plot):
    
    # Spawn figure
    # plt.figure(figsize=(10,6))
    plt.figure()
    plt.gca().set_aspect('equal', adjustable='box')

    # Accumulate data and plot
    data_to_plot = np.zeros_like([parameter_idx[p]][0])
    for l in location_list_to_plot:
        d_tmp = recon_data[parameter_idx[p]][l-1]
        d_tmp[np.where((xmesh > 1.8) | (xmesh < -1.8))] = 0
        d_tmp[np.where((ymesh > 1.8) | (ymesh < -1.8))] = 0
        d_tmp[np.where(d_tmp < np.max(d_tmp)/4)] = 0
        data_to_plot = data_to_plot + d_tmp
    plt.pcolormesh(xmesh, ymesh, np.flip(np.array(data_to_plot > 0).astype(int), axis=0), shading='auto', cmap=cmap)
    
    for j, l in enumerate(location_list_to_plot):
        
        # Data to plot
        data_to_plot = recon_data[parameter_idx[p]][l-1]
        
        # Real roi location
        roi_origin = np.array(location_to_position[l])
        
        # Mask the outer edges of the reconstructions
        data_to_plot[np.where((xmesh > 1.8) | (xmesh < -1.8))] = 0
        data_to_plot[np.where((ymesh > 1.8) | (ymesh < -1.8))] = 0
        
        # Apply a threshold
        data_to_plot[np.where(data_to_plot < np.max(data_to_plot)/4)] = 0
        
        # Find the centroid of the blob
        blob_coords_x = xmesh[np.where(np.flip(data_to_plot, axis=0)>0)]
        blob_coords_y = ymesh[np.where(np.flip(data_to_plot, axis=0)>0)]
        blob_coords = np.transpose(np.vstack((blob_coords_x, blob_coords_y)))
        centroid = np.array((np.mean(blob_coords_x), np.mean(blob_coords_y)))
        
        # Calculate localization error
        localization_error_2D = roi_origin - centroid
        localization_error = np.linalg.norm(localization_error_2D)
        
        # Calculate distances between points in blob and real_roi_origin
        distances = np.linalg.norm(np.subtract(blob_coords, roi_origin), axis=1)
        spatial_resolution_radius = np.max(distances)
        spatial_resolution_diameter = spatial_resolution_radius * 2
        
        # Store
        le[i][j] = localization_error
        sr[i][j] = spatial_resolution_diameter
        
        # Build plot title
        if p == 'E':
            title = 'Integrated Intensity'
        elif p == 'L':
            title = "Laplace"
        elif p == 'M':
            title = "Mean Time"
        
        # Plot
        alpha = 0.75
        plt.gca().add_patch(Circle(roi_origin, radius = .25, color='red', fill=False, linewidth=3, alpha=alpha))
        plt.gca().add_patch(Circle(roi_origin, radius=spatial_resolution_radius, color='yellow', fill=False, linewidth=3, alpha=alpha))
        plt.scatter(roi_origin[0], roi_origin[1], color='red', marker='o', label="ROI Center", alpha=alpha)
        plt.scatter(centroid[0], centroid[1], color='tomato', marker='o', label="Recon Center", alpha=alpha)
        plt.arrow(roi_origin[0], roi_origin[1], centroid[0] - roi_origin[0], centroid[1] - roi_origin[1], head_width = 0.05, length_includes_head=True, color='red', label="Localization Error", alpha=alpha)
        # plt.text(-2.4, -2.4, "Localization Error: {} mm\nSpatial Resolution: {} mm".format((round(localization_error*10, 1)), round(spatial_resolution_diameter*10, 1)), color='white', horizontalalignment="left", verticalalignment="bottom")
        # plt.savefig("C:\\Users\\Dell-User\\Downloads\\MOANA3_Visual_Cortex_Phantom\\figures\\" + title, bbox_inches='tight', dpi=600)
        
    # Finish plot
    # plt.colorbar()
    # plt.scatter(sx, sy, label="Source Location")
    plt.title(title, fontsize=16)
    plt.xlabel("X Position (cm)")
    plt.ylabel("Y Position (cm)")
    plt.xlim([-2, 2])
    plt.ylim([-2, 2])
    plt.xticks([-2, -1, 0, 1, 2])
    plt.yticks([-2, -1, 0, 1, 2])
    plt.grid(color='black', linestyle=':')
    # plt.text(-2.4, -2.4, "Localization Error: {} mm\nSpatial Resolution: {} mm".format((round(localization_error*10, 1)), round(spatial_resolution_diameter*10, 1)), color='white', horizontalalignment="left", verticalalignment="bottom")
    # plt.savefig("C:\\Users\\Dell-User\\Downloads\\MOANA3_Visual_Cortex_Phantom\\figures\\" + title, bbox_inches='tight', dpi=600)

le_std = np.std(le, axis=1) * 10
le = np.mean(le, axis=1) * 10
sr_std = np.std(sr, axis=1)*10
sr = np.mean(sr, axis=1) * 10

# Print results
for i, p in enumerate(parameter_list):
    print(p + " Mean Localization Error: " + str(round(le[i], 2)) + " mm")
    print(p + " Mean Spatial Resolution: " + str(round(sr[i], 2)) + " mm")