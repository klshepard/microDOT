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


#%% Determine how many reconstructions were included
field_idx = 0
field_dict = {}
for f in mat_file:
    if "roi" in f:
        
        # Interrogate the field name
        roi_number = int(f.split("_")[0].split("roi")[1])
        recon_type = f.split("_")[1]
        
        # Store
        field_dict[field_idx] = [f, roi_number, recon_type]
        
        # Increment
        field_idx += 1
        
# Figure out how many parameters and locations there are
number_of_locations = 0
number_of_parameters = 0
parameter_list = []
for idx in field_dict:
    
    # Update number of locations
    number_of_locations = max(number_of_locations, field_dict[idx][1])
    
    # Update number of parameters
    if field_dict[idx][2] not in parameter_list:
        parameter_list.append(field_dict[idx][2])
        number_of_parameters += 1
        
# Index order
parameter_idx = {}
for i in range(number_of_parameters):
    parameter_idx[parameter_list[i]] = i
        
# Create array
recon_data = np.hstack((np.array([number_of_parameters, number_of_locations+1]), np.array(np.shape(mat_file[field_dict[0][0]]))))
recon_data = np.zeros(tuple(recon_data))

# Fill array
for p in parameter_idx:
    
    # Get index to fill
    pi = parameter_idx[p]
    
    for l in range(number_of_locations+1):
        
        # Figure out the field that contains data for this parameter and location
        for f in field_dict:
            if (p in field_dict[f]) and (l in field_dict[f]):
                recon_data[pi][l] = mat_file[field_dict[f][0]]
                
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
    plt.figure()
    plt.gca().set_aspect('equal', adjustable='box')

    # Accumulate data and plot
    data_to_plot = np.zeros_like([parameter_idx[p]][0])
    for l in location_list_to_plot:
        data_to_plot = data_to_plot + recon_data[parameter_idx[p]][l]
    plt.pcolormesh(xmesh, ymesh, np.flip(np.array(data_to_plot > 0).astype(int), axis=0), shading='auto', cmap=cmap)
    
    for j, l in enumerate(location_list_to_plot):
        
        # Data to plot
        data_to_plot = recon_data[parameter_idx[p]][l]
        
        # Real roi location
        roi_origin = np.array(location_to_position[l])
        
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

    # Finish plot
    plt.title(title, fontsize=16)
    plt.xlabel("X Position (cm)")
    plt.ylabel("Y Position (cm)")
    plt.xlim([-2, 2])
    plt.ylim([-2, 2])
    plt.xticks([-2, -1, 0, 1, 2])
    plt.yticks([-2, -1, 0, 1, 2])
    plt.grid(color='black', linestyle=':')
    
# Find mean localization error and spatial resolution, and convert from cm to mm
le = np.mean(le, axis=1) * 10
sr = np.mean(sr, axis=1) * 10

# Print results
for i, p in enumerate(parameter_list):
    print(p + " Mean Localization Error: " + str(round(le[i], 2)) + " mm")
    print(p + " Mean Spatial Resolution: " + str(round(sr[i], 2)) + " mm")
