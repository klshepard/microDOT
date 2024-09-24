# -*- coding: utf-8 -*-
"""
@author: Kevin Renehan
"""

#%% Load the file 
import numpy as np
import matplotlib.pyplot as plt

# Map locations numbers to ROI positions
location_to_position = {
    1: (-8, 8), \
    2: (0, 8), \
    3: (8, 8), \
    4: (-8, 0), \
    5: (0, 0), \
    6: (8, 0), \
    7: (-8, -8), \
    8: (0, -8), \
    9: (8, -8), \
    }

# Create depth axis based on heights
number_of_locations = 9
locations = np.arange(number_of_locations)+1


#%% Load data

npz = np.load('moment_changes.npz')
contrast = npz['contrast']
mean_time = npz['mean_time']
laplace = npz['laplace']


#%% Figure 4B plots

def cylinder(center_x,center_y,radius,height_z):
    z = np.linspace(0, height_z, 50)
    theta = np.linspace(0, 2*np.pi, 50)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radius*np.cos(theta_grid) + center_x
    y_grid = radius*np.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid

from matplotlib.patches import Circle, Rectangle
import mpl_toolkits.mplot3d.art3d as art3d

no_background = True

plt.close('all')
plt.style.use('default')
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['lines.markersize'] = 6
legendsize=14

xpos_arr = np.array([np.arange(-12, 12+8, 8)[::-1] for i in range(4)]).reshape(16)
ypos_arr = np.array([12, 12, 12, 12, 4, 4, 4, 4, -4, -4, -4, -4, -12, -12, -12, -12])

# Sources and detectors to plot
s_to_plot = [0, 3]
d_to_plot = np.arange(16)
loc_to_plot = [1,9]
m_to_plot = [0,1,3]
w_to_plot = [0]

loc_to_plot = np.array([np.where(locations==l)[0][0] for l in loc_to_plot])

color_list = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

dx = 0.8
dy = 0.8
offset_x = dx + 0.05
offset_y = dy + 0.05
source_size = 1

# Plot
for l in loc_to_plot:
    for w in w_to_plot:
        for s in s_to_plot:
    
            # Create figure - 
            fig = plt.figure()
            ax = plt.axes(projection='3d')
            ax.set_xlabel('X Position (cm)')
            ax.set_ylabel('Y Position (cm)')
            ax.set_zlabel('Sensitivity Change')
            xticks=np.arange(-20, 20+10, 10)
            yticks=np.arange(-20, 20+10, 10)
            ax.set_xticks(xticks)
            ax.set_yticks(yticks)
            ax.set_xticklabels([str(xt/10) for xt in xticks])
            ax.set_yticklabels([str(yt/10) for yt in yticks])
            ax.set_xlim([min(xticks)*1.05, max(xticks)*1.05])
            ax.set_ylim([min(yticks)*1.05, max(yticks)*1.05])
            ax.view_init(elev=30, azim=53)
            
            for d in d_to_plot:
    
                # Norm
                norm0 = -1
                norm1 = -1
                norm3 = -1
                
                # Find the largest thing to plot
                m0p = contrast[l][s][w][d]/norm0 * 100
                m1p = mean_time[l][s][w][d]/norm1 * 100
                m3p = laplace[l][s][w][d]/norm3 * 100
                
                # Defaults
                m0_offset_x = offset_x
                m0_offset_y = -offset_y
                m1_offset_x = 0
                m1_offset_y = 0
                m3_offset_x = -offset_x
                m3_offset_y = offset_y
                
                if d == s:
                    
                    source_circle = Circle((xpos_arr[d], ypos_arr[d]), 2, facecolor='darkred', alpha=1, edgecolor='k')
                    ax.add_patch(source_circle)
                    art3d.pathpatch_2d_to_3d(source_circle, z=0, zdir="z")
                    
                else:
                
                    # Plot positive values
                    linestyle = '-' if (m0p>0) else ':'
                    ax.bar3d(xpos_arr[d]+m0_offset_x-dx, ypos_arr[d]+m0_offset_y-dy, 0, dx, dy, np.abs(m0p), alpha=1, color=color_list[0], edgecolor='k', shade=False, linestyle=linestyle)
                    m0_proxy = Rectangle((0, 0), 1, 1, fc=color_list[0], edgecolor='k')
                    linestyle = '-' if (m1p>0) else ':'
                    ax.bar3d(xpos_arr[d]+m1_offset_x-dx, ypos_arr[d]+m1_offset_y-dy, 0, dx, dy, np.abs(m1p), alpha=1, color=color_list[1], edgecolor='k', shade=False, linestyle=linestyle)
                    m1_proxy = Rectangle((0, 0), 1, 1, fc=color_list[1], edgecolor='k')
                    linestyle = '-' if (m3p>0) else ':'
                    ax.bar3d(xpos_arr[d]+m3_offset_x-dx, ypos_arr[d]+m3_offset_y-dy, 0, dx, dy, np.abs(m3p), alpha=1, color=color_list[2], edgecolor='k', shade=False, linestyle=linestyle)
                    m3_proxy = Rectangle((0, 0), 1, 1, fc=color_list[2], edgecolor='k')
                
            # Plot a cylinder
            ax.set_zlim((0, 5*1.1))
            cylinder_radius = 2.5
            cylinder_height = (ax.get_zlim()[1] - ax.get_zlim()[0]) / (ax.get_xlim()[1] - ax.get_xlim()[0]) * cylinder_radius*2
            cylinder_offset = ax.get_zlim()[1]
            Xc,Yc,Zc = cylinder(location_to_position[locations[l]][0],location_to_position[locations[l]][1],cylinder_radius,cylinder_height)
            ax.plot_surface(Xc, Yc, Zc+cylinder_offset, alpha=0.5, color='k')
            floor = plt.Circle((location_to_position[locations[l]][0], location_to_position[locations[l]][1]), cylinder_radius, color='k', alpha=0.5)
            ceiling = plt.Circle((location_to_position[locations[l]][0], location_to_position[locations[l]][1]), cylinder_radius, color='k', alpha=0.5)
            ax.add_patch(floor)
            ax.add_patch(ceiling)
            art3d.pathpatch_2d_to_3d(floor, z=cylinder_offset, zdir="z")
            art3d.pathpatch_2d_to_3d(ceiling, z=cylinder_offset+cylinder_height, zdir="z")
            
            # Plot a circular projection of the ROI onto the z=0 plane
            roi_projection_circle = Circle((location_to_position[locations[l]][0],location_to_position[locations[l]][1]), cylinder_radius, color='k', alpha=0.5)
            ax.add_patch(roi_projection_circle)
            art3d.pathpatch_2d_to_3d(roi_projection_circle, z=0, zdir="z")
            
            zt = np.arange(0, 5+1, 1)
            ax.set_zticks(zt)
            ax.set_zticklabels([str(tick)+"%" for tick in zt])
            
            # Add array outline
            device = plt.Rectangle((-40/2,-40/2), 40, 40, edgecolor=(0,0,0,1), facecolor=(1,1,1,0), linewidth=2)
            ax.add_patch(device)
            art3d.pathpatch_2d_to_3d(device, z=0, zdir="z")
            
            # Add the cable
            device = plt.Rectangle((-40/2-20, -10), 20, 20, edgecolor=(0,0,0,1), facecolor=(1,1,1,0), linewidth=2)
            ax.add_patch(device)
            art3d.pathpatch_2d_to_3d(device, z=0, zdir="z")
            
            # Add the pinholes
            for i, x, y in zip(range(len(xpos_arr)), xpos_arr, ypos_arr):
                if s == i:
                    continue
                device = plt.Rectangle((x-2, y-2), 4, 4, facecolor=(0,0,0,.15), edgecolor='k', linewidth=2)
                ax.add_patch(device)
                art3d.pathpatch_2d_to_3d(device, z=0, zdir="z")