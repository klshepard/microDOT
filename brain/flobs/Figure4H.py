# -*- coding: utf-8 -*-
"""
@author: Kevin Renehan
"""

import numpy as np
from scipy.io import loadmat
from scipy.interpolate import BSpline
import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize as op


#%% PLOT SETTINGS

plt.style.use('default')
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['lines.markersize'] = 6
legendsize=14

#%% LOAD FILE

# Time domain
mat = loadmat('micro_dot_recon.mat')
td_hbo = mat['slices_hbo']
td_hbr = mat['slices_hb']


#%% MANUALLY ADDED PARAMETERS

# Frame rate
clock_frequency_hz = 50e6
measurements_per_pattern = 178500
patterns_per_frame = 28

# Slice occurence (value of 2 implies that there was 1 slice for every 2 frames in the dataset)
slice_occurence = 2

# How long before the onset does the block average being
time_before_onset_s = 6

# Task duration
task_duration = 6

# Calculate FPS
fps = 1 / (1/clock_frequency_hz * measurements_per_pattern * patterns_per_frame)
tstep = 1/fps*slice_occurence

# Create time axis
time_axis_s = np.arange(len(td_hbo))*tstep - time_before_onset_s


# Basis times and onsets
basis_t = np.linspace(0, 281.8, 1410, endpoint=True)
basis_onsets = np.array([ 21.52,  47.52,  73.52,  99.52, 125.52, 151.52, 177.52, 203.52, 229.52, 255.52])


#%% Load flobs information


# Load flobs file
flobs_data = np.array(pd.read_table("hrfbasisfns.txt", sep='  ', dtype=float)) # 3 basis functions
flobs_data = np.transpose(flobs_data)

# Create a flobs time axis
flobs_tres = 0.05
flobs_time_s = np.arange(np.shape(flobs_data)[1]) * flobs_tres

# Plot basis functions
plt.figure()
colors = ['tab:blue', 'tab:green', 'tab:cyan', 'tab:pink', 'tab:purple']
for i in range(len(flobs_data)):
    plt.plot(flobs_time_s, flobs_data[i], label='Basis Function ' + str(i+1), color=colors[i])
plt.xticks(np.arange(0,28))
plt.grid(b=True)
plt.legend()


#%% FUNCTIONS


def intensity_in_region(row_col_coords, data):
    '''Averages the intensity of data in a square of size square_size centered at (center_row, center_column)'''
    
    # Move the frame axis to last
    d = np.transpose(data, axes=(1,2,0))
    
    # Get the unique rc coordinates
    rc = np.array(list(set(row_col_coords)))
    
    # Take from row
    intensity = np.empty( (len(rc), np.shape(d)[-1]) )
    for i, (r, c) in enumerate(rc):
        intensity[i] = d[r][c]
    
    # Mean
    d_mean = np.mean(intensity, axis=0)
    d_std = np.std(intensity, axis=0)
    
    # Return 
    return d_mean, d_std


def make_mask(row_col_coords, data, fill_val=-1):
    
    # Move the frame axis to last
    d = data.copy()
    d = np.transpose(d, axes=(1,2,0))
    
    # Take from row
    for i, (r, c) in enumerate(row_col_coords):
        d[r][c] = np.zeros(len(data)) + fill_val
    
    # Transpose back
    d = np.transpose(d, axes=(2,0,1))
    
    # Return 
    return d


def row_col_coords_from_meshgrid(meshgrid_output):
    
    return [(meshgrid_output[0].flat[i], meshgrid_output[1].flat[i]) for i in range(len(meshgrid_output[0].flat))]
    


    
#%% Load Petros' HRF

# Pull from files
petros_hrf = loadmat('fmri_hrf.mat')['pp']['hrfVis']
while petros_hrf.dtype == 'O':
    petros_hrf = petros_hrf[0]
petros_t = loadmat('fmri_hrf.mat')['pp']['time']
while petros_t.dtype == 'O':
    petros_t = petros_t[0]
petros_hrf = np.squeeze(petros_hrf)
petros_t = np.squeeze(petros_t)


# Pad and truncate
tres = (petros_t[1]-petros_t[0])
pad_length = int(abs(np.min(time_axis_s) / tres))

# Create a stim boxcar
petros_t = np.hstack([np.linspace(petros_t[0]-pad_length*tres, petros_t[0], num=pad_length), petros_t])
petros_hrf = np.hstack([np.zeros(pad_length), petros_hrf])

# Create boxcar
boxcar = np.zeros_like(petros_hrf)
boxcar[np.where((petros_t >= 0) & (petros_t <= 6))] = 1

# Convolution bounds
t_left = 2*petros_t[0]

# Convolve
conv_hrf = np.convolve(petros_hrf, boxcar, mode='full')
conv_t = np.arange(len(conv_hrf))*tres+t_left+basis_onsets[0]

# Spline
spl = BSpline(conv_t, conv_hrf, k=1)
spline_hrf = spl(basis_t)

# Spline to data
spl = BSpline(petros_t, petros_hrf, k=0)
petros_hrf_spline_to_data = spl(flobs_time_s)


#%% Region definition

# Diagonal band consistent with fMRI
region = \
    row_col_coords_from_meshgrid(np.meshgrid(np.arange(15, 20), np.arange(25,30))) + \
    row_col_coords_from_meshgrid(np.meshgrid(np.arange(23, 28), np.arange(32,37))) + \
    row_col_coords_from_meshgrid(np.meshgrid(np.arange(18, 26), np.arange(30,33))) + \
    row_col_coords_from_meshgrid(np.meshgrid(np.arange(20, 24), np.arange(28,33)))
    
    
#%% PLOT CHANGES OVER TIME IN REGION

# Choose data to plot
data_list = [td_hbo] 
title_list = ['TD HbO'] 

# Time limits
time_limits = np.array((-2, 18))

for data_to_plot, title in zip(data_list, title_list):
    
    # Create and fill mask
    mask = make_mask(region, data_to_plot)
    
    # Transpose
    time_trace_hbo, time_trace_hbo_std = intensity_in_region(region, td_hbo)
    time_trace_hbo = time_trace_hbo - np.mean(time_trace_hbo[np.where( (time_axis_s >= time_limits[0]) & (time_axis_s <= 0) )])
    time_trace_hb, time_trace_hb_std = intensity_in_region(region, td_hbr)
    time_trace_hb = time_trace_hb - np.mean(time_trace_hb[np.where( (time_axis_s >= time_limits[0]) & (time_axis_s <= 0) )])

    # Indices to define time ranges
    idcs = np.where((time_axis_s > time_limits[0]) & (time_axis_s < time_limits[1]))
    basis_idcs = np.where(( (basis_t-basis_onsets[0]) > time_limits[0]) & ((basis_t-basis_onsets[0]) < time_limits[1]))

    # Plot
    f1 = plt.figure()
    
    # Normal plot with fill between HbO
    f1.gca().plot(time_axis_s[idcs], time_trace_hbo[idcs]*1e6, color='tab:red', label='$HbO_2$', marker=None)
    f1.gca().fill_between(time_axis_s[idcs], (time_trace_hbo - time_trace_hbo_std)[idcs]*1e6, (time_trace_hbo + time_trace_hbo_std)[idcs]*1e6, color='tab:red', alpha=0.2)
    
    # Normal plot with fill between HbR
    f1.gca().plot(time_axis_s[idcs], time_trace_hb[idcs]*1e6, color='tab:blue', linestyle='-', label='$HbR$', marker=None)
    f1.gca().fill_between(time_axis_s[idcs], (time_trace_hb - time_trace_hb_std)[idcs]*1e6, (time_trace_hb + time_trace_hb_std)[idcs]*1e6, color='tab:blue', alpha=0.2)
    
    # Y label for HbO2/HbR Plot
    f1.gca().set_ylabel("$Hb(O_2|R)\ Change\ (\mu M)$")
    
    # Axis labels
    f1.gca().set_xlabel("Time Since Task Onset (s)")
    f1.gca().axvspan(0, task_duration, color='grey', alpha=0.25)#, label='Task Block')
    f1.gca().legend(fontsize=legendsize)
    f1.gca().grid(b=True)
    f1.gca().set_xticks(np.arange(-2, 20, 2))
    
#%% Prepare basis functions through convolution with boxcar

# Create a stim boxcar
boxcar = np.zeros_like(flobs_time_s)
boxcar[np.where((flobs_time_s >= 0) & (flobs_time_s <= 6))] = 1

# Convolve
flobs_data_conv = np.array([np.convolve(flobs_data[i], boxcar, mode='full') for i in range(len(flobs_data))])
flobs_time_s_conv = np.arange(len(flobs_data_conv[0]))*tres

# Spline# Spline
flobs_spline = np.zeros((len(flobs_data_conv), len(time_axis_s)))
for i in range(len(flobs_data_conv)):
    spl = BSpline(flobs_time_s_conv, flobs_data_conv[i], k=0)
    flobs_spline[i] = spl(time_axis_s)
    
# Plot indices
idcs = np.where( (flobs_time_s >= 0) & (flobs_time_s < 30) )
    
# Make a plot
fig, ax = plt.subplots(nrows=len(flobs_data), ncols=1, sharex=True)
ax = ax.flat
for i in range(len(flobs_data_conv)):
    ax[i].set_title("Basis Function " + str(i+1))
    ax[i].plot(flobs_time_s[idcs], flobs_data[i][idcs]/max(flobs_data[i][idcs]), label='Original ' + '$(b_{} )$'.format(i+1), color=colors[i], linestyle=':')
    ax[i].plot(flobs_time_s[idcs], boxcar[idcs], label='Stimulus Boxcar', color='k')
    ax[i].plot(flobs_time_s_conv[idcs], flobs_data_conv[i][idcs]/max(flobs_data_conv[i][idcs]), label='After Convolution ' + '$(b\'_{} )$'.format(i+1), color=colors[i])
    ax[i].legend(loc='upper right')
    ax[i].set_ylabel("A.U.")
    ax[i].grid(b=True)
ax[-1].set_xlabel("Time (s)")
    
    
#%% Fit measured data with basis functions

# Use this switch to control if hbo or hbr is plotted
do_hbo = True

if do_hbo:
    time_trace_for_fit = time_trace_hbo
else:
    time_trace_for_fit = -time_trace_hb

class BasisClass:
    
    basis_fxns = None
    idcs = None
    
    def __init__(self):
        pass
    
    def set_number_of_basis_functions(self):
        self.basis_fit = self.basis3_fit

    def basis3_fit(self, time, a, b, c):
        if len(self.basis_fxns) != 3:
            raise Exception("Number of basis functions ({}) does not match number of coefficients ({}) in basis_fit function".format(len(self.basis_fxns), 3))
        return a * self.basis_fxns[0][self.idcs] + b * self.basis_fxns[1][self.idcs] + c * self.basis_fxns[2][self.idcs]
    
# Indices for fit
idcs = np.where( (time_axis_s >= -2) & (time_axis_s <= 18) )

# Provide some information to the basis class
bc_inst = BasisClass()
bc_inst.basis_fxns = flobs_spline
bc_inst.basis_fxns = flobs_spline
bc_inst.set_number_of_basis_functions()
bc_inst.idcs = idcs

popt, pcov = op.curve_fit(bc_inst.basis_fit, time_axis_s[idcs], time_trace_for_fit[idcs])
perr = np.sqrt(np.diag(pcov))

# Plot the parameters
plt.figure()

# Plot the time trace
if do_hbo:
    plt.plot(time_axis_s[idcs], time_trace_for_fit[idcs]*1e6, color='tab:red', label='$\Delta [HbO_2 ]$', linestyle='--')
else:
    plt.plot(time_axis_s[idcs], time_trace_for_fit[idcs]*1e6, color='tab:blue', label='$- \Delta [HbR ]$', linestyle='--')

for i in range(len(flobs_spline)):
    
    # Things for this plot
    fit_data = (flobs_spline[i] * popt[i])[idcs]
    fit_std = perr[i]
    low_fit_data = (flobs_spline[i] * (popt[i] - fit_std))[idcs]
    high_fit_data = (flobs_spline[i] * (popt[i] + fit_std))[idcs]
    color = colors[i]
    
    # Make the plot
    plt.plot(time_axis_s[idcs], fit_data*1e6, color=color, label="$c_{}*b'_{}$".format(i, i))
    plt.fill_between(time_axis_s[idcs], low_fit_data*1e6, high_fit_data*1e6, color=color, alpha=0.2)
    
if do_hbo:
    plt.ylabel('$\Delta [HbO_2 ] \ (\mu M)$')
else:
    plt.ylabel('$- \Delta [HbR ] \ (\mu M)$')
plt.xlabel("Time Since Task Onset (s)")
plt.legend()
plt.grid(b=True)

# Plot coefficients and standard deviations
plt.figure()
plt.grid(b=True)
plt.bar(np.arange(len(flobs_spline))+1, popt[0:len(flobs_spline)]*1e6, yerr=perr[0:len(flobs_spline)]*1e6, color=colors)
plt.xlabel("Basis Function Number")
plt.xticks(np.arange(len(flobs_spline))+1)
plt.ylabel("Fit Coefficient")

    
# Plot the reconstructed HRF
plt.figure()
flobs_fit = np.array([flobs_spline[i]*popt[i] for i in range(len(flobs_spline))])
flobs_fit = np.sum(flobs_fit, axis=0)
if do_hbo:
    plt.plot(time_axis_s[idcs], time_trace_for_fit[idcs]*1e6, color='tab:red', linestyle='-', label='Measurement', alpha=0.8)
else:
    plt.plot(time_axis_s[idcs], time_trace_for_fit[idcs]*1e6, color='tab:blue', linestyle='-', label='Measurement', alpha=0.8)
plt.plot(time_axis_s[idcs], flobs_fit[idcs]*1e6, label='Fit', color='tab:green', alpha=0.8)
plt.legend(fontsize=legendsize)
plt.grid(b=True)
plt.xticks(np.arange(-2, 20, 2))
plt.xlabel("Time Since Task Onset (s)")
if do_hbo:
    plt.ylabel('$\Delta [HbO_2 ] \ (\mu M)$')
else:
    plt.ylabel('$- \Delta [HbR ] \ (\mu M)$')
plt.gca().axvspan(0, task_duration, color='grey', alpha=0.25)



# Calculate the R2 coefficient
ss_res = np.sum( (time_trace_for_fit[idcs] - flobs_fit[idcs]) ** 2)
ss_tot = np.sum( (time_trace_for_fit[idcs] - np.mean(time_trace_for_fit[idcs])) ** 2)
r2 = 1 - (ss_res / ss_tot)
print("R2 value is " + str(round(r2, 2)))

# Data HRF
data_hrf = np.array([popt[i]*flobs_data[i] for i in range(len(flobs_data))])

# Sum basis functions
data_hrf = np.sum(data_hrf, axis=0)


plt.figure()
if do_hbo:
    plt.plot(flobs_time_s, data_hrf/max(data_hrf), color='tab:red', label='$HbO_2$')
else:
    plt.plot(flobs_time_s, data_hrf/max(data_hrf), color='tab:blue', label='$HbR$')
plt.plot(petros_t, petros_hrf/max(petros_hrf), label='$BOLD$', color='black')
plt.xlabel("Time Since Task Onset (s)")
plt.ylabel("Normalized Response")
plt.grid(b=True)
plt.legend(fontsize=legendsize)
plt.gca().axvspan(0, 0.5, color='grey', alpha=0.25)

# Figure out the delay between the two peaks
peak_delay =flobs_time_s[np.argmax(petros_hrf_spline_to_data)] -  flobs_time_s[np.argmax(data_hrf)]
print("Peak delay is " + str(round(peak_delay, 2)) + " s")

# Calculate the correlation
cor = np.corrcoef(data_hrf, petros_hrf_spline_to_data)
print("Correlation coefficient is " + str(round(cor[0][1], 2)))