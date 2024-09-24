#%% Handle imports

import numpy as np
import matplotlib.pyplot as plt
import os
from distance import distance

#%% Load previously calculated data

# Initialize top level dictionary
dlist = []
dist = []
sd_pair = []

# Depth range
depth_range = (6, 30)

# Laplace coefficient
Ls = -1e-4

# Source 0, detector 8
source, detector = 0, 8
d = distance(source, detector)
sd_pair.append((source, detector))
data = np.load('neuropt_S0D8/neuropt_S0D8' + '_Ls' + "{:.0e}".format(Ls) + '.npz')
idcs = np.where((data['depth'] >= depth_range[0]) & (data['depth'] <= depth_range[1]))
dict_temp = {'depth': data['depth'][idcs]}
keys = list(data.keys())
keys.remove('depth')
keys.remove('laplace_coefficient')
for k in keys:
    dict_temp[k] = np.take(data[k], indices=idcs[0], axis=-1)
dist.append(d)
dlist.append(dict_temp)

# Source 0, detector 9
source, detector = 0, 9
d = distance(source, detector)
sd_pair.append((source, detector))
data = np.load('neuropt_S0D9/neuropt_S0D9' + '_Ls' + "{:.0e}".format(Ls) + '.npz')
idcs = np.where((data['depth'] >= depth_range[0]) & (data['depth'] <= depth_range[1]))
dict_temp = {'depth': data['depth'][idcs]}
keys = list(data.keys())
keys.remove('depth')
keys.remove('laplace_coefficient')
for k in keys:
    dict_temp[k] = np.take(data[k], indices=idcs[0], axis=-1)
dist.append(d)
dlist.append(dict_temp)

# Source 0, detector 10
source, detector = 0, 10
d = distance(source, detector)
sd_pair.append((source, detector))
data = np.load('neuropt_S0D10/neuropt_S0D10' + '_Ls' + "{:.0e}".format(Ls) + '.npz')
idcs = np.where((data['depth'] >= depth_range[0]) & (data['depth'] <= depth_range[1]))
dict_temp = {'depth': data['depth'][idcs]}
keys = list(data.keys())
keys.remove('depth')
keys.remove('laplace_coefficient')
for k in keys:
    dict_temp[k] = np.take(data[k], indices=idcs[0], axis=-1)
dist.append(d)
dlist.append(dict_temp)

# Source 0, detector 7
source, detector = 0, 7
d = distance(source, detector)
sd_pair.append((source, detector))
data = np.load('neuropt_S0D7/neuropt_S0D7' + '_Ls' + "{:.0e}".format(Ls) + '.npz')
idcs = np.where((data['depth'] >= depth_range[0]) & (data['depth'] <= depth_range[1]))
dict_temp = {'depth': data['depth'][idcs]}
keys = list(data.keys())
keys.remove('depth')
keys.remove('laplace_coefficient')
for k in keys:
    dict_temp[k] = np.take(data[k], indices=idcs[0], axis=-1)
dist.append(d)
dlist.append(dict_temp)

# Source 0, detector 15
source, detector = 0, 15
d = distance(source, detector)
sd_pair.append((source, detector))
data = np.load('neuropt_S0D15/neuropt_S0D15' + '_Ls' + "{:.0e}".format(Ls) + '.npz')
idcs = np.where((data['depth'] >= depth_range[0]) & (data['depth'] <= depth_range[1]))
dict_temp = {'depth': data['depth'][idcs]}
keys = list(data.keys())
keys.remove('depth')
keys.remove('laplace_coefficient')
for k in keys:
    dict_temp[k] = np.take(data[k], indices=idcs[0], axis=-1)
dist.append(d)
dlist.append(dict_temp)

# Depth axis should be the same
depth = dlist[0]['depth']


#%% Plot settings

plt.style.use('default')
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['lines.markersize'] = 6
legendsize=14



#%% Figure 3G - Histogram of SDS separations with single color
    
# Source detector separations
sds = []
for s in range(16):
    for d in range(16):
        sds.append(distance(s, d))
sds = np.array(sds)

# Make histogram
plt.figure()
scalp_total = 0
brain_total = 0
sds_axis = np.arange(-1, 35+2, 2)
N, bins, patches = plt.hist(sds, bins=len(sds_axis), range=(min(sds_axis), max(sds_axis)+2), edgecolor='Black')
for i, p in enumerate(patches):
    p.set_facecolor("tan")
    scalp_total += N[i]
plt.xticks(sds_axis[0:-1]+1)
plt.grid(axis='y')
plt.ylim([0, 60])
plt.xlim([-2, 36])
plt.ylabel("Number of Channels")
plt.xlabel("Source Detector Separation (mm)")


#%% Figure 3H - Plot showing relative sensitivity of mean time

plot_pairs = \
    ( \
        ('mean_time', 1), \
        ('mean_time', 2), \
        ('mean_time', 3), \
        ('contrast', 4), \
     )

# Load data
laplace = -1*np.array([dlist[k]['laplace'][sd_pair[k][0]][sd_pair[k][1]] for k in range(len(dlist))])
variance = np.array([dlist[k]['variance'][sd_pair[k][0]][sd_pair[k][1]] for k in range(len(dlist))])
mean_time = np.array([dlist[k]['mean_time'][sd_pair[k][0]][sd_pair[k][1]] for k in range(len(dlist))])
contrast = np.array([dlist[k]['contrast'][20][sd_pair[k][0]][sd_pair[k][1]] for k in range(len(dlist))])

# Ratio of deeper/shallower data structure
laplace_ratio = np.zeros_like(laplace)
variance_ratio = np.zeros_like(variance)
mean_time_ratio = np.zeros_like(mean_time)
contrast_ratio = np.zeros_like(contrast)

for k in range(len(dlist)):
    contrast_ratio[k] = contrast[k] / np.max(contrast[k])
    mean_time_ratio[k] = mean_time[k] / np.max(mean_time[k])
    variance_ratio[k] = variance[k] / np.max(variance[k])
    laplace_ratio[k] = laplace[k] / np.max(laplace[k])

plt.figure()
for k, d in plot_pairs:
    
    # Color
    if d == 4:
        color = '#1f77b4'
    elif d == 0:
        color = '#ff7f0e'
    elif d == 1:
        color = '#2ca02c'
    elif d == 2:
        color = '#d62728'
    elif d == 3:
        color = '#9467bd'
    
    if k == 'contrast':
        plt.plot(depth, contrast_ratio[d], label= 'E ' + str(round(dist[d], 1)) + ' mm', marker='*', markersize=10, markeredgecolor='black')
    elif k == 'mean_time':
        plt.plot(depth, mean_time_ratio[d], label= 'M ' + str(round(dist[d], 1)) + ' mm', marker='o', markeredgecolor='black')
    elif k == 'laplace':
        plt.plot(depth, laplace_ratio[d], label= 'L ' + str(round(dist[d], 1)) + ' mm', marker='o', markeredgecolor='black')
plt.legend(title='SDS')
plt.xlabel("Depth of ROI (mm)")
plt.ylim([-0.05, 1.05])
plt.xticks(ticks=depth)
plt.xlim([5, 27])
plt.ylabel("Relative Sensitivity")
plt.grid()
plt.legend()


#%% Figure 3I - Plot showing relative sensitivity of mean time

plot_pairs = \
    ( \
        ('laplace', 1), \
        ('laplace', 2), \
        ('laplace', 3), \
        ('contrast', 4), \
     )

# Load data
laplace = -1*np.array([dlist[k]['laplace'][sd_pair[k][0]][sd_pair[k][1]] for k in range(len(dlist))])
variance = np.array([dlist[k]['variance'][sd_pair[k][0]][sd_pair[k][1]] for k in range(len(dlist))])
mean_time = np.array([dlist[k]['mean_time'][sd_pair[k][0]][sd_pair[k][1]] for k in range(len(dlist))])
contrast = np.array([dlist[k]['contrast'][20][sd_pair[k][0]][sd_pair[k][1]] for k in range(len(dlist))])

# Ratio of deeper/shallower data structure
laplace_ratio = np.zeros_like(laplace)
variance_ratio = np.zeros_like(variance)
mean_time_ratio = np.zeros_like(mean_time)
contrast_ratio = np.zeros_like(contrast)

for k in range(len(dlist)):
    contrast_ratio[k] = contrast[k] / np.max(contrast[k])
    mean_time_ratio[k] = mean_time[k] / np.max(mean_time[k])
    variance_ratio[k] = variance[k] / np.max(variance[k])
    laplace_ratio[k] = laplace[k] / np.max(laplace[k])

plt.figure()
for k, d in plot_pairs:
    
    # Color
    if d == 4:
        color = '#1f77b4'
    elif d == 0:
        color = '#ff7f0e'
    elif d == 1:
        color = '#2ca02c'
    elif d == 2:
        color = '#d62728'
    elif d == 3:
        color = '#9467bd'
    
    if k == 'contrast':
        plt.plot(depth, contrast_ratio[d], label= 'E ' + str(round(dist[d], 1)) + ' mm', marker='*', markersize=10, markeredgecolor='black')
    elif k == 'mean_time':
        plt.plot(depth, mean_time_ratio[d], label= 'M ' + str(round(dist[d], 1)) + ' mm', marker='o', markeredgecolor='black')
    elif k == 'laplace':
        plt.plot(depth, laplace_ratio[d], label= 'L ' + str(round(dist[d], 1)) + ' mm', marker='o', markeredgecolor='black')
plt.legend(title='SDS')
plt.xlabel("Depth of ROI (mm)")
plt.ylim([-0.05, 1.05])
plt.xticks(ticks=depth)
plt.xlim([5, 27])
plt.ylabel("Relative Sensitivity")
plt.grid()
plt.legend()