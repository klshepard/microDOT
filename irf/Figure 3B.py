# -*- coding: utf-8 -*-
"""
@author: Kevin Renehan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Load data
nir_irf = np.load('nir_irf_s4_0.8v.npz')['irf']
ir_irf = np.load('ir_irf_s3_0.7v.npz')['irf']

# Plot
plt.style.use('default')
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['lines.markersize'] = 3
legendsize=14

# Time axis
t = np.arange(150) * 0.065

# Plt
plt.figure()
nir_to_plot = nir_irf[0][1]
nir_to_plot = nir_to_plot / np.max(nir_to_plot)
nir_to_plot[0:8] = 0
ir_to_plot = ir_irf[0][3]
ir_to_plot[0:6] = 0
ir_to_plot = ir_to_plot / np.max(ir_to_plot) *10
# plt.semilogy(nir_to_plot, color='darkred', label='NIR')
# plt.semilogy(ir_to_plot, color='purple', label='IR')
plt.semilogy(t - t[np.argmax(nir_to_plot)], nir_to_plot, color='darkred', label='680 nm')
plt.semilogy(t - t[np.argmax(ir_to_plot)], ir_to_plot, color='purple', label='850 nm')
plt.legend(fontsize=legendsize)
plt.xlabel("Time (ns)")
plt.ylabel("Normalized Counts")
plt.grid(which='both')





