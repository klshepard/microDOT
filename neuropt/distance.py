# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 17:09:43 2023

@author: Dell-User
"""

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