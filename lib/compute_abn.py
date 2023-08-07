#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 12:14:21 2023

@author: c2056366
"""
import numpy as np

# a function to compute the bandpower abnormalities
def computeAbnormalities(controls,patients):

    
    #these create 2D matrices of size nRegions x frequency
   
        
    control_mean=np.mean(controls,2)

    control_std=np.std(controls,2)


    #reshape to create 3D matrices of size nRegions x frequency x number of patients
    control_mean_3D=np.repeat(control_mean[:, :, np.newaxis], patients.shape[2], axis=2)
    control_std_3D=np.repeat(control_std[:, :, np.newaxis], patients.shape[2], axis=2)

    #z-score every element in the patient array
    return (patients-control_mean_3D)/control_std_3D