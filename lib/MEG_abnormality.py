#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 09:17:11 2023
This script is designed to calculate  abnormality and  Drs based on MEG bandpower.
1. Load MEG bandpower for Complete power spectrum/PEriodic band power
2. Load Resection and outcome tables
3. Load aperiodic components
4. Calcualte abnormality for each possible power spectrum/exponent
5. Calcualtge ROC_AUC and store metadata info (outcome, ID) and Drs values in a table
6. Save the results to use it in further analysis


@author: Csaba Kozma CNNP lab 2023 
"""

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import scipy.io as sio
import compute_abn

#Complete band power data
ucl_ids=np.load('/Users/c2056366/Documents/Fooof_paper_material/MEG_data/Bp_raw_data/ucl_ids.npy')
ucl_bandpower_complete=np.load('/Users/c2056366/Documents/Fooof_paper_material/MEG_data/Bp_raw_data/ucl_bandpower.npy')
cardiff_bandpower_complete=np.load('/Users/c2056366/Documents/Fooof_paper_material/MEG_data/Bp_raw_data/cardiff_bandpower.npy')

#Load resection and outcome tables
resection_data=pd.read_csv('/Users/c2056366/Documents/Fooof_paper_material/MEG_data/Bp_raw_data/meg_resection_table.csv')
outcome_data=pd.read_csv('/Users/c2056366/Documents/Fooof_paper_material/MEG_data/Bp_raw_data/meg_outcome_table.csv',dtype={'ID':str})


#Load periodic band power data
UCLH=sio.loadmat('/Users/c2056366/Documents/Fooof_paper_material/MEG_data/Bp_raw_data/ucl_bp_flattened_30.mat')
ucl_bandpower_periodic=UCLH['ucl_bp']
Cardiff=sio.loadmat('/Users/c2056366/Documents/Fooof_paper_material/MEG_data/Bp_raw_data/cardiff_bp_flattened_30.mat')
cardiff_bandpower_periodic=Cardiff['cardiff_bp']

#Aperiodic data
cardiff_bandpower_aperiodic=pd.read_csv('/Users/c2056366/Documents/Fooof_paper_material/MEG_data/Bp_raw_data/Cardiff_aperiodic_30hz_MEG.csv')
ucl_bandpower_aperiodic=pd.read_csv('/Users/c2056366/Documents/Fooof_paper_material/MEG_data/Bp_raw_data/UCLH_aperiodic_30hz_MEG.csv')

###########Analyse psd type ############

type_psd='periodic' #(aperiodic, complete, periodic)

if type_psd=='aperiodic':
    "APERIODIC Band Power"
    
    #these create 1D matrices of size nRegions x frequency
    controls=np.array(cardiff_bandpower_aperiodic.aperiodic_cmps_2).reshape(114,88)
    control_mean=np.mean(controls,1)
    control_std=np.std(controls,1)
        
        
    #reshape to create 2D matrices of size nRegions x frequency x number of patients
    patients=np.array(ucl_bandpower_aperiodic.aperiodic_cmps_2).reshape(114,55)
    controls_mean_2D= np.repeat(control_mean[:,  np.newaxis], patients.shape[1], axis=1)
    control_std_2D=np.repeat(control_std[:,  np.newaxis], patients.shape[1], axis=1)
    
    #z-score every element in the patient array
    ucl_z_scores=(patients-controls_mean_2D)/control_std_2D
    #retain the maximum across the frequency bands
    max_abnormalities=np.abs(ucl_z_scores)
    
    #store the drs results in a list
    patient_drs=[]
    patient_outcome=[]
    patient_ids=[]
elif type_psd=='complete':
    
    "COMPLETE Band Power"
    #compute patient abnormalities
    ucl_z_scores=compute_abn.computeAbnormalities(cardiff_bandpower_complete,ucl_bandpower_complete)
    #retain the maximum across the frequency bands
    max_abnormalities=np.max(np.abs(ucl_z_scores),axis=1)
    
    #store the drs results in a list
    patient_drs=[]
    patient_outcome=[]
    patient_ids=[]
else:
    "PERIODIC Band Power"
    #compute patient abnormalities
    ucl_z_scores=compute_abn.computeAbnormalities(cardiff_bandpower_periodic,ucl_bandpower_periodic)
    #retain the maximum across the frequency bands
    max_abnormalities=np.max(np.abs(ucl_z_scores),axis=1)
    
    #store the drs results in a list
    patient_drs=[]
    patient_outcome=[]
    patient_ids=[]

"ABNOTMALITY"
for i,subject in enumerate(ucl_ids):

    #patients to ignore due to missing data
    if str(subject) in ['368','863','984']:
        continue

    #get the bandpower abnormalities
    #ignore the first 14 data points as these are subcortical
    patient_abnormalities=max_abnormalities[:,i]

    #get the patient resection data
    #for simplicity only doing resected/spared. So replace any unknown (-1) to spared (0)
    patient_resection=[1 if x>0 else 0 for x in resection_data[str(subject)][14:]]

    #store the patient drs in the list
    patient_drs.append(1-roc_auc_score(patient_resection,patient_abnormalities))

    #store the patients outcome in the list
    patient_outcome.append(outcome_data[outcome_data['ID']==str(subject)]['outcome'].values[0])

    #store the id

    patient_ids.append(str(subject))
    

#store all the results in a final table
final_results=pd.DataFrame({'ID':patient_ids,'outcome':patient_outcome,'drs':patient_drs})

###Save results for further analysis
final_results.to_csv('/Users/c2056366/Documents/Fooof_paper_material/MEG_data/Bp_raw_data/MEG_auc_cs_flat_30hz.csv', sep='\t')
