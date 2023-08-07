function [aperiodic_components, flattened_psd, fooof_settings ]=meas_fooof_indiv(pxx, freq, fooof_features, python_env_path)
%% MEAS_FOOOF_FEATURES
%
% This MATLAB wrapper script loads the fooof python environment and
% parametarize periodic and aperiodic components.
%
% Export offset, aperiodic exponent, peak frequencies, aperiodic activity,
% perioidc activity (flattened PSD)
%
% Csaba Kozma
% CNNP Lab, Newcastle University
% July 2023

%% PYTHON wrapper for FOOOF
% Python wrapper code
python_wrapper_code = [

"import matplotlib.pyplot as plt"
"from fooof import FOOOF"
"from fooof.sim.gen import gen_aperiodic"
"import numpy as np"
"import pandas as pd"
"import os"
"fm = FOOOF(peak_width_limits, max_n_peaks, min_peak_height)"
"Offsets_channel = []"
"Slope_channel = []"
"periodic_fit = []"
"fooof_aperiodic_components = {'Offsets': [], 'Slopes': []}"
"fooof_aperiodic_components_output = pd.DataFrame(fooof_aperiodic_components)"
"data=np.asarray(power_spectrum_density)"

"for i in range(power_spectrum_density.shape[0]):"
"   spectrum = data[i,:]"
"   fm.fit(np.array(freqs), spectrum, [freq_range_start, freq_range_end]) #generate the fooof model to run"
"   init_ap_fit = gen_aperiodic(fm.freqs, fm._robust_ap_fit(fm.freqs, fm.power_spectrum)) #fit the slope"
"   init_flat_spec = fm.power_spectrum - init_ap_fit #calculate the flattened/periodic psd"
"   periodic_fit.append(init_flat_spec) #appened all the periodic fits for every channel"
"   Offs_pt = fm.aperiodic_params_[0] #extract offset per channel"
"   Slp_pt = fm.aperiodic_params_[1] #extract aperiodic exponent per channel"
"   Offsets_channel.append(Offs_pt) #appened offsets for every channel"
"   Slope_channel.append(Slp_pt) #appened aperiodic exponents for every channel"

"fooof_aperiodic_components_output['Offsets']=Offsets_channel #generate one dataframe with all offsets and exponent values"
"fooof_aperiodic_components_output['Slopes']=Slope_channel"];
pyenv("Version",python_env_path);
%terminal message

disp('Paramaterizing periodic and aperiodic components')

%Prep frequency range to run FOOOF on
freq_range_start=fooof_features.freq_range(1);
freq_range_end = fooof_features.freq_range(end);

%Run the python code thorugh MATLAB pyrun function
%"fooof_aperiodic_components_output" = aperiodic_components as MATLAB
%output
% "periodic_fit" = flattened_psd as MATLAB output
[aperiodic_components,flattened_psd]=pyrun(python_wrapper_code,["fooof_aperiodic_components_output" "periodic_fit"],power_spectrum_density=pxx,freqs=freq,freq_range_start=freq_range_start,...
freq_range_end=freq_range_end,peak_width_limits=fooof_features.peak_width_limits,...
max_n_peaks=fooof_features.max_n_peaks,min_peak_height=fooof_features.min_peak_height);

% transform the python variable to MATLAB structure
aperiodic_components = df2t(aperiodic_components);
aperiodic_components=aperiodic_components{:,:};

% save FOOOF settings
fooof_settings = struct();
fooof_settings.peak_width_limits=fooof_features.peak_width_limits;
fooof_settings.max_n_peaks=fooof_features.max_n_peaks;
fooof_settings.min_peak_height=fooof_features.min_peak_height;
fooof_settings.freq_range = fooof_features.freq_range;
end