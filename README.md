# FOOOF-Based Brain Abnormality Analysis

This repository contains MATLAB and Python scripts for analyzing brain abnormalities using FOOOF (Fitting Oscillations & One Over F) decomposition of power spectral density data from intracranial EEG (iEEG) and Magnetoencephalography (MEG) recordings.

## Overview

This project implements a comprehensive pipeline for:
- Decomposing power spectral density data into periodic (oscillatory) and aperiodic (1/f) components using FOOOF
- Calculating normative maps from control populations
- Identifying brain abnormalities in patient populations through z-score comparisons
- Visualizing abnormalities on brain surfaces
- Analyzing the relationship between abnormalities and clinical outcomes

## Project Structure

```
├── Figure1.m              # Generates normative brain maps
├── Figure2.m              # Individual patient abnormality analysis
├── Figure3.m              # Group-level abnormality analysis
├── Figure4.m              # MEG-specific analysis
├── define_paths.m         # Path configuration script
├── atlas/                 # Brain atlas data
│   ├── ATLAS.mat
│   ├── atlasinfo.mat
│   └── brain_surf/
├── data/                  # Input data files
│   ├── MEG_*.csv         # MEG bandpower data
│   ├── RAM_*.csv         # Control population data
│   ├── UCLH_*.csv        # Patient population data
│   └── metadata.*        # Patient metadata
└── lib/                   # Utility functions
    ├── calc_abnormality.m           # Core abnormality calculation
    ├── calc_abnormality_aperiodic.m # Aperiodic component analysis
    ├── MEG_abnormality.py          # MEG-specific analysis
    ├── plot_abnormalities.m        # Brain visualization
    ├── load_*.m                    # Data loading utilities
    └── [other utility functions]
```

## Key Features

### 1. FOOOF Decomposition
- **Complete Power Spectrum**: Traditional bandpower analysis
- **Periodic Components**: Oscillatory activity after removing 1/f background
- **Aperiodic Components**: 1/f slope exponent analysis

### 2. Abnormality Detection
- Z-score based comparison against normative population
- Support for multiple frequency bands (Delta, Theta, Alpha, Beta)
- Resection-aware analysis for epilepsy patients

### 3. Visualization
- Brain surface mapping of abnormalities
- Individual patient plots
- Group-level statistical analysis
- Custom colormaps for different frequency bands

### 4. Multi-Modal Support
- **iEEG**: High-resolution intracranial recordings
- **MEG**: Non-invasive magnetoencephalography

## Prerequisites

### MATLAB Dependencies
- MATLAB R2019b or later
- [FieldTrip](https://www.fieldtriptoolbox.org/) (for MEG analysis)
- Statistics and Machine Learning Toolbox
- Image Processing Toolbox

### Python Dependencies
- Python 3.7+
- NumPy
- Pandas
- SciPy
- scikit-learn

## Installation

1. Clone this repository:
```bash
git clone [repository-url]
cd [repository-name]
```

2. Download and install FieldTrip:
```bash
# Follow instructions at https://www.fieldtriptoolbox.org/download/
```

3. Update paths in `define_paths.m`:
```matlab
analysis_location = '/path/to/this/repository';
fieldtrip_location = '/path/to/fieldtrip';
```

## Usage

### Basic Workflow

1. **Setup Paths** (required first):
```matlab
define_paths
```

2. **Generate Normative Maps** (Figure 1):
```matlab
Figure1
```

3. **Analyze Individual Patients** (Figure 2):
```matlab
Figure2
```

4. **Group Analysis** (Figure 3):
```matlab
Figure3
```

5. **MEG Analysis** (Figure 4):
```matlab
Figure4
```

### Python MEG Analysis

For MEG-specific abnormality calculations:
```bash
cd lib
python MEG_abnormality.py
```

## Data Format

### Input Data Structure
- **Control Data**: Normalized bandpower values per brain region
- **Patient Data**: Individual patient recordings with metadata
- **Atlas Data**: Brain parcellation schemes (2-region and 4-region)
- **Metadata**: Patient information including outcomes and resection data

### Output Files
- Brain surface visualizations saved to `Figures/` directory
- Statistical results and abnormality scores
- ROC-AUC calculations for clinical outcome prediction

## Configuration

### Parcellation Schemes
- **2-region**: Coarse brain division
- **4-region**: Finer anatomical parcellation

### Frequency Bands
- Delta: 1-4 Hz
- Theta: 4-8 Hz
- Alpha: 8-13 Hz
- Beta: 13-30 Hz

### Analysis Types
- `complete`: Full power spectrum analysis
- `periodic`: Oscillatory components only
- `aperiodic`: 1/f slope analysis

## Citation

If you use this code in your research, please cite:

```
[Add your citation here]
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## License

[Add your license information here]

## Contact

For questions or issues, please contact:
- **Author**: Csaba Kozma
- **Lab**: CNNP Lab, Newcastle University
- **Email**: [c.a.kozma2@newcastle.ac.uk]

## Acknowledgments

- FieldTrip toolbox developers
- FOOOF developers
- Clinical collaborators and data contributors 