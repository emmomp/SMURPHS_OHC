# SMURPHS_OHC
[![DOI](https://zenodo.org/badge/461815488.svg)](https://zenodo.org/badge/latestdoi/461815488)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

This repositary contains python code and notebooks to accompany the manuscript "Ocean Heat Content responses to changing Anthropogenic Aerosol Forcing Strength: regional and multi-decadal variability" Boland et al 2022 ([pre-print](https://doi.org/10.1002/essoar.10511062.1)). The contents will allow for the reproduction of all tables and figures in the paper, as well the reproduction of the data files necessary for the figures. See below for more details.

Feel free to use or reproduce the code and figures but please attribute as outlined in the license.

For more details on the SMURPHS ensemble, see Dittus et al. 2020 (https://doi.org/10.1029/2019GL085806)

E Boland Feb 2022 [emmomp@bas.ac.uk](mailto:emmomp@bas.ac.uk)

## Requirements

To reproduce the paper's figures and tables:
- Cartopy==0.17.0
- cftime==1.0.3.4
- matplotlib==3.0.2
- numpy==1.15.4
- scikit_learn==1.0.2
- scipy==1.1.0
- xarray==0.11.0

To reproduce the data from model output:
- cftime==1.0.3.4
- numpy==1.15.4
- pandas==0.23.4
- pyresample==1.16.0
- scipy==1.1.0
- statsmodels==0.9.0
- xarray==0.11.0

## Steps to reproduce the paper's figures and tables

To reproduce the paper's figures and tables, follow these steps:
- Download the data required for the figures from the Figshare repository https://figshare.com/articles/dataset/data_in/19281761/2 and place in a directory named 'data_in'. Which data is required for which figures is listed below. Alternatively this data can be re-generated from the original model output using the python files in the [code](code/) directory - see "Steps to reproduce the paper's analysis from model output".
- Install necessary libraries (see requirements above or [figure_notebooks/requirements,txt](figure_notebooks/requirements.txt)).
- Clone the [figure_notebooks](figure_notebooks/) directory into the same directory that contains 'data_in'.
- Run the notebooks.

### Further Details

To reproduce the figures, the notebooks will look for the following directories/files in a directory called 'data_in':
- Figure 1 & Table 2: ohc_tseries, pic_data
- Figure 2: ohc_tseries, pic_data, other_model_data
- Figure 3: ohc_trends
- Figures 4, S1, S2: ohc_xy
- Figures 5, S3-S6: ohc_yz
- Figure 6: ohc_xy
- Figure 7: ohc_yz
- Figure S7: amoc_tseries
- Figure S8: SIE_SH.nc

## Steps to reproduce the paper's analysis from model output

To reproduce the data files required to produce the figures, follow these steps:

- Access the model output: details coming
- Install necessary libraries (see requirements above or [code/requirements,txt](code/requirements.txt))
- Clone the code directory [code](code/) which containts python code to create the data for the figures. Automatically writes to ../data_in/. See below for further details for which scripts write to which sub-directories.
- Run the python scripts.

### Further details

The data files themselves are written to '../data_in' using the python scripts in [code](code/) as follows, listed by sub-directory:
- ohc_tseries: ohc_by_basin_depth.py
- pic_data: ohc_by_basin_depth_pic.py, ohc_pic_drift.py, ohc_xy_pic_drift.py, ohc_yz_pic_drift.py, ohc_xy_pic.py, ohc_yz_pic.py
- ohc_trends: ohc_weightedtrends_obs.py, ohc_weightedtrends.py
- ohc_xy: ohc_xy.py, ohc_xy_trends.py
- ohc_yz: ohc_yz.py, ohc_yz_trends.py
- amoc_tseries: calculate_AMOC.py
- SIE_SH.nc : calc_SH_SIE.py
