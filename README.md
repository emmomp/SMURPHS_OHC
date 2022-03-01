# SMURPHS_OHC

# UNDER DEVELOPMENT

figure_notebooks: Jupyter notebooks to reproduce figures (and one table) from "Ocean Heat Content responses to changing Anthropogenic Aerosol Forcing Strength: regional and multi-decadal variability" Boland et al 2022 (doi coming).

Each notebook loads data from ../data_in/. To work, data_in should contain the files from [figshare details]. Alternatively this data can be re-generated using the python files in the code directory - see below.

code: Python code to create the data for the figures. Automatically writes to ../data_in/. Requires access to raw model output data which can be accessed [reference] and a machine to run the analyis. 

Feel free to use or reproduce the code and figures but please attribute as outlined in the license.

E Boland Feb 2022 [emmomp@bas.ac.uk](mailto:emmomp@bas.ac.uk)

## Further Details

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

The data files themselves were created using the python scripts in [code][code/] as follows:
- ohc_tseries: ohc_by_basin_depth.py
- pic_data: ohc_by_basin_depth_pic.py, ohc_pic_drift.py, ohc_xy_pic_drift.py, ohc_yz_pic_drift.py, ohc_xy_pic.py, ohc_yz_pic.py
- ohc_trends: ohc_weightedtrends_obs.py, ohc_weightedtrends.py
- ohc_xy: ohc_xy.py, ohc_xy_trends.py
- ohc_yz: ohc_yz.py, ohc_yz_trends.py
- amoc_tseries: calculate_AMOC.py
- SIE_SH.nc : calc_SH_SIE.py
