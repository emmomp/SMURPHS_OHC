#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calculate_volcanic_index.py

Code to calculate a volcanic forcing index time series for the CMIP6 historical period.

Required to reproduce data for Boland et al. 2022 (preprint https://www.essoar.org/doi/10.1002/essoar.10511062.3)
See https://github.com/emmomp/SMURPHS_OHC for details

For the CMIP6 historical HadGEM ensemble, see Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

Updated Nov 2022

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import xarray as xr
from datetime import date

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'CMIP6 historical volcanic forcing from Boland et al 2022 (https://www.essoar.org/doi/10.1002/essoar.10511062.3)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of the SMURPHS ensemble, See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806'}

data_dir='/gws/nopw/j04/cmip6_prep_vol1/cmip6_ancils/data/ancils/n96e/timeseries_1850-2014/VolcanicAod/v2/' # Holding accessible via JASMIN
save_dir = '../data_in/' #Directory to save data to

volc=xr.open_dataset(data_dir+'volc_aer_absorption_lw.nc')
easy=volc['easy_absorption_lw'].sum(dim=['latitude','model_level_number','waveband'])

index_out=easy.copy()
index_out.attrs.update(volc.attrs)
index_out.attrs.update(attrs)
index_out.attrs['description']='Volcanic Aerosol LW Absorption summed over latitude, model level, and waveband'

index_out.to_netcdf(save_dir+'volc_index.nc')