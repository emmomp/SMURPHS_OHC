#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ohc_xy_pic_drift.py

Code to load ocean heat content from the PIC of HadEM3-GC31-LL and calculate 
linear drift to de-drift SMURPHS ensemble. Heat content produced by 
ohc_xy_pic.py

Required to reproduce data for Boland et al. 2023 (https://doi.org/10.1029/2022JC018725)
See https://github.com/emmomp/SMURPHS_OHC for details

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) 
and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

Updated Jun 2023

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import xarray as xr
import utils
import numpy as np
from datetime import date

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'HadGem3-GC31-LL PIC OHC data from Boland et al. 2023 (https://doi.org/10.1029/2022JC018725)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of HadGem3-GC31-LL CMIP6 PIC data, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995'}

save_dir = '../data_in/' #Directory to load from and save data to

ohc_xy=xr.open_dataarray(save_dir+'pic_data/ohc_xy_pic.nc')
nt=ohc_xy['time'].size
ohc_xy['time_mths']=('time',np.arange(0,nt))
drift=utils.lin_regress(ohc_xy.time_mths[0:6000],ohc_xy[0:6000,:,:],[['time',],['time',]])

OHC_drift=drift.sel(parameter='slope')
OHC_drift.name='OHC_xy_drift'
del OHC_drift['parameter']
OHC_drift.attrs['long_name']='Drift in depth integrated OHC from PIC'
OHC_drift.attrs['units']='J/m^2/month'
OHC_drift.attrs.update(attrs)
OHC_drift.to_netcdf(save_dir+'pic_data/ohc_xy_pic_drift.nc')
