#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ohc_yz_pic_drift.py

Code to load ocean heat content from the PIC of HadEM3-GC31-LL and calculate linear drift to de-drift SMURPHS ensemble. Heat content produced by  ohc_yz_pic.py

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

See https://github.com/emmomp/SMURPHS_OHC for details

Created on Tue Nov 30 13:59:56 2021

@author: emmomp
"""
import xarray as xr
import utils
import numpy as np
from datetime import date

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'HadGem3-GC31-LL PIC OHC data from Boland et al (in prep)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of HadGem3-GC31-LL CMIP6 PIC data, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995'}

save_dir = '../data_in/' #Directory to load from and save data to

basins=['global','so','atl','ind','pac']
print('loading yz pic data')
ohc_yz=[]
for basin in basins:
    ohc_yz_basin=xr.open_dataarray(save_dir+'pic_data/ohc_yz_{}_pic.nc'.format(basin))
    ohc_yz_basin['basin']=basin
    ohc_yz.append(ohc_yz_basin.isel(time=slice(0,6000)))
ohc_yz=xr.concat(ohc_yz,'basin')
print('data loaded, calculating drifts')
ohc_yz['time_mths']=('time',np.arange(0,ohc_yz['time'].size))
drift=utils.lin_regress(ohc_yz.time_mths,ohc_yz,[['time',],['time',]])
print('drifts calculated, writing to file')
OHC_drift=drift.sel(parameter='slope')
OHC_drift.name='OHC_yz_drift'
del OHC_drift['parameter']
OHC_drift.attrs['long_name']='Drift in zonally integrated OHC from PIC'
OHC_drift.attrs['units']='J/m^2/month'
OHC_drift.attrs.update(attrs)
OHC_drift.to_netcdf(save_dir+'pic_data/ohc_yz_pic_drift.nc')
print('done')
