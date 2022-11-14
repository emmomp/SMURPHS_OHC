#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ohc_pic_drift.py

Code to load ocean heat content from the PIC of HadEM3-GC31-LL and calculate 
linear drift to de-drift SMURPHS ensemble. Heat content produced by 
ohc_by_basin_depth_pic.py

Required to reproduce data for Boland et al. 2022 (preprint https://www.essoar.org/doi/10.1002/essoar.10511062.3)
See https://github.com/emmomp/SMURPHS_OHC for details

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) 
and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

Updated Nov 2022

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import xarray as xr
import .utils
import numpy as np
from datetime import date

basins=['global','atl','pac','so','ind']

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'HadGem3-GC31-LL PIC OHC data from Boland et al 2022 (https://www.essoar.org/doi/10.1002/essoar.10511062.3)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of HadGem3-GC31-LL CMIP6 PIC data, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995'}

save_dir = '../data_in/' #Directory to load from and save data to

pic_all=[]
for basin in basins:
    pic=xr.open_dataset('{}pic_data/ohc_pic_bydepth_{}.nc'.format(save_dir,basin))
    pic_fd=xr.open_dataset('{}pic_data/ohc_pic_{}.nc'.format(save_dir,basin))
    pic_fd=pic_fd.assign_coords({'lev_bins':'Full Depth'})
    pic=xr.concat([pic_fd,pic],'lev_bins')
    pic['basin']=basin
    pic_all.append(pic)
pic_all=xr.concat(pic_all,'basin')

pic_all['time_months']=(('time',),np.arange(0,6000))
drift=utils.lin_regress(pic_all.time_months,pic_all.ohc/1e18,[['time',],['time',]])

OHC_drift=drift.sel(parameter='slope')
OHC_drift.name='OHC_drift'
del OHC_drift['parameter']
OHC_drift.attrs['long_name']='Drift in OHC from PIC'
OHC_drift.attrs['units']='1e18 J/m**3/month'
OHC_drift.attrs.update(attrs) 
OHC_drift.to_netcdf(save_dir+'pic_data/ohc_pic_all_drift.nc')
