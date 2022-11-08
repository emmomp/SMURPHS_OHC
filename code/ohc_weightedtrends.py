#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ohc_weightedtrends.py

Code to load global ocean heat content time series and calculate weighted linear trends.
Global ocean heat content data from the SMURPHS ensemble, produce using ohc_by_basin_depth.py

Required to reproduce data for Boland et al. 2022 (preprint https://www.essoar.org/doi/10.1002/essoar.10511062.3)
See https://github.com/emmomp/SMURPHS_OHC for details

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

Updated Nov 2022

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import xarray as xr
import numpy as np
import utils
from datetime import date

save_dir = '../data_in/' #Directory to load from and save data to

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'Analysed OHC data from Boland et al 2022 (https://www.essoar.org/doi/10.1002/essoar.10511062.3)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of the SMURPHS ensemble, See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806'}

exps=['hist-0p2','hist-0p4','hist-0p7','hist-1p0','hist-1p5']

volcello=xr.open_dataset(data_dir+'/volcello/gn/v20190628/volcello_Omon_HadGEM3-GC31-LL_piControl_r1i1p1f1_gn_185001-189912.nc')
vol=volcello['volcello'][0]

# Open model ohc time series
ohc_global=[]
for exp in exps:
    ohc_global.append(xr.open_mfdataset(save_dir+'ohc_tseries/ohc_bydepth_'+exp+'*global*.nc',concat_dim='run',combine='nested',data_vars=['ohc',]))
ohc_global=xr.combine_nested(ohc_global,concat_dim='exp',data_vars=['ohc',]) 
ohc_global['time_mths']=('time',np.arange(0,1980))

# Open drift from PIC experiment
ohc_pic_drift=xr.open_dataarray(save_dir+'pic_data/ohc_pic_all_drift.nc')
time_months=xr.DataArray(dims=['time',],data=np.arange(0,1980))

# Remove drift
ohc_global['ohc']=(ohc_global['ohc']-(ohc_pic_drift.sel(basin='global')*1e18*time_months))/1e22

# Create two summed ranges
ohc_0700=ohc_global.ohc.sel(lev_bins=['0-300m','300-700m']).sum(dim='lev_bins')
ohc_02000=ohc_global.ohc.sel(lev_bins=['0-300m','300-700m','700-2000m']).sum(dim='lev_bins')

# Multiply by vols to retrieve total OHC
ohc_0700=ohc_0700*volcello.sel(lev=slice(None,700)).sum()
ohc_02000=ohc_02000*volcello.sel(lev=slice(None,2000)).sum()

# Linear fit to 1955-2015
long_fit_0700=utils.lin_regress(ohc_global.time_mths.sel(time=slice('1955-01-01',None)),ohc_0700.sel(time=slice('1955-01-01',None)),[["time"], ["time"]])
long_fit_0700['drange']='OHC0-700m'
long_fit_02000=utils.lin_regress(ohc_global.time_mths.sel(time=slice('1955-01-01',None)),ohc_02000.sel(time=slice('1955-01-01',None)),[["time"], ["time"]])
long_fit_02000['drange']='OHC0-2000m'
long_fit=xr.concat([long_fit_0700,long_fit_02000],'drange')

long_fit.attrs.update(attrs)  
long_fit.name='stats'
long_fit.attrs['description']='Linear fit statistics for 1955-2015 SMURPHS 0-700m and 0-2000m OHC time series'
long_fit.attrs['slope_units']='x10^22 J/mth'
long_fit.to_netcdf(save_dir+'ohc_trends/longfit_model.nc')

# Running 30 year fits

all_fit=[]
for tchunk in range(0,1608):
    lin_fit=utils.lin_regress(ohc_global.time_mths[tchunk:tchunk+12*31],ohc_0700.isel(time=slice(tchunk,tchunk+12*31)),[["time"], ["time"]])
    lin_fit=lin_fit.assign_coords(time_mths=tchunk)
    all_fit.append(lin_fit)
all_fits_0700=xr.combine_nested(all_fit,concat_dim='time_mths')
all_fits_0700['time_years']=1850+(all_fits_0700['time_mths']+15.5*12)/12
all_fits_0700=all_fits_0700.swap_dims({'time_mths':'time_years'})
all_fits_0700['drange']='OHC0-700m'

all_fit=[]
for tchunk in range(0,1608):
    lin_fit=utils.lin_regress(ohc_global.time_mths[tchunk:tchunk+12*31],ohc_02000.isel(time=slice(tchunk,tchunk+12*31)),[["time"], ["time"]])
    lin_fit=lin_fit.assign_coords(time_mths=tchunk)
    all_fit.append(lin_fit)
all_fits_02000=xr.combine_nested(all_fit,concat_dim='time_mths')
all_fits_02000['time_years']=1850+(all_fits_02000['time_mths']+15.5*12)/12
all_fits_02000=all_fits_02000.swap_dims({'time_mths':'time_years'})
all_fits_02000['drange']='OHC0-2000m'

all_fits=xr.concat([all_fits_0700,all_fits_02000],'drange')

all_fits.attrs.update(attrs)  
all_fits.name='stats'
all_fits.attrs['description']='Linear fit statistics for running 30 year sections of SMURPHS 0-700m and 0-2000m OHC time series'
all_fits.attrs['slope_units']='x10^22 J/mth'
all_fits.to_netcdf(save_dir+'ohc_trends/allfit_model.nc')
