#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ohc_yz_trends.py

Code to load zonally integrated ocean heat content time series and calculate linear trends.
Global ocean heat content data from the SMURPHS ensemble, produced using ohc_yz.py
De-drifted by PIC data, produced using ohc_yz_pic_drift.py

Required to reproduce data for Boland et al. 2023 (https://doi.org/10.1029/2022JC018725)
See https://github.com/emmomp/SMURPHS_OHC for details

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

Updated Jun 2023

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import xarray as xr
import utils
import numpy as np
from datetime import date

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'SMURPHS OHC data from Boland et al. 2023 (https://doi.org/10.1029/2022JC018725)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of the SMURPHS ensemble, See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806'}

save_dir = '../data_in/' #Directory to save data to

exps=['hist-0p2','hist-0p4','hist-0p7','hist-1p0','hist-1p5']
runs=['r1i1p1f2','r2i1p1f2','r3i1p1f2','r4i1p1f2','r5i1p1f2']
basins=['global','so','pac','atl','ind']

startdate='1955-01-01'
enddate='2015-01-01'

trend_startdates=['1960-01-16','1980-01-16']
trend_enddates=['1990-12-16','2010-12-16']
trend_middates=['1975-06-16','1995-06-16']
 
foo=[]
for ie,exp in enumerate(exps):
    foo.append([])
    for basin in basins:
        ds=xr.open_mfdataset(save_dir+'ohc_yz/ohc_yz_'+basin+'_'+exp+'_*_'+str(startdate)+'_'+str(enddate)+'.nc',combine='nested',concat_dim='run')
        ds['basin']=(('basin'),[basin,])
        foo[ie].append(ds)
ohc_yz=xr.combine_nested(foo,concat_dim=['exp','basin']) 
ohc_yz['time_mths']=('time',np.arange(0,720)) 

print('ohc data loaded, running linear regressions')

run_fit=[]
all_fit=[]
for tchunk in range(0,2):
    #fit all runs
    tslice=slice(trend_startdates[tchunk],trend_enddates[tchunk])
    xx=ohc_yz['time_mths'].sel(time=tslice)
    yy=ohc_yz['ohc'].sel(time=tslice)
    lin_fit=utils.lin_regress(xx,yy,[["time"], ["time"]])
    lin_fit=lin_fit.assign_coords(time=trend_middates[tchunk])
    run_fit.append(lin_fit)
    #fit run mean
    yy=yy.mean(dim='run')
    lin_fit=utils.lin_regress(xx,yy,[["time"], ["time"]])
    lin_fit=lin_fit.assign_coords(time=trend_middates[tchunk])
    all_fit.append(lin_fit)
    
run_fit=xr.concat(run_fit,'time')
all_fit=xr.concat(all_fit,'time')

# Load drift
drift=xr.open_dataarray(save_dir+'pic_data/ohc_yz_pic_drift.nc')

# Remove drift
run_slope=run_fit.loc[dict(parameter='slope')]-drift
all_slope=all_fit.loc[dict(parameter='slope')]-drift

print('regressions finished, writing to file')
all_slope.attrs.update(attrs)
all_slope.name='stats'
all_slope.attrs['description']='Linear fit slope for 30 year sections of SMURPHS 0-700m and 0-2000m OHC time series'
all_slope.attrs['units']='x10^22 J/m^2/mth'
run_slope.attrs.update(attrs)
run_slope.name='stats'
run_slope.attrs['description']='Linear fit slope for 30 year sections of SMURPHS 0-700m and 0-2000m OHC time series'
run_slope.attrs['units']='x10^22 J/m^2/mth'

all_slope.to_netcdf(save_dir+'ohc_yz/ohc_yz_trend_runmean.nc')
run_slope.to_netcdf(save_dir+'ohc_yz/ohc_yz_trend_byrun.nc')    

print('all done')
