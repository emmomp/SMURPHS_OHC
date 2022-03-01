#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ohc_yz_trends.py

Code to load zonally integrated ocean heat content time series and calculate linear trends, as part of analysis for Boland et al (in prep)

Global ocean heat content data from the SMURPHS ensemble, produced using ohc_yz.py
De-drifted by PIC data, produced using ohc_yz_pic_drift.py

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

Created on Wed Dec  9 16:01:17 2020

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import xarray as xr
import utils
import numpy as np
from datetime import date

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'SMURPHS OHC data from Boland et al (in prep)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of the SMURPHS ensemble, See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806'}

save_dir = '../data_in/' #Directory to save data to

exps=['historical0p2','historical0p4','historical0p7','historical1p0','historical1p5']
runs=['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1']

masks=xr.open_dataset(save_dir+'other_model_data/subbasins_eORCA1-GO6-Daley.nc')

basin_masks={
        'so':(masks.nav_lat<=-35.0),
        'atl':masks.tmaskatl.astype(bool),
        'ind':masks.tmaskind.astype(bool),
        'pac':masks.tmaskpac.astype(bool)
            }

basins=list(basin_masks.keys())
basins=basins+['global']

startdate=19550101
enddate=20141231

startdate_str='1955-01-01'
enddate_str='2014-12-30'

trend_startdates=['1960-01-16','1980-01-16']
trend_enddates=['1990-12-16','2010-12-16']
trend_middates=['1975-06-16','1995-06-16']

y_start=[1960,1980]
y_end=[1991,2011]
y_mid=[1975.5,1995.5]
 
foo=[]
for ie,exp in enumerate(exps):
    foo.append([])
    for basin in basins:
        ds=xr.open_mfdataset(save_dir+'ohc_yz/ohc_yz_'+basin+'_'+exp+'_*_'+str(startdate)+'_'+str(enddate)+'.nc',combine='nested',concat_dim='run')
        ds['basin']=(('basin'),[basin,])
        foo[ie].append(ds)
ohc_yz=xr.combine_nested(foo,concat_dim=['exp','basin']) 
ohc_yz['time_mths']=('time_centered',np.arange(0,721)) 

print('ohc data loaded, running linear regressions')

run_fit=[]
all_fit=[]
for tchunk in range(0,2):
    #fit all runs
    tslice=slice(trend_startdates[tchunk],trend_enddates[tchunk])
    xx=ohc_yz['time_mths'].sel(time_centered=tslice)
    yy=ohc_yz['ohc'].sel(time_centered=tslice)
    lin_fit=utils.lin_regress(xx,yy,[["time_centered"], ["time_centered"]])
    lin_fit=lin_fit.assign_coords(time=trend_middates[tchunk])
    run_fit.append(lin_fit)
    #fit run mean
    yy=yy.mean(dim='run')
    lin_fit=utils.lin_regress(xx,yy,[["time_centered"], ["time_centered"]])
    lin_fit=lin_fit.assign_coords(time=trend_middates[tchunk])
    all_fit.append(lin_fit)
    
run_fit=xr.concat(run_fit,concat_dim='time')
all_fit=xr.concat(all_fit,concat_dim='time')

# Load drift
drift=xr.open_dataarray(save_dir+'pic_data/ohc_yz_pic_drift.nc')

# Remove drift
ny=run_fit['y'].size
run_slope= run_fit.isel(y=slice(1,ny-1),parameter=0)-drift
all_slope=all_fit.isel(y=slice(1,ny-1),parameter=0)-drift

print('regressions finished, writing to file')
all_slope.attrs.update(attrs)
all_slope.name='stats'
all_slope.attrs['description']='Linear fit statistics for 30 year sections of SMURPHS 0-700m and 0-2000m OHC time series'
all_slope.attrs['slope_units']='x10^22 J/mth'
run_slope.attrs.update(attrs)
run_slope.name='stats'
run_slope.attrs['description']='Linear fit statistics for 30 year sections of SMURPHS 0-700m and 0-2000m OHC time series'
run_slope.attrs['slope_units']='x10^22 J/mth'

all_slope.to_netcdf(save_dir+'ohc_yz/ohc_yz_trend_runmean.nc')
run_slope.to_netcdf(save_dir+'ohc_yz/ohc_yz_trend_byrun.nc')    

print('all done')