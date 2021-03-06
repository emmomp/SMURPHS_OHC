#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ohc_xy_trends.py

Code to load depth integrated ocean heat content time series and calculate linear trends, as part of analysis for Boland et al (in prep)

Global ocean heat content data from the SMURPHS ensemble, produced using ohc_xy.py
De-drifted by PIC data, produced using ohc_xy_pic_drift.py

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

Created on Wed Nov 24 12:23:15 2021

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import xarray as xr
import numpy as np
import utils
from datetime import date

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'SMURPHS OHC data from Boland et al (in prep)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of the SMURPHS ensemble, See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806'}

save_dir = '../data_in/' #Directory to save data to

exps=['historical0p2','historical0p4','historical0p7','historical1p0','historical1p5']
runs=['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1']

startdate=19550101
enddate=20141231

trend_startdates=['1960-01-16','1980-01-16']
trend_enddates=['1990-12-16','2010-12-16']
trend_middates=['1975-06-16','1995-06-16']

y_start=[1960,1980]
y_end=[1991,2011]
y_mid=[1975.5,1995.5]
    
foo=[]
for exp in exps:
    ds=xr.open_mfdataset(save_dir+'ohc_xy/ohc_xy_'+exp+'_*_'+str(startdate)+'_'+str(enddate)+'.nc',combine='nested',concat_dim='run')
    foo.append(ds)
ohc_xy=xr.concat(foo,'exp') 
ohc_xy['time_mths']=('time_centered',np.arange(0,721))
print('ohc data loaded, running linear regressions')

run_fit=[]
all_fit=[]
for tchunk in range(0,2):
    #fit all runs
    tslice=slice(trend_startdates[tchunk],trend_enddates[tchunk])
    xx=ohc_xy['time_mths'].sel(time_centered=tslice)
    yy=ohc_xy['ohc'].sel(time_centered=tslice)
    lin_fit=utils.lin_regress(xx,yy,[["time_centered"], ["time_centered"]])
    lin_fit=lin_fit.assign_coords(time=trend_middates[tchunk])
    run_fit.append(lin_fit)
    #fit run mean
    yy=yy.mean(dim='run')
    lin_fit=utils.lin_regress(xx,yy,[["time_centered"], ["time_centered"]])
    lin_fit=lin_fit.assign_coords(time=trend_middates[tchunk]) 
    all_fit.append(lin_fit)
    
run_fit=xr.concat(run_fit,'time')
all_fit=xr.concat(all_fit,'time')

# Load drift
drift=xr.open_dataarray(save_dir+'pic_data/ohc_xy_pic_drift.nc')
drift=drift.rename({'j':'y','i':'x','latitude':'nav_lat','longitude':'nav_lon'})

# Remove drift
ny=run_fit['y'].size
nx=run_fit['x'].size
run_fit = run_fit.isel(y=slice(1,ny-1),x=slice(1,nx-1))-drift
all_fit=all_fit.isel(y=slice(1,ny-1),x=slice(1,nx-1))-drift

print('regressions finished, regridding for plots')
resample_data=utils.setup_regrid(run_fit.nav_lon,run_fit.nav_lat)

for tchunk in range(0,2):
    print('Trends '+str(y_start[tchunk])+' '+str(y_end[tchunk]))
    all_fit_regrid=utils.repeat_regrid(all_fit.isel(parameter=0,time=tchunk),False,resample_data)
    all_fit_regrid['time']=all_fit['time'][tchunk]
    all_fit_regrid.attrs.update(attrs)
    all_fit_regrid.name='stats'  
    all_fit_regrid.attrs['description']='Linear fit statistics for 30 year sections of SMURPHS 0-700m and 0-2000m OHC time series'
    all_fit_regrid.attrs['slope_units']='x10^22 J/m^2/mth'
    all_fit_regrid.to_netcdf(save_dir+'ohc_xy/ohc_xy_trend_regrid_runmean_'+str(y_start[tchunk])+str(y_end[tchunk])+'.nc')
    
    run_fit_regrid=[]
    for run in runs:
        dplot=utils.repeat_regrid(run_fit.sel(run=run).isel(parameter=0,time=tchunk),False,resample_data)
        dplot=dplot.assign_coords(run=run)
        dplot['time']=run_fit['time'][tchunk]
        run_fit_regrid.append(dplot)            
    run_fit_regrid=xr.concat(run_fit_regrid,'run')
    run_fit_regrid.attrs.update(attrs)
    run_fit_regrid.name='stats'
    run_fit_regrid.attrs['description']='Linear fit statistics for 30 year sections of SMURPHS 0-700m and 0-2000m OHC time series'
    run_fit_regrid.attrs['slope_units']='x10^22 J/m^2/mth'
    run_fit_regrid.to_netcdf(save_dir+'ohc_xy/ohc_xy_trend_regrid_byrun_'+str(y_start[tchunk])+str(y_end[tchunk])+'.nc')
                      
print('all done')

