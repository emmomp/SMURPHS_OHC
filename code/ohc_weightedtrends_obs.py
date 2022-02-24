#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ohc_weightedtrends_obs.py

Code to load observed ocean heat content time series and calculate weighted
linear trends, as part of analysis for Boland et al (in prep)

Observation file IAP_OHC_estimate_update300m_700m_2000m.txt downloaded from http://159.226.119.60/cheng/ in Jan 2021
see Cheng et al 2017 (https://doi.org/10.1126/sciadv.1601545) for details.

Observation files heat_content_anomaly_0-700_yearly.nc and heat_content_anomaly_0-2000_pentad.nc 
downloaded from https://www.ncei.noaa.gov/products/world-ocean-atlas in Jul 2021
see Boyer et al 2018 (https://accession.nodc.noaa.gov/NCEI-WOA18) for details.

See https://github.com/emmomp/SMURPHS_OHC for details

Created on Tue Feb 22 10:51:20 2022

@author: emmomp
"""
import xarray as xr
import numpy as np
import pandas as pd
import utils
from datetime import date

save_dir = '../data_in/' #Directory to load from and save data to

attrs_iap={'contact':'emmomp@bas.ac.uk',
       'references':'Analysed OHC data from Boland et al (in prep)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of IAP data, see http://www.ocean.iap.ac.cn/'}

attrs_noaa={'contact':'emmomp@bas.ac.uk',
       'references':'Analysed OHC data from Boland et al (in prep)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis WOA 2018 data, see https://www.ncei.noaa.gov/products/world-ocean-atlas'}

print('loading observations')

# Open IAP observation data
obs_global_file=save_dir+'observations/IAP_OHC_estimate_update300m_700m_2000m.txt'
f=open(obs_global_file,'r')
datax=f.readlines()
f.close()
names=datax[17].split('] [')
names[0]=names[0].strip('[')
names[-1]=names[-1][:-2]
names[4]=names[2]+' Error'
names[7]=names[5]+' Error'
names[10]=names[8]+' Error'
obs_global = np.loadtxt(obs_global_file,delimiter='\t',skiprows=18,usecols=range(2,11))
varlist=list(zip(['time',]*9,obs_global.T))
dateindex=pd.date_range('1939-12-01',freq='M',periods=972)+pd.DateOffset(days=15)
ohc_iap=xr.Dataset(coords={'time':('time',dateindex)},data_vars=dict(zip(names[2:],varlist)))
ohc_iap=ohc_iap.sel(time=slice(None,'2014-12-31'))
dateindex=dateindex[:ohc_iap.time.shape[0]]
ohc_iap['time_mths']=('time',np.arange(0,ohc_iap.time.shape[0]))
ohc_iap['time_yrs']=dateindex.year+ (dateindex.month-1)/12.

# Load NOAA observation data
obs_noaa_700=xr.open_dataset(save_dir+'observations/heat_content_anomaly_0-700_yearly.nc',decode_times=False)
newt=np.zeros(66,dtype='datetime64[M]')
for t in range(0,66):
    newt[t]=np.datetime64('1955-01')+np.timedelta64(obs_noaa_700.time.data[t].astype(int),'M')
obs_noaa_700['datetime']=(('time'),newt)   

obs_noaa_2000=xr.open_dataset(save_dir+'observations/heat_content_anomaly_0-2000_pentad.nc',decode_times=False)
newt=np.zeros(60,dtype='datetime64[M]')
for t in range(0,60):
    newt[t]=np.datetime64('1955-01')+np.timedelta64(obs_noaa_2000.time.data[t].astype(int),'M')
obs_noaa_2000['datetime']=(('time'),newt)   

print('Data loaded, performing regressions')

# Perform weighted regressions

# IAP 
# One fit 1955-2015
long_fit_iap=[]
tslice=slice('1955-01-01',None)
xx=ohc_iap['time_mths'].sel(time=tslice)
for drange in ['OHC0-700m','OHC0-2000m']:
    yy=ohc_iap[drange].sel(time=tslice)
    yy_err=ohc_iap[drange+' Error'].sel(time=tslice)
    weights=1/(yy_err**2)
    lstats=utils.new_weighted_linregress(xx,yy,yy_err)
    fit=xr.DataArray(data=lstats,dims=('parameter',),coords={'parameter':["slope","intercept","r_value","p_value","std_err"],'drange':drange}) 
    long_fit_iap.append(fit)
long_fit_iap=xr.concat(long_fit_iap,'drange')
long_fit_iap.attrs.update(attrs_iap)  
long_fit_iap.attrs['description']='Linear fit statistics for 1955-2015 for IAP 0-700m and 0-2000m OHC time series'
long_fit_iap.attrs['slope_units']='x10^22 J/mth'
long_fit_iap.name='stats'
long_fit_iap.to_netcdf(save_dir+'ohc_trends/longfit_obs_iap.nc')

# 30 year running trends
all_fit_iap=[]
for tchunk in range(0,ohc_iap.time.shape[0]-12*31):
    xx=ohc_iap.time_mths[tchunk:tchunk+12*31]   
    tslice=slice(tchunk,tchunk+12*31)
    all_fit_iap.append([])
    for drange in ['OHC0-700m','OHC0-2000m']:
        yy=ohc_iap[drange].isel(time=tslice)
        yy_err=ohc_iap[drange+' Error'].isel(time=tslice)
        lstats=utils.new_weighted_linregress(xx,yy,yy_err)
        fit=xr.DataArray(data=lstats,dims=('parameter',),coords={'parameter':["slope","intercept","r_value","p_value","std_err"],'drange':drange})
        fit=fit.assign_coords(time_yrs=ohc_iap.time_yrs[tchunk].data+15.5)
        all_fit_iap[tchunk].append(fit)
all_fit_iap=xr.combine_nested(all_fit_iap,concat_dim=['time_yrs','drange'])
all_fit_iap.attrs.update(attrs_iap)  
all_fit_iap.attrs['description']='Linear fit statistics for running 30 year sections of IAP OHC time series'
all_fit_iap.attrs['slope_units']='x10^22 J/mth'
all_fit_iap.name='stats'
all_fit_iap.to_netcdf(save_dir+'ohc_trends/allfit_obs_iap.nc')

# NOAA
# One fit 1955-2020/2016

xx=obs_noaa_700['time']
yy=obs_noaa_700['yearl_h22_WO']
yy_err=obs_noaa_700['yearl_h22_se_WO']
lstats=utils.new_weighted_linregress(xx,yy,yy_err)
long_fit_noaa_700=xr.DataArray(data=lstats,dims=('parameter',),coords={'parameter':["slope","intercept","r_value","p_value","std_err"],'drange':'OHC0-700m'}) 

xx=obs_noaa_2000['time']
yy=obs_noaa_2000['pent_h22_WO']
yy_err=obs_noaa_2000['pent_h22_se_WO']
lstats=utils.new_weighted_linregress(xx,yy,yy_err)
long_fit_noaa_2000=xr.DataArray(data=lstats,dims=('parameter',),coords={'parameter':["slope","intercept","r_value","p_value","std_err"],'drange':'OHC0-2000m'}) 

long_fit_noaa=xr.concat([long_fit_noaa_700,long_fit_noaa_2000],'drange')
long_fit_noaa.attrs.update(attrs_noaa)  
long_fit_noaa.attrs['description']='Linear fit statistics for WOA OHC time-series: 1955-2020 for 0-700m, 1957-2016 for 0-2000m'
long_fit_noaa.attrs['slope_units']='x10^22 J/mth'
long_fit_noaa.name='stats'
long_fit_noaa.to_netcdf(save_dir+'ohc_trends/longfit_obs_noaa.nc')

all_fit_noaa=[]
for tchunk in range(0,obs_noaa_700.time.shape[0]-30):
    xx=obs_noaa_700['time'][tchunk:tchunk+30]
    yy=obs_noaa_700['yearl_h22_WO'][tchunk:tchunk+30]
    yy_err=obs_noaa_700['yearl_h22_se_WO'][tchunk:tchunk+30]
    lstats=utils.new_weighted_linregress(xx,yy,yy_err)
    lin_fit=xr.DataArray(data=lstats,dims=('parameter',),coords={'parameter':["slope","intercept","r_value","p_value","std_err"]})
    lin_fit=lin_fit.assign_coords(time_yrs=1955+obs_noaa_700['time'][tchunk+15]/12)
    all_fit_noaa.append(lin_fit)
all_fit_noaa_700=xr.concat(all_fit_noaa,'time_yrs')
all_fit_noaa_700=all_fit_noaa_700.assign_coords(drange='OHC0-700m')

all_fit_noaa=[]
for tchunk in range(0,obs_noaa_2000.time.shape[0]-30):
    xx=obs_noaa_2000['time'][tchunk:tchunk+30]
    yy=obs_noaa_2000['pent_h22_WO'][tchunk:tchunk+30]
    yy_err=obs_noaa_2000['pent_h22_WO'][tchunk:tchunk+30]
    lstats=utils.new_weighted_linregress(xx,yy,yy_err)
    lin_fit=xr.DataArray(data=lstats,dims=('parameter',),coords={'parameter':["slope","intercept","r_value","p_value","std_err"]})
    lin_fit=lin_fit.assign_coords(time_yrs=1955+obs_noaa_2000['time'][tchunk+15]/12)
    all_fit_noaa.append(lin_fit)
all_fit_noaa_2000=xr.concat(all_fit_noaa,'time_yrs')
all_fit_noaa_2000=all_fit_noaa_2000.assign_coords(drange='OHC0-2000m')

all_fit_noaa=xr.concat([all_fit_noaa_700,all_fit_noaa_2000],'drange')
all_fit_noaa.attrs.update(attrs_noaa)  
all_fit_noaa.attrs['description']='Linear fit statistics for running 30 year sections of WOA 2018 OHC time-series'
all_fit_noaa.attrs['slope_units']='x10^22 J/mth'
all_fit_noaa.name='stats'
all_fit_noaa.to_netcdf(save_dir+'ohc_trends/allfit_obs_noaa.nc')

print('all done')