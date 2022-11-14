#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
utils.py

Functions required for the analysis of SMURPHS OHC data.

Required to reproduce data for Boland et al. 2022 (preprint https://www.essoar.org/doi/10.1002/essoar.10511062.3)
See https://github.com/emmomp/SMURPHS_OHC for details

For more details on the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806).

Updated Nov 2022

@author: emmomp@bas.ac.uk Emma J D Boland
"""
from scipy import signal
import numpy as np
import xarray as xr
from scipy import stats
import statsmodels.api as sm
import os
import pyresample

##Create 4th order Bworth filter
def butter_lowpass(data,cut,order=4,sample_freq=1) :
    nyq = 0.5*sample_freq
    pass_freq =1./cut/nyq
    sos=signal.butter(order,pass_freq,'low',output='sos')
    filt=signal.sosfiltfilt(sos,data)
    return filt

def butter_ufunc(data,cut,tdim,order=4,sample_freq=1):
    filt = xr.apply_ufunc(
        butter_lowpass,
        data,
        cut,
        order,
        sample_freq,
        input_core_dims=[[tdim],[],[],[]],
        output_core_dims=[[tdim]],
        vectorize=True,
        dask='parallelized',
        output_dtypes=['float64']
        )
    return filt

def new_linregress(x,y):
    # Wrapper around scipy linregress to use in apply_ufunc
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    return np.array([slope, intercept, r_value, p_value, std_err])

def linregress_2d(x,y):
    y=y.broadcast_like(x).data.ravel()
    x=x.data.ravel()
    mask=np.isfinite([x,y]).all(axis=0)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x[mask],y[mask])
    return np.array([slope, intercept, r_value, p_value, std_err])

def lin_regress(time,ohc,indims,dtype='float64'):    
    stats = xr.apply_ufunc(
    new_linregress,  # first the function
    time,  # now arguments in the order expected by func
    ohc,  # as above
    input_core_dims=indims,  # list with one entry per arg
    output_core_dims=[['parameter']],  # returned data has one dimension
    vectorize=True,  # loop over non-core dims
    dask='parallelized',
    output_dtypes=[dtype],  # one per output
    dask_gufunc_kwargs={'output_sizes':{"parameter": 5}},
#    output_sizes={"parameter": 5},
    )
    stats['parameter']=('parameter',["slope","intercept","r_value","p_value","std_err"])
    return stats

def new_weighted_linregress(x,y,w):
    # Wrapper around scipy linregress to use in apply_ufunc
    mod_wls = sm.WLS(y.data,sm.add_constant(x.data), weights=1.0 / (w.data ** 2))
    res_wls=mod_wls.fit()
    slope=res_wls.params[1]
    intercept=res_wls.params[0]
    r_value=np.sqrt(res_wls.rsquared)
    p_value=res_wls.f_pvalue
    std_err=res_wls.bse[1]
    return np.array([slope, intercept, r_value, p_value, std_err])

def polyfit_xr(ohc_global,nd=4):
    model = xr.apply_ufunc(
    np.polyfit,  # first the function
    ohc_global.time_days, 
    ohc_global.ohc,  
    nd, 
    input_core_dims=[["time"], ["time"], []],  # list with one entry per arg
    output_core_dims=[["coeffs",]],  # returned data has one dimension
    vectorize=True,  # loop over non-core dims
    dask="parallelized",
    output_dtypes=["float64"],  # one per output
    dask_gufunc_kwargs={'output_sizes':{"coeffs": nd+1}}
    )
    return model

def polyval_xr(model,x):
    fit = xr.apply_ufunc(
    np.polyval,  # first the function
    model,  # now arguments in the order expected by 'interp1_np'
    x,  # as above
    input_core_dims=[["coeffs"], ["time"]],  # list with one entry per arg
    output_core_dims=[["time"]],  # returned data has one dimension
    vectorize=True,  # loop over non-core dims
    dask="parallelized",
    output_dtypes=["float64"],  # one per output
    )
    return fit

def setup_regrid(nav_lon,nav_lat,new_lon=np.linspace(-179,180,360),new_lat=np.linspace(-89.5,89.5,180)):
    orig_grid = pyresample.geometry.SwathDefinition(lons=nav_lon.data, lats=nav_lat.data)
    yi,xi=np.meshgrid(new_lat,new_lon)
    new_grid  = pyresample.geometry.GridDefinition(lons=xi,lats=yi)
    resample_data = pyresample.kd_tree.get_neighbour_info(orig_grid, new_grid, 100000, neighbours=1)
    return resample_data   

def repeat_regrid(ds,islogical,resample_data,new_lon=np.linspace(-179,180,360), new_lat=np.linspace(-89.5,89.5,180),loop_dim='exp'):    
    grid_shape=[new_lon.size,new_lat.size]
    if islogical:
        foo = pyresample.kd_tree.get_sample_from_neighbour_info('nn', grid_shape, ds.transpose(...,loop_dim).astype(int).values,
                                              resample_data[0], resample_data[1],resample_data[2])    
        
        foo=np.ma.masked_equal(foo, 0)
        foo[foo>0]=1
    else:  
        foo = pyresample.kd_tree.get_sample_from_neighbour_info('nn', grid_shape, ds.transpose(...,loop_dim).values,
                                              resample_data[0], resample_data[1],resample_data[2])    
    ds2=xr.DataArray(foo,dims=['lon','lat',loop_dim],coords={'lon':(('lon'),new_lon),'lat':(('lat'),new_lat),loop_dim:(loop_dim,ds[loop_dim].data)})
    return ds2