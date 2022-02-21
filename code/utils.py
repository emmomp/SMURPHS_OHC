#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 15:59:43 2020

@author: emmomp
"""
from scipy import signal
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import xarray as xr
from scipy import stats
import cartopy.crs as ccrs
import cartopy.feature as cf
import pyresample
#import nc_time_axis

figs_dir = '../plots/smurphs_ensemble/'

y_start=[1960,1980]
y_end=[1991,2011]
y_mid=[1975.5,1995.5]
depthlabels=['0-700m','0-2000m']

griddata = xr.open_dataset('~/smurphs_ensemble/for_Dan/nemo_at491o_1m_19991201-20000101_grid-T.nc')

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

def polyfit_xr(ohc_global,nd=4):
    model = xr.apply_ufunc(
    np.polyfit,  # first the function
    ohc_global.time_days, 
    ohc_global.ohc,  
    nd, 
    input_core_dims=[["time_centered"], ["time_centered"], []],  # list with one entry per arg
#    output_core_dims=[["coeffs",],["time_centered",]],  # returned data has one dimension
    output_core_dims=[["coeffs",]],  # returned data has one dimension
   # exclude_dims=set(("lat",)),  # dimensions allowed to change size. Must be a set!
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
    input_core_dims=[["coeffs"], ["time_centered"]],  # list with one entry per arg
    output_core_dims=[["time_centered"]],  # returned data has one dimension
    vectorize=True,  # loop over non-core dims
    dask="parallelized",
    output_dtypes=["float64"],  # one per output
    )
    return fit

        
