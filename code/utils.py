#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 15:59:43 2020

@author: emmomp
"""
from scipy import signal
import numpy as np
import xarray as xr
from scipy import stats
import statsmodels.api as sm
import os

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

# From Scott Hosking github.com/scotthosking/baspy/examples/
def get_file_date_ranges(fnames, filename_structure):
    ### Get start and end dates from file names
    ind = filename_structure.split('_').index('StartDate-EndDate')
    start_dates, end_dates = np.array([]), np.array([])

    for fname in list(fnames):

        fname = os.path.splitext(fname)[0] # rm extention

        ### Is this file time-varying?
        if '_fx_' in fname:
            ### Fixed variable (e.g., land-mask)
            start_date, end_date = 0, 0

        else:
            ### Time-varying (e.g., temperature)
            date_str = fname.split('_')[ind].split('-')

            ### Interpreting the dates
            if len(date_str) == 2:
                ### e.g.,'19900101-20000101'
                start_date, end_date = int(date_str[0]), int(date_str[1])

            elif (len(date_str) == 3):
                if (date_str[2] == 'clim'):
                    ### e.g.,'186001-187912-clim'
                    ### see /badc/cmip5/data/cmip5/output1/MOHC/HadGEM2-ES/piControl/mon/atmos/Amon/r1i1p1/latest/pfull
                    start_date, end_date = int(date_str[0]+'01'), int(date_str[1]+'31')
                else:
                    print('Cannot identify dates '+fname)

            elif (len(date_str) == 1) & (int(date_str[0]) >= 1800) & (int(date_str[0]) <= 2300):
                ### e.g., '1990'
                ### see /badc/cmip5/data/cmip5/output1/ICHEC/EC-EARTH/amip/subhr/atmos/cfSites/r3i1p1/latest/ccb
                start_date, end_date = int(date_str[0]), int(date_str[0])

            else:
                ### Can't define date
                ###### To do: if no date_range then get from ncdump? !!
                print('Cannot identify dates '+fname)
                start_date, end_date = np.nan, np.nan

        start_dates = np.append( start_dates, start_date )
        end_dates   = np.append( end_dates,   end_date )   

    return start_dates, end_dates

def get_date_range_files(filename_structure,netcdfs,start_date,end_date):

    ### List file names that appear to be within our date range
    if (start_date != None) | (end_date != None):
        file_date_range = get_file_date_ranges(netcdfs, filename_structure)
        if (start_date == None): start_date = np.min(file_date_range[0])
        if (end_date == None):   end_date   = np.max(file_date_range[1])
    
        #keep_ind1 = np.where( (end_date >= file_date_range[0])   & (end_date <= file_date_range[1])   )[0]
        keep_ind1 = np.where( (end_date >= file_date_range[0])   & (start_date <= file_date_range[0])   )[0]
        keep_ind2 = np.where( (end_date >= file_date_range[1]) & (start_date <= file_date_range[1]) )[0]
        #keep_ind2 = np.where( (start_date >= file_date_range[0]) & (start_date <= file_date_range[1]) )[0]
        keep_ind  = np.sort(np.unique(np.append(keep_ind1,keep_ind2)))
        
        if (len(keep_ind) == 0): 
                raise ValueError("No netcdf files within requested date range")
        else:
                netcdfs = np.array(netcdfs)[keep_ind]
    
    return netcdfs 
