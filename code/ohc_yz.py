#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ohc_yz.py

Code to load ocean temperature data from the SMURPHS ensemble model output and calculate time series of zonally-integrated Ocean Heat Content from 1955-2014. 

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

See https://github.com/emmomp/SMURPHS_OHC for details

Created on Mon Nov  8 16:33:35 2021

@author: emmomp
"""
import baspy as bp
import xarray as xr
import date_range
from datetime import date

rho_0 = 1.027e3 
c_p = 3850

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'SMURPHS OHC data from Boland et al (in prep)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of the SMURPHS ensemble, See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806'}

save_dir = '../data_in/' #Directory to save data to

griddata = xr.open_dataset(save_dir+'other_model_data/nemo_grid-T.nc')
dx=griddata.e1t

# Generate dataframe with all temp data
df = bp.catalogue(dataset='smurphs',Var=['thetao'])
nt = len(bp.get_files(df.iloc[0]))

exps = set(df.Experiment)
runs = set(df.RunID)

masks=xr.open_dataset(save_dir+'other_model_data/subbasins_eORCA1-GO6-Daley.nc')

basin_masks={
        'so':(masks.nav_lat<=-35.0),
        'atl':masks.tmaskatl.astype(bool),
        'ind':masks.tmaskind.astype(bool),
        'pac':masks.tmaskpac.astype(bool)
            }

fname_struct='Var_Model_Experiment_RunID_Frequency_StartDate-EndDate_grid'
startdate=19550101
enddate=20141231

startdate_str='1955-01-01'
enddate_str='2014-12-30'

def preproc_data(ds):
    ds=ds.swap_dims({'time_counter':'time_centered'})
    return ds.drop('deptht_bounds')

for index, row in df.iterrows():
    exp = row[0]
    run = row[1]
    print('Opening row '+exp+' '+run)
    files = bp.get_files(row)    
    sfiles = date_range.get_date_range_files(fname_struct,files,startdate,enddate)
    print('Found ',len(sfiles),' files ')

    with xr.open_mfdataset(sfiles,preprocess=preproc_data,concat_dim='time_centered',combine='nested',data_vars='minimal', coords='minimal', compat='override') as data:
        print('Got data, calculating ohc ')
        data_weighted = data.thetao*dx
        ohc_yz=data_weighted.sum(dim='x')*rho_0*c_p
        ohc_yz.coords['lat']=data.nav_lat.mean(dim='x')
        ohc_yz['exp']=exp
        ohc_yz['run']=run
        ohc_yz.name='ohc'
        ohc_yz.attrs['long_name']='Ocean Heat Content, zonally integrated'
        ohc_yz.attrs['units']='J/m^2'   
        ohc_yz.attrs.update(attrs)               
        ohc_yz.ohc.to_netcdf(save_dir+'ohc_yz_global_'+exp+'_'+run+'_'+str(startdate)+'_'+str(enddate)+'.nc')
    
        for basin in basin_masks.keys():
                print(basin)
                data_masked = data_weighted.where(basin_masks[basin])
                ohc_yz=data_masked.sum(dim='x')*rho_0*c_p
                ohc_yz.coords['lat']=data.nav_lat.mean(dim='x')
                ohc_yz['exp']=exp
                ohc_yz['run']=run
                ohc_yz.name='ohc'
                ohc_yz.attrs.update(attrs)    
                ohc_yz.attrs['long_name']='Zonally Integrated Heat content'
                ohc_yz.attrs['units']='J/m^2'                     
                ohc_yz.ohc.to_netcdf(save_dir+'ohc_yz_'+basin+'_'+exp+'_'+run+'_'+str(startdate)+'_'+str(enddate)+'.nc')
                    
print('all done')