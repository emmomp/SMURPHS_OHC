#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ohc_xy.py

Code to load ocean temperature data from the SMURPHS ensemble model output and 
calculate time series of depth-integrated Ocean Heat Content from 1955-2014. 

Required to reproduce data for Boland et al. 2022 (preprint https://www.essoar.org/doi/10.1002/essoar.10511062.3)
See https://github.com/emmomp/SMURPHS_OHC for details

For the SMURPHS ensemble, see Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) 
and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

Updated Nov 2022

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import utils
import baspy as bp
import xarray as xr
from datetime import date

rho_0 = 1.027e3 
c_p = 3850

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'SMURPHS OHC data from Boland et al 2022 (https://www.essoar.org/doi/10.1002/essoar.10511062.3)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of the SMURPHS ensemble, See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806'}

save_dir = '../data_in/' #Directory to save data to

griddata = xr.open_dataset(save_dir+'other_model_data/nemo_grid-T.nc')

depthlabels=['0-700m','0-2000m']
depthbins = [[0,701],[0,2001]]
nd = len(depthbins)

# Generate dataframe with all temp data
df = bp.catalogue(dataset='smurphs',Var=['thetao'])
nt = len(bp.get_files(df.iloc[0]))

exps = set(df.Experiment)
runs = set(df.RunID)

fname_struct='Var_Model_Experiment_RunID_Frequency_StartDate-EndDate_grid'
startdate=19550101
enddate=20141231

def preproc_data(ds):
    ds=ds.swap_dims({'time_counter':'time_centered'})
    return ds.drop('deptht_bounds')

for index, row in df.iterrows():
    exp = row[0]
    run = row[1]
    print('Opening row '+exp+' '+run)
    files = bp.get_files(row)
    # Pull out 1955 onwards
    sfiles = utils.get_date_range_files(fname_struct,files,startdate,enddate)
    print('Found ',len(sfiles),' files ')

    with xr.open_mfdataset(sfiles,preprocess=preproc_data,concat_dim='time_centered',combine='nested',data_vars='minimal', coords='minimal', compat='override') as data:
        print('Got data, calculating ohc ')
        data_weighted = data.thetao*griddata.thkcello.squeeze(drop=True)
        ohc_xy=data_weighted.sum(dim='deptht')*rho_0*c_p
        ohc_xy['exp']=exp
        ohc_xy['run']=run
        ohc_xy.name='ohc'
        ohc_xy.attrs['long_name']='Ocean Heat Content, full depth integrated'
        ohc_xy.attrs['units']='J/m^2'   
        ohc_xy.attrs.update(attrs)   
        
        ohc_xy_bybins=[]
        for s in range(0,2):
            foo=data_weighted.sel(deptht=slice(depthbins[s][0],depthbins[s][1])).sum(dim='deptht')*rho_0*c_p
            foo['depth_range']=depthlabels[s]
            ohc_xy_bybins.append(foo)
        ohc_xy_bybins=xr.concat(ohc_xy_bybins,'depth_range')  
        ohc_xy_bybins['exp']=exp
        ohc_xy_bybins['run']=run
        ohc_xy_bybins.name='ohc'
        ohc_xy_bybins.attrs['long_name']='Heat content, depth integrated'
        ohc_xy_bybins.attrs['units']='J/m^2'   
        ohc_xy_bybins.attrs.update(attrs)     
   
        ohc_xy.to_netcdf(save_dir+'ohc_xy/ohc_xy_'+exp+'_'+run+'_'+str(startdate)+'_'+str(enddate)+'.nc')
        ohc_xy_bybins.to_netcdf(save_dir+'ohc_xy/ohc_xy_bydepth_'+exp+'_'+run+'_'+str(startdate)+'_'+str(enddate)+'.nc')

print('All done')
    
