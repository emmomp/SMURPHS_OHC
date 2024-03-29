#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ohc_by_basin_depth.py

Code to load ocean temperature data from the SMURPHS ensemble model output and 
calculate Ocean Heat Content integrals by basin and depth range. 

Required to reproduce data for Boland et al. 2023 (https://doi.org/10.1029/2022JC018725)
See https://github.com/emmomp/SMURPHS_OHC for details

thetao files can be found as part of the CEDA archive here: https://catalogue.ceda.ac.uk/uuid/5808b237bdb5485d9bc3595f39ce85e3
volcello can be accessed as part of the badc CMIP6 holdings, see https://help.ceda.ac.uk/article/4801-cmip6-data for details

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) 
and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

Updated Jun 2023

@author: emmomp@bas.ac.uk
"""

import xarray as xr
import glob
import dask
from datetime import date

rho_0 = 1.027e3 # Standard density
c_p = 3850 # Heat capacity of sea water

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'SMURPHS OHC data from Boland et al. 2023 (https://doi.org/10.1029/2022JC018725)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of the SMURPHS ensemble, See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806'}

save_dir = '../data_in/' #Directory to save data to

data_dir = '/gws/nopw/j04/smurphs/E/adittus/MASS_download/ocean/data_andrea/' # Path to thetao files

masks=xr.open_dataset(save_dir+'other_model_data/subbasins_eORCA1-GO6-Daley.nc')
# Match up grid formats by removing NEMO halo
masks=masks.isel(x=slice(1,-1),y=slice(1,-1)) 
masks=masks.rename({'x':'i','y':'j'})
volcello=xr.open_dataset('/badc/cmip6/data/CMIP6/CMIP/MOHC/HadGEM3-GC31-LL/historical/r1i1p1f3/Omon/volcello/gn/v20190624/volcello_Omon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_185001-189912.nc')
vol=volcello['volcello'][0]

depthlabels=['0-300m','300-700m','700-2000m','2km+']
depthbins = [0,300,700,2000,6001]
nd = len(depthbins)

basin_masks={
        'global':xr.full_like(vol,True),
        'so':(vol.latitude<=-35.0),
        'atl':masks.tmaskatl.astype(bool),
        'ind':masks.tmaskind.astype(bool),
        'pac':masks.tmaskpac.astype(bool)
            }

basin_longname={
        'global':'Global',
        'so':'Southern Ocean',
        'atl':'Atlantic',
        'ind':'Indian',
        'pac':'Pacific',
        }

exps=['hist-0p2','hist-0p4','hist-0p7','hist-1p0','hist-1p5']
runs=['r1i1p1f2','r2i1p1f2','r3i1p1f2','r4i1p1f2','r5i1p1f2']

for exp in exps:
    for run in runs:        
        # Check if output already exists
        files = glob.glob(save_dir+'ohc_tseries/ohc_*'+exp+'_'+run+'*nc')
        if len(files)==10:
            print('Skipping '+exp+' '+run)
        else:        
            print('Opening '+exp+' '+run)
        
            files = glob.glob('{}/{}/{}/Omon/thetao/gn/v20190213/thetao_*.nc'.format(data_dir,exp,run))
            #Initialise lists to hold data
            ohc_datasets={}
            ohc_bydepth_datasets={}   
            print('Found {} files'.format(len(files)))

            with xr.open_mfdataset(files,concat_dim='time',combine='nested') as data:
                print('Got data, calculating ohc')
                data_weighted = data['thetao']*vol

                for basin in basin_masks.keys():
                    print(basin)
                    with dask.config.set(**{'array.slicing.split_large_chunks': True}):
                        data_masked = data_weighted.where(basin_masks[basin])
                        vol_masked=vol.where(basin_masks[basin])
                        ohc_datasets[basin]=data_masked.sum(dim=['i','j','lev'])/vol_masked.sum(dim=['i','j','lev'])*rho_0*c_p

                        data_masked_binned = data_masked.groupby_bins(data.lev,depthbins,labels=depthlabels)
                        vol_masked_binned=vol_masked.groupby_bins(data.lev,depthbins,labels=depthlabels)
                        ohc_bydepth_datasets[basin]=data_masked_binned.sum(dim=['i','j','lev'])/vol_masked_binned.sum(dim=['i','j','lev'])*rho_0*c_p            

            print('Loading and writing to file')
            for basin in basin_masks.keys():
                print(basin)
                ohc=ohc_datasets[basin]
                ohc.load()
                ohc.name='ohc'
                ohc.attrs['long_name']=basin_longname[basin]+' depth integrated ocean heat content, volume averaged'
                ohc.attrs['units']='J/m**3'
                ohc['basin']=basin
                ohc['exp']=exp
                ohc['run']=run
                ohc.attrs.update(attrs)
                ohc.to_netcdf(save_dir+'ohc_tseries/ohc_'+exp+'_'+run+'_'+basin+'.nc')

                ohc_bydepth=ohc_bydepth_datasets[basin]
                ohc_bydepth.load()         
                ohc_bydepth.name='ohc'
                ohc_bydepth.attrs['long_name']=basin_longname[basin]+' heat content by depth bin, volume averaged'
                ohc_bydepth.attrs['units']='J/m**3'   
                ohc_bydepth['basin']=basin
                ohc_bydepth['exp']=exp
                ohc_bydepth['run']=run      
                ohc_bydepth.attrs.update(attrs)       
                ohc_bydepth.to_netcdf(save_dir+'ohc_tseries/ohc_bydepth_'+exp+'_'+run+'_'+basin+'.nc')                
            
            print('Done '+exp+' '+run)
         
print('All Done')
