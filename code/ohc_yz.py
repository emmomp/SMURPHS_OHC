#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ohc_yz.py

Code to load ocean temperature data from the SMURPHS ensemble model output and calculate time series of zonally-integrated Ocean Heat Content from 1955-2014. 

Required to reproduce data for Boland et al. 2022 (preprint https://www.essoar.org/doi/10.1002/essoar.10511062.3)
See https://github.com/emmomp/SMURPHS_OHC for details

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

Updated Nov 2022

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import xarray as xr
from datetime import date
import glob

rho_0 = 1.027e3 
c_p = 3850

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'SMURPHS OHC data from Boland et al 2022 (https://www.essoar.org/doi/10.1002/essoar.10511062.3)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of the SMURPHS ensemble, See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806'}

save_dir = '../data_in/' #Directory to save data to

griddata = xr.open_dataset(save_dir+'other_model_data/nemo_grid-T.nc')
dx=griddata.e1t
# Match up grid formats by removing NEMO halo and renaming coords
dx=dx.isel(x=slice(1,-1),y=slice(1,-1)) 
dx=dx.rename({'x':'i','y':'j','nav_lon':'longitude','nav_lat':'latitude'})

exps=['hist-0p2','hist-0p4','hist-0p7','hist-1p0','hist-1p5']
runs=['r1i1p1f2','r2i1p1f2','r3i1p1f2','r4i1p1f2','r5i1p1f2']

masks=xr.open_dataset(save_dir+'other_model_data/subbasins_eORCA1-GO6-Daley.nc')
# Match up grid formats by removing NEMO halo
masks=masks.isel(x=slice(1,-1),y=slice(1,-1)) 
masks=masks.rename({'x':'i','y':'j'})

basin_masks={
        'so':(dx.latitude<=-35.0),
        'atl':masks.tmaskatl.astype(bool),
        'ind':masks.tmaskind.astype(bool),
        'pac':masks.tmaskpac.astype(bool)
            }

startdate='1955-01-01'
enddate='2015-01-01'

for exp in exps:
    for run in runs:   
        # Check if output already exists
        files = glob.glob(save_dir+'ohc_yz/ohc_yz_*'+exp+'_'+run+'*nc')
        if len(files)==2:
            print('Skipping '+exp+' '+run)
        else:        
            print('Opening '+exp+' '+run)
            
            files = glob.glob('{}/{}/{}/Omon/thetao/gn/v20190213/thetao_*195001*.nc'.format(data_dir,exp,run))+ \
                    glob.glob('{}/{}/{}/Omon/thetao/gn/v20190213/thetao_*200001*.nc'.format(data_dir,exp,run)) 
            
            with xr.open_mfdataset(files,concat_dim='time',combine='nested',data_vars='minimal', coords='minimal', compat='override') as data:
                print('Got data, calculating ohc ')
                data=data.sel(time=slice(startdate,enddate))
                data_weighted = data.thetao*dx
                ohc_yz=data_weighted.sum(dim='i')*rho_0*c_p
                ohc_yz.coords['lat']=data.latitude.mean(dim='x')
                ohc_yz['exp']=exp
                ohc_yz['run']=run
                ohc_yz.name='ohc'
                ohc_yz.attrs['long_name']='Ocean Heat Content, zonally integrated'
                ohc_yz.attrs['units']='J/m^2'   
                ohc_yz.attrs.update(attrs)               
                ohc_yz.to_netcdf(save_dir+'ohc_yz/ohc_yz_global_'+exp+'_'+run+'_'+str(startdate)+'_'+str(enddate)+'.nc')
    
                for basin in basin_masks.keys():
                    print(basin)
                    data_masked = data_weighted.where(basin_masks[basin])
                    ohc_yz=data_masked.sum(dim='i')*rho_0*c_p
                    ohc_yz.coords['lat']=data.latitude.mean(dim='x')
                    ohc_yz['exp']=exp
                    ohc_yz['run']=run
                    ohc_yz.name='ohc'
                    ohc_yz.attrs.update(attrs)    
                    ohc_yz.attrs['long_name']='Zonally Integrated Heat content'
                    ohc_yz.attrs['units']='J/m^2'                     
                    ohc_yz.to_netcdf(save_dir+'ohc_yz/ohc_yz_'+basin+'_'+exp+'_'+run+'_'+str(startdate)+'_'+str(enddate)+'.nc')
                    
print('all done')
