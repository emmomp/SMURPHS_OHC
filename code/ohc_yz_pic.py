#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ohc_yz_pic.py

Code to load ocean temperature data from the PIC of HadEM3-GC31-LL calculate time series of zonally-integrated Ocean Heat Content for de-drifting SMURPHS ensemble. 

Required to reproduce data for Boland et al. 2022 (preprint https://www.essoar.org/doi/10.1002/essoar.10511062.3)
See https://github.com/emmomp/SMURPHS_OHC for details

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

Created on Mon Nov 29 15:52:09 2021

@author: emmomp@bas.ac.uk Emma J D Boland
"""
from datetime import date
import xarray as xr
import glob

rho_0 = 1.027e3 
c_p = 3850

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'HadGem3-GC31-LL PIC OHC data from Boland et al 2022 (https://www.essoar.org/doi/10.1002/essoar.10511062.3)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of HadGem3-GC31-LL CMIP6 PIC data, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995'}

save_dir = '../data_in/' #Directory to save data to
data_dir = '/badc/cmip6/data/CMIP6/CMIP/MOHC/HadGEM3-GC31-LL/piControl/r1i1p1f1/Omon' # Holding accessible via JASMIN

griddata = xr.open_dataset(save_dir+'other_model_data/nemo_grid-T.nc')
dx=griddata.e1t
# Match up grid formats by removing NEMO halo and renaming coords
dx=dx.isel(x=slice(1,-1),y=slice(1,-1)) 
dx=dx.rename({'x':'i','y':'j','nav_lon':'longitude','nav_lat':'latitude'})

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

files = glob.glob('{}/thetao/gn/v20190628/thetao_*.nc'.format(data_dir))

with xr.open_mfdataset(files,concat_dim='time',combine='nested') as data:
    
    print('global')
    data_weighted=data.thetao*dx
    ohc_yz=data_weighted.sum(dim='i')*rho_0*c_p
    ohc_yz['latitude']=data_weighted['latitude'].mean(dim='i')
    ohc_yz.name='ohc'
    ohc_yz.attrs['long_name']='Ocean Heat Content, zonally integrated'
    ohc_yz.attrs['units']='J/m^2'    
    ohc_yz.attrs.update(attrs)     
    ohc_yz.to_netcdf(save_dir+'pic_data/ohc_yz_global_pic.nc')
    
    for basin in basin_masks.keys():
        print(basin)
        data_masked = data_weighted.where(basin_masks[basin])
        ohc_yz=data_masked.sum(dim='i')*rho_0*c_p
        ohc_yz['nav_lat']=data_weighted['latitude'].mean(dim='i')
        ohc_yz.name='ohc'
        ohc_yz.attrs['long_name']='Ocean Heat Content, zonally integrated'
        ohc_yz.attrs['units']='J/m^2'    
        ohc_yz.attrs.update(attrs)     
        ohc_yz.to_netcdf(save_dir+'pic_data/ohc_yz_'+basin+'_pic.nc')

print('all done')
