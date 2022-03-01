#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ohc_yz_pic.py

Code to load ocean temperature data from the PIC of HadEM3-GC31-LL calculate time series of zonally-integrated Ocean Heat Content for de-drifting SMURPHS ensemble. 

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

See https://github.com/emmomp/SMURPHS_OHC for details

Created on Mon Nov 29 15:52:09 2021

@author: emmomp@bas.ac.uk Emma J D Boland
"""
from datetime import date
import baspy as bp
import xarray as xr

rho_0 = 1.027e3 
c_p = 3850

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'HadGem3-GC31-LL PIC OHC data from Boland et al (in prep)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of HadGem3-GC31-LL CMIP6 PIC data, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995'}

save_dir = '../data_in/' #Directory to save data to

griddata = xr.open_dataset(save_dir+'other_model_data/nemo_grid-T.nc')
dx=griddata.e1t
ny=dx['y'].size
nx=dx['x'].size
dx=dx.isel(y=slice(1,ny-1),x=slice(1,nx-1))

df = bp.catalogue(dataset='cmip6',Var=['thetao'],Model='HadGEM3-GC31-LL',Experiment='piControl')
data= bp.open_dataset(df)

masks=xr.open_dataset(save_dir+'other_model_data/subbasins_eORCA1-GO6-Daley.nc')
masks=masks.isel(y=slice(1,ny-1),x=slice(1,nx-1))

basin_masks={
        'so':(masks.nav_lat<=-35.0),
        'atl':masks.tmaskatl.astype(bool),
        'ind':masks.tmaskind.astype(bool),
        'pac':masks.tmaskpac.astype(bool)
            }


print('Got data, calculating ohc')
data_slice = data.thetao.isel(time=slice(0,6000))
data_slice=data_slice.rename({'j':'y','i':'x','latitude':'nav_lat','longitude':'nav_lon','lev':'deptht'})
data_weighted=data_slice*dx.squeeze(drop=True)

ohc_yz=data_weighted.sum(dim='x')*rho_0*c_p
ohc_yz['nav_lat']=data_weighted['nav_lat'].mean(dim='x')
ohc_yz.name='ohc'
ohc_yz.attrs['long_name']='Ocean Heat Content, zonally integrated'
ohc_yz.attrs.update(attrs)    
ohc_yz.attrs['units']='J/m^2'       
   
print('writing global to file')
ohc_yz.to_netcdf(save_dir+'pic_data/ohc_yz_global_pic.nc')

for basin in basin_masks.keys():
    print(basin)
    data_masked = data_weighted.where(basin_masks[basin])
    ohc_yz=data_masked.sum(dim='x')*rho_0*c_p
    ohc_yz['nav_lat']=data_weighted['nav_lat'].mean(dim='x')
    ohc_yz.name='ohc'
    ohc_yz.attrs['long_name']='Ocean Heat Content, zonally integrated'
    ohc_yz.attrs['units']='J/m^2'    
    ohc_yz.attrs.update(attrs)     
    ohc_yz.to_netcdf(save_dir+'pic_data/ohc_yz_'+basin+'_pic.nc')

print('all done')
