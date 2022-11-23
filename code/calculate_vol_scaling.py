#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calculate_vol_scaling.py

Code to calculate the volume of the HegGEM3-GC3.1-LL model by basin and depth bin, and the zonal width of each basin.

Required to reproduce data for Boland et al. 2022 (preprint https://www.essoar.org/doi/10.1002/essoar.10511062.3)
See https://github.com/emmomp/SMURPHS_OHC for details

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

Updated Nov 2022

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import xarray as xr
from datetime import date

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'SMURPHS ocean volumes from Boland et al 2022 (https://www.essoar.org/doi/10.1002/essoar.10511062.3)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of the SMURPHS ensemble, See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806'}

data_dir = '../data_in/'
save_dir = '../data_in/' #Directory to save data to

masks=xr.open_dataset(save_dir+'other_model_data/subbasins_eORCA1-GO6-Daley.nc')
# Match up grid formats by removing NEMO halo
masks=masks.isel(x=slice(1,-1),y=slice(1,-1)) 
masks=masks.rename({'x':'i','y':'j'})
volcello=xr.open_dataset('/badc/cmip6/data/CMIP6/CMIP/MOHC/HadGEM3-GC31-LL/historical/r1i1p1f3/Omon/volcello/gn/v20190624/volcello_Omon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_185001-189912.nc')
vol=volcello['volcello'][0]
areacello=xr.open_dataset('/badc/cmip6/data/CMIP6/CMIP/MOHC/HadGEM3-GC31-LL/piControl/r1i1p1f1/Ofx/areacello/gn/v20190709/areacello_Ofx_HadGEM3-GC31-LL_piControl_r1i1p1f1_gn.nc')
area=areacello['areacello']

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

griddata = xr.open_dataset(save_dir+'other_model_data/nemo_grid-T.nc')
dx=griddata.e1t
dx=dx.where(dx>0)
# Match up grid formats by removing NEMO halo and renaming coords
dx=dx.isel(x=slice(1,-1),y=slice(1,-1)) 
dx=dx.rename({'x':'i','y':'j','nav_lon':'longitude','nav_lat':'latitude'})

# Calculate volume by basin and depth range for scaling
basin_vol=[]
basin_width=[]
for basin in basin_masks.keys():
    vol_fulld=vol.where(basin_masks[basin]).sum()
    vol_fulld['basin']=basin
    vol_fulld['lev_bins']='Full Depth'
    vol_binned = vol.where(basin_masks[basin]).groupby_bins(vol.lev,depthbins,labels=depthlabels).sum(['i','j','lev'])    
    basin_vol.append(xr.concat([vol_fulld,vol_binned],'lev_bins'))
    widthx=dx.where(basin_masks[basin]).sum(dim='i')
    widthx['basin']=basin
    basin_width.append(widthx)    
basin_vol=xr.concat(basin_vol,'basin')
basin_vol.name='basin_volume'
basin_vol.attrs['units']='m^3'
basin_vol.attrs['long_name']='Ocean Volume'
basin_vol.attrs['description']='Summed volume of ocean cells by depth range and basin'
surface_area=area.where(vol.isel(lev=0).squeeze()>0).sum(['i','j'])
surface_area.name='surface_area'
surface_area.attrs['units']='m^2'
surface_area.attrs['long_name']='Ocean Surface Area'

scalings = xr.merge([basin_vol,surface_area],combine_attrs="drop_conflicts")
scalings.attrs.update(attrs)
scalings.to_netcdf(save_dir+'other_model_data/vol_scalings.nc')

basin_width=xr.concat(basin_width,'basin')
ny=basin_width['y'].size
basin_width=basin_width.rename({'z':'lev'})
basin_width.attrs['units']='m'
basin_width.attrs['long_name']='Basin Zonal Width'
basin_width.attrs['description']='Zonally summed x-width by depth, latitude and basin'
basin_width.attrs.update(attrs)
basin_width.to_netcdf(save_dir+'other_model_data/width_scalings.nc')

print('Finished')
