#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calculate_vol_scaling.py

Code to calculate the volume of the HegGEM3-GC3.1-LL model by basin and depth bin.

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

See https://github.com/emmomp/SMURPHS_OHC for details

Created in Jun 2022

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import xarray as xr
from datetime import date

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'SMURPHS ocean volumes from Boland et al (in prep)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of the SMURPHS ensemble, See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806'}

data_dir = '../data_in/'
save_dir = '../data_in/' #Directory to save data to

basins=['global','atl','pac','so','ind']
depthlabels=['0-300m','300-700m','700-2000m','2km+']
depthbins = [0,300,700,2000,6001]
masks=xr.open_dataset(data_dir+'other_model_data/subbasins_eORCA1-GO6-Daley.nc')
griddata = xr.open_dataset(data_dir+'other_model_data/nemo_grid-T.nc')
vol=(griddata.area*griddata.thkcello).squeeze()
basin_masks={
        'global':xr.full_like(vol,True),
        'so':(vol.nav_lat<=-35.0),
        'atl':masks.tmaskatl.astype(bool),
        'ind':masks.tmaskind.astype(bool),
        'pac':masks.tmaskpac.astype(bool)
            }

# Calculate volume by basin and depth range for scaling
basin_vol=[]
for basin in basins:
    vol_fulld=vol.where(basin_masks[basin]).sum()
    vol_fulld['basin']=basin
    vol_fulld['deptht_bins']='Full Depth'
    vol_binned = vol.where(basin_masks[basin]).groupby_bins(vol.deptht,depthbins,labels=depthlabels).sum(['y','x','deptht'])    
    basin_vol.append(xr.concat([vol_fulld,vol_binned],'deptht_bins'))
basin_vol=xr.concat(basin_vol,'basin')
basin_vol.name='basin_volume'
basin_vol.attrs['units']='m^3'
basin_vol.attrs['long_name']='Ocean Volume'
basin_vol.attrs['description']='Summed volume of ocean cells by depth range and basin'
surface_area=griddata['area'].where(griddata.thkcello.isel(deptht=0).squeeze()>0).sum(['y','x'])
surface_area.name='surface_area'
surface_area.attrs['units']='m^2'
surface_area.attrs['long_name']='Ocean Surface Area'

scalings = xr.merge([basin_vol,surface_area],combine_attrs="drop_conflicts")
scalings.attrs.update(attrs)
scalings.to_netcdf(save_dir+'other_model_data/vol_scalings.nc')

print('Finished')
