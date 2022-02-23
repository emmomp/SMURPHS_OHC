#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ohc_by_basin_depth_pic.py

Code to load ocean temperature data from the PIC of HadEM3-GC31-LL
calculate Ocean Heat Content integrals by basin and depth range for de-drifting
SMURPHS ensemble. 

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) 
and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

See https://github.com/emmomp/SMURPHS_OHC for details

Created on 17 Mar 2020

@author: emmomp@bas.ac.uk
"""

import baspy as bp
import xarray as xr
from datetime import date

rho_0 = 1.027e3 
c_p = 3850

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'HadGem3-GC31-LL PIC OHC data from Boland et al (in prep)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of HadGem3-GC31-LL CMIP6 PIC data, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995'}

save_dir = '../data_in/' #Directory to save data to

masks=xr.open_dataset(save_dir+'other_model_data/subbasins_eORCA1-GO6-Daley.nc')
# Match up grid formats
masks=masks.isel(x=slice(1,-1),y=slice(1,-1)) 
masks=masks.rename({'x':'i','y':'j'})

df = bp.catalogue(dataset='cmip6',Var=['volcello'],Model='HadGEM3-GC31-LL',Experiment='piControl')
data=bp.open_dataset(df)
vol=data.volcello

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

# Generate dataframe with all temp data
df = bp.catalogue(dataset='cmip6',Var=['thetao'],Model='HadGEM3-GC31-LL',Experiment='piControl')
nt = len(bp.get_files(df.iloc[0]))

print('Opening data')
ds = bp.open_dataset(df)
data=ds.thetao
print('Got data, calculating ohc')
data_weighted = data*vol

for basin in basin_masks.keys():
    print(basin)
    data_masked = data_weighted.where(basin_masks[basin])
    ohc=data_masked.sum(dim=['i','j','lev'])*rho_0*c_p
    
    data_masked_binned = data_masked.groupby_bins(data.lev,depthbins,labels=depthlabels)
    ohc_bydepth=data_masked_binned.sum(dim=['i','j','lev'])*rho_0*c_p     
   
    ohc.name='ohc'
    ohc.attrs['long_name']=basin_longname[basin]+' depth integrated ocean heat content'
    ohc.attrs['units']='J'
    ohc['basin']=basin
    ohc.attrs.update(attrs)     
    ohc.to_netcdf(save_dir+'pic_data/ohc_pic_'+basin+'.nc')
          
    ohc_bydepth.name='ohc'
    ohc_bydepth.attrs['long_name']=basin_longname[basin]+' heat content by depth bin'
    ohc_bydepth.attrs['units']='J'   
    ohc_bydepth['basin']=basin
    ohc_bydepth.attrs.update(attrs)     
    ohc_bydepth.to_netcdf(save_dir+'pic_data/ohc_pic_bydepth_'+basin+'.nc')

print('All Done')
