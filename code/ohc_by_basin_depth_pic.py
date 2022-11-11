#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ohc_by_basin_depth_pic.py

Code to load ocean temperature data from the PIC of HadEM3-GC31-LL
calculate Ocean Heat Content integrals by basin and depth range for de-drifting
SMURPHS ensemble. 

Required to reproduce data for Boland et al. 2022 (preprint https://www.essoar.org/doi/10.1002/essoar.10511062.3)
See https://github.com/emmomp/SMURPHS_OHC for details

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) 
and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

Updated Nov 2022

@author: emmomp@bas.ac.uk Emma J D Boland
"""

import baspy as bp
import xarray as xr
import glob
from datetime import date

rho_0 = 1.027e3 
c_p = 3850

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'HadGem3-GC31-LL PIC OHC data from Boland et al 2022 (https://www.essoar.org/doi/10.1002/essoar.10511062.3)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of HadGem3-GC31-LL CMIP6 PIC data, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995'}

save_dir = '../data_in/' #Directory to save data to

data_dir = '/badc/cmip6/data/CMIP6/CMIP/MOHC/HadGEM3-GC31-LL/piControl/r1i1p1f1/Omon' # Holding accessible via JASMIN

masks=xr.open_dataset(save_dir+'other_model_data/subbasins_eORCA1-GO6-Daley.nc')
# Match up grid formats by removing NEMO halo
masks=masks.isel(x=slice(1,-1),y=slice(1,-1)) 
masks=masks.rename({'x':'i','y':'j'})
volcello=xr.open_dataset(data_dir+'/volcello/gn/v20190628/volcello_Omon_HadGEM3-GC31-LL_piControl_r1i1p1f1_gn_185001-189912.nc')
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

files = glob.glob('{}/thetao/gn/v20190628/thetao_*.nc'.format(data_dir))


with xr.open_mfdataset(files,concat_dim='time',combine='nested') as data:
    print('Got data, calculating ohc')
    data_weighted = data['thetao']*vol

    for basin in basin_masks.keys():        
    # Check if output already exists
    files = glob.glob(save_dir+'ohc_tseries/pic_data/ohc_pic_'+basin+'.nc')
    if len(files)==2:
        print('Skipping {}'.format(basin))
    else:    
        print(basin)
        data_masked = data_weighted.where(basin_masks[basin])
        vol_masked = vol.where(basin_masks[basin])
        ohc=data_masked.sum(dim=['i','j','lev'])/vol_masked.sum(dim=['i','j','lev'])*rho_0*c_p

        data_masked_binned = data_masked.groupby_bins(data.lev,depthbins,labels=depthlabels)
        vol_masked_binned = vol_masked.groupby_bins(data.lev,depthbins,labels=depthlabels)
        ohc_bydepth=data_masked_binned.sum(dim=['i','j','lev'])/vol_masked_binned.sum(dim=['i','j','lev'])*rho_0*c_p     

        ohc.name='ohc'
        ohc.attrs['long_name']=basin_longname[basin]+' depth integrated ocean heat content, volume averaged'
        ohc.attrs['units']='J/m**3'
        ohc['basin']=basin
        ohc.attrs.update(attrs)     
        ohc.to_netcdf(save_dir+'pic_data/ohc_pic_'+basin+'.nc')

        ohc_bydepth.name='ohc'
        ohc_bydepth.attrs['long_name']=basin_longname[basin]+' heat content by depth bin, volume averaged'
        ohc_bydepth.attrs['units']='J/m**3'   
        ohc_bydepth['basin']=basin
        ohc_bydepth.attrs.update(attrs)     
        ohc_bydepth.to_netcdf(save_dir+'pic_data/ohc_pic_bydepth_'+basin+'.nc')

print('All Done')
