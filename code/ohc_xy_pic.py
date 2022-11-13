#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ohc_xy_pic.py

Code to load ocean temperature data from the PIC of HadEM3-GC31-LL
calculate time series of depth-integrated Ocean Heat Content for de-drifting
SMURPHS ensemble. 

Required to reproduce data for Boland et al. 2022 (preprint https://www.essoar.org/doi/10.1002/essoar.10511062.3)
See https://github.com/emmomp/SMURPHS_OHC for details

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) 
and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

Updated Nov 2022

@author: emmomp@bas.ac.uk Emma J D Boland
"""
from datetime import date
import xarray as xr

rho_0 = 1.027e3 
c_p = 3850

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'HadGem3-GC31-LL PIC OHC data from Boland et al 2022 (https://www.essoar.org/doi/10.1002/essoar.10511062.3)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of HadGem3-GC31-LL CMIP6 PIC data, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995'}

save_dir = '../data_in/' #Directory to save data to
data_dir = '/badc/cmip6/data/CMIP6/CMIP/MOHC/HadGEM3-GC31-LL/piControl/r1i1p1f1/Omon' # Holding accessible via JASMIN

thkcello=xr.open_dataset('{}/thkcello/gn/v20190628/thkcello_Omon_HadGEM3-GC31-LL_piControl_r1i1p1f1_gn_185001-189912.nc'.format(data_dir))
dz=thkcello['thkcello'][0]

files = glob.glob('{}/thetao/gn/v20190628/thetao_*.nc'.format(data_dir))

with xr.open_mfdataset(files,concat_dim='time',combine='nested') as data:

    print('Got data, calculating ohc ')
    data_weighted = data.thetao*dz.squeeze(drop=True)
    ohc_xy=data_weighted.sum(dim='lev')*rho_0*c_p
    ohc_xy.name='ohc'
    ohc_xy.attrs['long_name']='Ocean Heat Content, full depth integrated'
    ohc_xy.attrs['units']='J/m^2'     
    ohc_xy.attrs.update(attrs)

    print('writing to file')
    ohc_xy.to_netcdf(save_dir+'pic_data/ohc_xy_pic.nc')
    
print('All done')
