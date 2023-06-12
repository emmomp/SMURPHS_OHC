#!/usr/bin/env python
# coding: utf-8
"""
calculate_AMOC.py

Code to load atlantic streamfunction data from the SMURPHS ensemble model output and calculate time series AMOC at 45N.

Required to reproduce data for Boland et al. 2023 (https://doi.org/10.1029/2022JC018725)
See https://github.com/emmomp/SMURPHS_OHC for details

zomsfatl files can be found as part of the CEDA archive here: https://catalogue.ceda.ac.uk/uuid/5808b237bdb5485d9bc3595f39ce85e3

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

Updated Jun 2023

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import xarray as xr
import numpy as np
from datetime import date
import glob

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'SMURPHS AMOC data from Boland et al. 2023 (https://doi.org/10.1029/2022JC018725)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of the SMURPHS ensemble, See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806'}

save_dir = '../data_in/' #Directory to save data to

data_dir = '/gws/nopw/j04/smurphs/E/adittus/MASS_download/ocean/data_andrea/' # Directory for zomsfatl data

exps=['hist-0p2','hist-0p4','hist-0p7','hist-1p0','hist-1p5']
runs=['r1i1p1f2','r2i1p1f2','r3i1p1f2','r4i1p1f2','r5i1p1f2']

for ie,exp in enumerate(exps):
    for run in runs:
        print('loading streamfunctions {} {}'.format(exp,run))
        files = glob.glob('{}/{}/{}/Omon/zomsfatl/gn/v20190213/thetao_*.nc'.format(data_dir,exp,run))
        with xr.open_mfdataset(files,concat_dim='time',combine='nested') as strf:
            print('loaded and combined')
            # Find 45N
            i45=np.abs(strf.nav_lat-45).argmin(dim='y')
            strf_45N=strf.zomsfatl.squeeze().isel(y=i45.load())
            amoc=strf_45N.max(dim='depthw').squeeze()
            amoc.attrs['units']='Sv'
            amoc.attrs['long_name']='AMOC strength at 45N'
            amoc.attrs['description']='AMOC strength at 45N, determined from depth maximum of zomsfatl, zonally averaged meridional overturning streamfn in the Atlantic'
            amoc_depth=strf_45N.idxmax(dim='depthw').squeeze()
            amoc_depth.attrs['units']='m'
            amoc_depth.attrs['long_name']='Depth of AMOC max at 45N'        
            amoc_all=xr.Dataset(data_vars={'AMOC_strength':amoc,'AMOC_depth':amoc_depth})
            amoc_all.attrs.update(attrs)   
            amoc_all['exp']=(('exp',),[new_exp])
            amoc_all['run']=(('run',),[run])             

            print('writing amoc to file')
            amoc_all.to_netcdf(save_dir+'amoc_tseries/amoc_45N_{}_{}.nc'.format(new_exp,run))

print('finished')





