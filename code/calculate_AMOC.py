#!/usr/bin/env python
# coding: utf-8
"""
calculate_AMOC.py

Code to load atlantic streamfunction data from the SMURPHS ensemble model output and calculate time series AMOC at 45N.

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

See https://github.com/emmomp/SMURPHS_OHC for details

Created in Jan 2021

@author: emmomp
"""
import xarray as xr
import numpy as np
from datetime import date

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'SMURPHS AMOC data from Boland et al (in prep)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of the SMURPHS ensemble, See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806'}

def pre_proc(ds):
    return ds.swap_dims({'time_counter':'time_centered'})

save_dir = '../data_in/' #Directory to save data to

alt_exps=['hist-0p2','hist-0p4','hist-0p7','hist-1p0','hist-1p5']
exps=['historical0p2','historical0p4','historical0p7','historical1p0','historical1p5']
exp_map=dict(zip(alt_exps,exps))
runs=['r1i1p1f2','r2i1p1f2','r3i1p1f2','r4i1p1f2','r5i1p1f2']

for ie,exp in enumerate(alt_exps):
    for run in runs:
        print('loading streamfunctions {} {}'.format(exp,run))
        strf=xr.open_mfdataset('/gws/nopw/j04/smurphs/E/adittus/MASS_download/ocean/data_andrea/{}/{}/Omon/zomsfatl/gn/v20190213/zomsfatl_nemo_{}_{}_1y_*.nc'.format(exp,run,exp,run),
                      preprocess=pre_proc,combine='by_coords')
        
        new_exp=exp_map[exp]
        strf['exp']=(('exp',),[new_exp])
        strf['run']=(('run',),[run])
        
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

        print('writing amoc to file')
        amoc_all.to_netcdf(save_dir+'other_model_data/amoc_45N_{}_{}.nc'.format(new_exp,run))

print('finished')





