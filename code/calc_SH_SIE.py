#!/usr/bin/env python
# coding: utf-8
"""
calculate_SH_SIE.py

Code to calculate Southern Hemisphere Sea Ice Extent from the SMURPHS ensemble model output.

siconca files can be found directly from https://gws-access.ceda.ac.uk/public/smurphs/SMURPHS/ or as part of the CEDA archive here: https://catalogue.ceda.ac.uk/uuid/5808b237bdb5485d9bc3595f39ce85e3

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

See https://github.com/emmomp/SMURPHS_OHC for details

Created in Jan 2021

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import xarray as xr
from datetime import date

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'SMURPHS SIE data from Boland et al (in prep)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of the SMURPHS ensemble, See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806'}

save_dir = '../data_in/' #Directory to save data to

si_dir = '/gws/nopw/j04/smurphs/public/SMURPHS/coupled/siconca/' # Directory to find sea ice data

alt_exps=['hist-0p2','hist-0p4','hist-0p7','hist-1p0','hist-1p5']
exps=['historical0p2','historical0p4','historical0p7','historical1p0','historical1p5']
exp_map=dict(zip(alt_exps,exps))
runs=['r1i1p1f2','r2i1p1f2','r3i1p1f2','r4i1p1f2','r5i1p1f2']

areacella=xr.open_dataset(save_dir+'other_model_data/areacella.nc')

print('Loading SI conc')
si=[[],[],[],[],[]]
for ie,exp in enumerate(alt_exps):    
    for run in runs:
        ds=xr.open_dataset(si_dir+'siconca_Amon_HadGEM3-GC3.1_{}_{}_gn_185001-201412.nc'.format(exp,run))
        si[ie].append(ds.siconca)
si_all=xr.combine_nested(si,concat_dim=['exp','run'],coords='minimal')
si_all['run']=(('run'),runs)
si_all['exp']=(('exp'),exps)

print('Calculating SI extent')
si_so = si_all.where(si_all.latitude<-40,drop=True)
siarea=areacella.cell_area.where(si_so>0.15)/1e6
siarea=siarea.sum(dim=['latitude','longitude'])
siarea.attrs={'long_name':'Sea Ice Area (Southern Hemisphere)','units':'km^2'}
siarea.name='siarea_sh'
siarea.attrs.update(attrs)
siarea.to_netcdf(save_dir+'SIE_SH.nc')
print('All done')
