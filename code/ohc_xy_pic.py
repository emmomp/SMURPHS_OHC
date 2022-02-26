#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ohc_xy_pic.py

Code to load ocean temperature data from the PIC of HadEM3-GC31-LL
calculate time series of depth-integrated Ocean Heat Content for de-drifting
SMURPHS ensemble. 

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) 
and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

See https://github.com/emmomp/SMURPHS_OHC for details

Created on Wed Nov 24 12:27:40 2021

@author: emmomp
"""
from datetime import date
import baspy as bp

rho_0 = 1.027e3 
c_p = 3850

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'HadGem3-GC31-LL PIC OHC data from Boland et al (in prep)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of HadGem3-GC31-LL CMIP6 PIC data, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995'}

save_dir = '../data_in/' #Directory to save data to

df = bp.catalogue(dataset='cmip6',Var=['thkcello'],Model='HadGEM3-GC31-LL',Experiment='piControl')
data=bp.open_dataset(df)
dz=data.thkcello

df = bp.catalogue(dataset='cmip6',Var=['thetao'],Model='HadGEM3-GC31-LL',Experiment='piControl')
data= bp.open_dataset(df)

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
