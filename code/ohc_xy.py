#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ohc_xy.py

Code to load ocean temperature data from the SMURPHS ensemble model output and 
calculate time series of depth-integrated Ocean Heat Content from 1955-2014. 

Required to reproduce data for Boland et al. 2022 (preprint https://www.essoar.org/doi/10.1002/essoar.10511062.3)
See https://github.com/emmomp/SMURPHS_OHC for details

For the SMURPHS ensemble, see Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) 
and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

Updated Nov 2022

@author: emmomp@bas.ac.uk Emma J D Boland
"""
import glob
import xarray as xr
from datetime import date

rho_0 = 1.027e3 
c_p = 3850

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'SMURPHS OHC data from Boland et al 2022 (https://www.essoar.org/doi/10.1002/essoar.10511062.3)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of the SMURPHS ensemble, See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806'}

save_dir = '../data_in/' #Directory to save data to

data_dir = '/gws/nopw/j04/smurphs/E/adittus/MASS_download/ocean/data_andrea/' # Holding accessible via JASMIN

thkcello=xr.open_dataset('/badc/cmip6/data/CMIP6/CMIP/MOHC/HadGEM3-GC31-LL/historical/r1i1p1f3/Omon/thkcello/gn/v20190624/thkcello_Omon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_185001-189912.nc')
dz=thkcello['thkcello'][0]

depthlabels=['0-700m','0-2000m']
depthbins = [[0,701],[0,2001]]
nd = len(depthbins)

exps=['hist-0p2','hist-0p4','hist-0p7','hist-1p0','hist-1p5']
runs=['r1i1p1f2',' r2i1p1f2','r3i1p1f2','r4i1p1f2','r5i1p1f2']

startdate='1955-01-01'
enddate='2015-01-01'

for exp in exps:
    for run in runs:               
        files = glob.glob('{}/{}/{}/Omon/thetao/gn/v20190213/thetao_*195001*.nc'.format(data_dir,exp,run))+ \
                glob.glob('{}/{}/{}/Omon/thetao/gn/v20190213/thetao_*200001*.nc'.format(data_dir,exp,run)) 

        with xr.open_mfdataset(files,concat_dim='time',combine='nested',data_vars='minimal', coords='minimal', compat='override') as data:
            print('Got data, calculating ohc ')
            data=data.sel(time=slice(startdate,enddate))
            data_weighted = data.thetao*dz
            ohc_xy=data_weighted.sum(dim='lev')*rho_0*c_p
            ohc_xy['exp']=exp
            ohc_xy['run']=run
            ohc_xy.name='ohc'
            ohc_xy.attrs['long_name']='Ocean Heat Content, full depth integrated'
            ohc_xy.attrs['units']='J/m^2'   
            ohc_xy.attrs.update(attrs)   

            ohc_xy_bybins=[]
            for s in range(0,2):
                foo=data_weighted.sel(lev=slice(depthbins[s][0],depthbins[s][1])).sum(dim='lev')*rho_0*c_p
                foo['depth_range']=depthlabels[s]
                ohc_xy_bybins.append(foo)
            ohc_xy_bybins=xr.concat(ohc_xy_bybins,'depth_range')  
            ohc_xy_bybins['exp']=exp
            ohc_xy_bybins['run']=run
            ohc_xy_bybins.name='ohc'
            ohc_xy_bybins.attrs['long_name']='Heat content, depth integrated'
            ohc_xy_bybins.attrs['units']='J/m^2'   
            ohc_xy_bybins.attrs.update(attrs)     

            ohc_xy.to_netcdf(save_dir+'ohc_xy/ohc_xy_'+exp+'_'+run+'_'+str(startdate)+'_'+str(enddate)+'.nc')
            ohc_xy_bybins.to_netcdf(save_dir+'ohc_xy/ohc_xy_bydepth_'+exp+'_'+run+'_'+str(startdate)+'_'+str(enddate)+'.nc')

print('All done')
    
