#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ohc_by_basin_depth.py

Code to load ocean temperature data from the SMURPHS ensemble model output and 
calculate Ocean Heat Content integrals by basin and depth range. 

For the SMURPHS ensemble, see (See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806) 
and for the CMIP6 historical ensemble, Andrews et al. (2020) https://doi.org/10.1029/2019MS001995.

See https://github.com/emmomp/SMURPHS_OHC for details

Created on 17 Mar 2020

@author: emmomp@bas.ac.uk
"""

import baspy as bp
import xarray as xr
import glob
import cftime
from datetime import date

rho_0 = 1.027e3 # Standard density
c_p = 3850 # Heat capacity of sea water

attrs={'contact':'emmomp@bas.ac.uk',
       'references':'SMURPHS OHC data from Boland et al (in prep)',
       'date':'Created on '+date.today().strftime("%d/%m/%Y"),
       'notes':'Data produced by analysis of the SMURPHS ensemble, See Dittus et al. 2020 https://doi.org/10.1029/2019GL085806'}

save_dir = '../data_in/' #Directory to load from and save data to

masks=xr.open_dataset(save_dir+'subbasins_eORCA1-GO6-Daley.nc')
griddata = xr.open_dataset(save_dir+'/nemo_at491o_1m_19991201-20000101_grid-T.nc')
vol = (griddata.area*griddata.thkcello).squeeze()
del vol['time_centered']

depthlabels=['0-300m','300-700m','700-2000m','2km+']
depthbins = [0,300,700,2000,6001]
nd = len(depthbins)

basin_masks={
        'global':xr.full_like(vol,True),
        'so':(vol.nav_lat<=-35.0),
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
        'na':'North Atlantic region (Latitude: [45°N, 67°N] Longitude: [100°W, 20°E] Depth: top 1056m)'
        }

# Generate dataframe with all temp data
df = bp.catalogue(dataset='smurphs',Var=['thetao'])
nt = len(bp.get_files(df.iloc[0]))

exps = set(df.Experiment)
runs = set(df.RunID)
                
def preproc_data(ds):
    return ds.drop('deptht_bounds')

for index, row in df.iterrows():
    exp = row[0]
    run = row[1]
    files = glob.glob(save_dir+'ohc_*'+exp+'_'+run+'*nc')
    # Check if output already exists
    if len(files)==10:
        print('Skipping '+exp+' '+run)
    else:        
        print('Opening row '+exp+' '+run)
        files = bp.get_files(row)
        print('Found ',len(files),' files')
        if len(files)==1980:
            # Split file list into slices of 50 to stop massive slowdown
            nf=50
            file_slices=[files[i:i + nf] for i in range(0, len(files), nf)]
            #Initialise lists to hold data
            ohc_datasets={}
            ohc_bydepth_datasets={}
            for basin in basin_masks.keys():
                ohc_datasets[basin]=[]
                ohc_bydepth_datasets[basin]=[]
            
            for a,afiles in enumerate(file_slices):
                print('Slice ',a)
                try: 
                    with xr.open_mfdataset(afiles,preprocess=preproc_data,concat_dim='time_counter') as data:
                        print('Got data, calculating ohc')
                        data_weighted = data*vol
                    
                        for basin in basin_masks.keys():
                            print(basin)
                            data_masked = data_weighted.where(basin_masks[basin])
                            ohc_datasets[basin].append(data_masked.thetao.sum(dim=['x','y','deptht'])*rho_0*c_p)
                            
                            data_masked_binned = data_masked.groupby_bins(data.deptht,depthbins,labels=depthlabels)
                            ohc_bydepth_datasets[basin].append(data_masked_binned.sum(dim=['x','y','deptht'])*rho_0*c_p)  

                except: #Some of the time data doesn't read nicely
                    print('Reading files individually')
                    for file in afiles:
                        with xr.open_dataset(file) as data:
                            data=data.drop('deptht_bounds')
                            if file==afiles[0]:
                                deptht=data['deptht']
                            if data['deptht'][0]>10:
                                data['deptht']=deptht
                                lastdate=ohc_datasets[basin][-1].time_centered.compute().item()
                                data['time_centered']=cftime.Datetime360Day(lastdate.year,lastdate.month+1,lastdate.day)
                            print(file)
                            print('Got data, calculating ohc')
                            data_weighted = data*vol
                    
                            for basin in basin_masks.keys():
                                print(basin)
                                data_masked = data_weighted.where(basin_masks[basin])
                                ohc_datasets[basin].append(data_masked.thetao.sum(dim=['x','y','deptht'])*rho_0*c_p)
                                
                                data_masked_binned = data_masked.groupby_bins(data.deptht,depthbins,labels=depthlabels)
                                ohc_bydepth_datasets[basin].append(data_masked_binned.sum(dim=['x','y','deptht'])*rho_0*c_p)  

                        
            print('Combining, loading, and writing to file')
            for basin in basin_masks.keys():
                print(basin)
                ohc=xr.concat(ohc_datasets[basin],'time_counter')
                ohc.load()
                ohc.name='ohc'
                ohc.attrs['long_name']=basin_longname[basin]+' depth integrated ocean heat content'
                ohc.attrs['units']='J'
                ohc['basin']=basin
                ohc['exp']=exp
                ohc['run']=run
                ohc.attrs.update(attrs)
                ohc.to_netcdf(save_dir+'ohc_'+exp+'_'+run+'_'+basin+'.nc')
                
                ohc_bydepth=xr.concat(ohc_bydepth_datasets[basin],'time_counter')
                ohc_bydepth.load()         
                ohc_bydepth.name='ohc'
                ohc_bydepth.attrs['long_name']=basin_longname[basin]+' heat content by depth bin'
                ohc_bydepth.attrs['units']='J'   
                ohc_bydepth['basin']=basin
                ohc_bydepth['exp']=exp
                ohc_bydepth['run']=run      
                ohc_bydepth.attrs.update(attrs)       
                ohc_bydepth.to_netcdf(save_dir+'ohc_bydepth_'+exp+'_'+run+'_'+basin+'.nc')                
            
            print('Done '+exp+' '+run)
        else: 
            print('Skipping '+exp+' '+run)
         
print('All Done')
