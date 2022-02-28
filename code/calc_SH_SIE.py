#!/usr/bin/env python
# coding: utf-8

# In[1]:


import xarray as xr
import matplotlib.pyplot as plt
import numpy as np


# In[2]:


exps=['historical0p2','historical0p4','historical0p7','historical1p0','historical1p5']
exp_names=['0.2','0.4','0.7','1.0','1.5']
runs=['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1']
smurphs_cmap= [(255/255, 0 , 0),(255/255 ,165/255,   0),(190/255, 190/255, 190/255),(0   ,0 ,255/255),(160 /255, 32/255 ,240/255)]


# In[3]:


exp=exps[0]
run=runs[0]
ds_test=xr.open_dataset('~/smurphs_ensemble/SMURPHS_coupled_ensemble/{}/{}/Amon/siconca/gn/v20190213/siconca_Amon_HadGEM3-GC3.1_{}_{}_gn_185001-201412.nc'.format(exp,run,exp,run))


# In[4]:


ds_test


# In[3]:


si=[[],[],[],[],[]]
for ie,exp in enumerate(exps):    
    for run in runs:
        ds=xr.open_dataset('~/smurphs_ensemble/SMURPHS_coupled_ensemble/{}/{}/Amon/siconca/gn/v20190213/siconca_Amon_HadGEM3-GC3.1_{}_{}_gn_185001-201412.nc'.format(exp,run,exp,run))
        si[ie].append(ds.siconca)
                


# In[4]:


si_all=xr.combine_nested(si,concat_dim=['exp','run'],coords='minimal')


# In[5]:


si_all['run']=(('run'),runs)
si_all['exp']=(('exp'),exps)
si_all


# In[6]:


si_so = si_all.where(si_all.latitude<-40,drop=True)


# In[7]:


areacella=xr.open_dataset('/gws/nopw/j04/smurphs/E/adittus/areacella.nc')


# In[8]:


areacella.cell_area


# In[9]:


siarea=areacella.cell_area.where(si_so>0.15)/1e6


# In[10]:


siarea.attrs={'long_name':'Sea Ice Area (Southern Hemisphere)','units':'km^2'}


# In[11]:


siarea.name='siarea_sh'


# In[12]:


siarea_sum=siarea.sum(dim=['latitude','longitude'])


# In[13]:


siarea_sum.to_netcdf('SIA_SH.nc')


# In[14]:


siarea_march=siarea_sum.where(siarea_sum.time.dt.month==3,drop=True)


# In[21]:


siarea_sep=siarea_sum.where(siarea_sum.time.dt.month==9,drop=True)


# In[17]:


tt=siarea_march.time.sel(time=slice('1900-01-01',None))
tt_index=tt.dt.year+(tt.dt.dayofyear-1)/365.25


# In[22]:


#dplot=siarea_march-siarea_march.sel(time=slice('1850-01-01','1900-01-01')).mean(dim=['time','run'])
dplot=siarea_sep.sel(time=slice('1900-01-01',None))
dplot['time']=tt_index

get_ipython().run_line_magic('matplotlib', 'inline')
fig=plt.figure(figsize=[15,8])
cmap=smurphs_cmap
pp=[]
for ie,exp in enumerate(exps):
    ax=plt.subplot(2,3,ie+1)
    pp.append((dplot/1e6).sel(exp=exp).plot(x='time',hue='run',color=cmap[ie],linewidth=1,add_legend=False,alpha=0.6)[0])
    (dplot/1e6).mean(dim='run').sel(exp=exp).plot(x='time',color=cmap[ie],linewidth=2,add_legend=False)
  #  plt.axhline(0,color='k')
    plt.xlim([1900,2013])
   # plt.ylim([-4,7])
    if np.mod(ie,3)==0:
        plt.ylabel('x10$^6$ km$^2$')
    else:
        plt.ylabel('')
    if ie>2:
        plt.xlabel('Y')
    else:
        plt.xlabel('')
    plt.title(exp_names[ie],fontweight='bold')
plt.suptitle('Sept Antarctic Sea Ice Extent relative to 1850-1900 mean',fontweight='bold',y=0.94)
plt.savefig('./Sept_SIE_byexp.png',bbox_inches='tight')
plt.savefig('./Sept_SIE_byexp.pdf',bbox_inches='tight')


# In[37]:


#dplot=siarea_march-siarea_march.sel(time=slice('1850-01-01','1900-01-01')).mean(dim=['run','exp'])
dplot=siarea_march.sel(time=slice('1900-01-01',None))
dplot['time']=tt_index
dplot2=siarea_sep.sel(time=slice('1900-01-01',None))
dplot2['time']=tt_index

get_ipython().run_line_magic('matplotlib', 'inline')
fig=plt.figure(figsize=[10,5])
cmap=smurphs_cmap
pp=[]
ax=fig.add_subplot(1,2,1)
for ie,exp in enumerate(exps):
    (dplot/1e6).sel(exp=exp).plot(ax=ax,x='time',hue='run',color=cmap[ie],linewidth=1,add_legend=False,alpha=0.2)
    pp.append((dplot/1e6).mean(dim='run').sel(exp=exp).plot(x='time',color=cmap[ie],linewidth=2,add_legend=False)[0])
#plt.axhline(0,color='k')
plt.xlim([1900,2013])
plt.ylabel('x10$^6$ km$^2$')
plt.xlabel('Y')
plt.title('March Antarctic Sea Ice Extent',fontweight='bold')
ax=fig.add_subplot(1,2,2)
for ie,exp in enumerate(exps):
    (dplot2/1e6).sel(exp=exp).plot(ax=ax,x='time',hue='run',color=cmap[ie],linewidth=1,add_legend=False,alpha=0.2)
    pp.append((dplot2/1e6).mean(dim='run').sel(exp=exp).plot(x='time',color=cmap[ie],linewidth=2,add_legend=False)[0])
#plt.axhline(0,color='k')
plt.xlim([1900,2013])
plt.ylabel('')
plt.xlabel('Y')
plt.legend(pp,exp_names)
plt.title('Sept Antarctic Sea Ice Extent',fontweight='bold')

plt.savefig('./SIE_all.png',bbox_inches='tight')
plt.savefig('./SIE_all.pdf',bbox_inches='tight')


# In[ ]:


get_ipython().system('jupyter nbconvert --to script SIE.ipynb')


# In[ ]:




