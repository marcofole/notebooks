#!/usr/bin/env python
# coding: utf-8

# <img src='https://gitlab.eumetsat.int/eumetlab/oceans/ocean-training/tools/frameworks/-/raw/main/img/CMSAF_Name_Colour.png' align='left' width='25%'/>
# <img src='https://gitlab.eumetsat.int/eumetlab/oceans/ocean-training/tools/frameworks/-/raw/main/img/eumet_logo.png' align='right' width='30%'/>

# # **CLAAS - A Climate Data Record of Cloud Properties**
# # _19 April 2023 - ONLINE Short Course_
# ## <span style="color:blue">**Using the data - Jupyter Notebook** </span>

# <hr>

# ## **Introduction**

# This notebook
# * demonstrates how to load CLAAS-3 data into a jupyter notebook;
# * shows how to handle CLAAS-3 products and extract parameters from the data.
# * provides examples for plotting and analysing CLAAS-3 products.

# ### **Products in use**

# The following **CLAAS-3** products will be used in this notebook:
# 
# | **Product Name**                               | Product Family | Area | Temp. Resolution | Spat. Resolution (degrees) | Stat. |
# |--------------------------------------------| ----| ---  |:---: | :---: | --- |
# | [<span style="color:red">**CFC - Fractional Cloud Cover**</span>](https://wui.cmsaf.eu/safira/action/viewProduktDetails?eid=22211_22232&fid=38) | <span style="color:blue">**CLAAS ed. 3.0**</span> | <span style="color:blue">**METEOSAT full disk**</span> | <span style="color:blue">**Daily**</span> | <span style="color:blue">**0.05 x 0.05**</span> | <span style="color:blue">**Mean**</span> |
# | [<span style="color:red">**CTO - Cloud Top Pressure**</span>](https://wui.cmsaf.eu/safira/action/viewProduktDetails?eid=22220_22241&fid=38) | <span style="color:blue">**CLAAS ed. 3.0**</span> | <span style="color:blue">**METEOSAT full disk**</span> | <span style="color:blue">**Monthly**</span> | <span style="color:blue">**0.05 x 0.05**</span> | <span style="color:blue">**Mean**</span> |
# | [<span style="color:red">**CPH - Liquid Cloud Fraction**</span>](https://wui.cmsaf.eu/safira/action/viewProduktDetails?eid=22217_22238&fid=38) | <span style="color:blue">**CLAAS ed. 3.0**</span> | <span style="color:blue">**METEOSAT full disk**</span> | <span style="color:blue">**Monthly**</span> | <span style="color:blue">**0.05 x 0.05**</span> | <span style="color:blue">**Mean diurnal-cycle**</span> |

# ### Link to this notebook, used files and list of python pakages:
# https://public.cmsaf.dwd.de/data/isolodov/CLAAS3_workshop

# ## **Importing required libraries**

# In[1]:


import os #launching Terminal instruction in the notebook
import tarfile #unpacking .tar files
import xarray as xr # reading NetCDF files
import numpy as np # computations
import pandas as pd # computations
from scipy import stats # statistics
import matplotlib.pyplot as plt # plotting
import cartopy.crs as ccrs # georeferencing
import cartopy.feature as cfeature
import calendar
import warnings
warnings.filterwarnings("ignore")


# <hr>

# ## **Load CLAAS-3 products**

# **Option 1** *quick, small file size*: Load CLAAS-3 products using URL (these are pre-modified files containing pre-specified region and one parameter per NetCDF file)
# 
# **Option 2** *working with original products, files can be large*: Load CLAAS-3 products which were directly ordered through CM SAF Web User Interface
# 
# We will load three different products:
# 1) level-3 daily mean of <u>cloud fraction</u> (CFC) for the whole SEVIRI disk
# 2) level-3 monthly means of <u>cloud top pressure</u> (CTO)
# 3) level-3 monthly mean <u>diurnal cycle of liquid cloud fraction</u> (CPH)

# #### **Option 1**
# 
# The selected products are made available via permanent link to a CM SAF server where the data are stored.

# In[2]:


# 1)
ds_l3_dm_cfc = xr.open_dataset('https://public.cmsaf.dwd.de/data/perm/training/CLAAS_SC/cfc_2023-04-06-2023-04-06.nc'+'#mode=bytes')


# In[10]:


#2)
ds_l3_mm_ctp=  xr.open_dataset('https://public.cmsaf.dwd.de/data/perm/training/CLAAS_SC/ctp_2004-01-01-2023-03-01.nc'+'#mode=bytes')


# In[18]:


#3)
ds_l3_cph_md = xr.open_dataset('https://public.cmsaf.dwd.de/data/perm/training/CLAAS_SC/cph_mmdc_2004-01-01-2023-03-01.nc'+'#mode=bytes')


# look into the loaded data set

# In[19]:


ds_l3_cph_md


# #### **Option 2**

# The user has the possibility to download directly the products ordered in the [**CM SAF Web User Interface**](https://wui.cmsaf.eu/safira/action/viewProduktSearch) into the Jupyter environment, using one of the options that provided in the email reiceved by the user, after ordering the product. For example using the **https** option the download instruction is the following, for the ICDR CTO product:
# 
# <span style="color:red">**Important:**</span>
# **Do not run the instruction here below during the hands-on session of the short course, data are already available for the users. This part is just to show how a user can get the data once he/she has ordered them in the official CM SAF web portal, for further applications/analysis.**

# In[ ]:


os.system('wget -r -np -nH --cut-dirs=1 --reject="index.html" --user=routcm --password=4gVdHUdpq8UhHcIJIP https://cmsaf.dwd.de/data/ORD49357/')


# Once the data have been downloaded, in the **notebooks** section the user will find a **ORD49357** folder: inside the folder a **.tar** file is present. To actually use the ordered data it is necessary to unpack this file:
# 
# `/home/jovyan/notebooks/ORD49357/ORD49357.tar`
# 
# and with the following instruction the data will be unpacked and stored in the **work** directory

# In[24]:


my_tar = tarfile.open('./ORD49357/ORD49357.tar')
my_tar.extractall('../work') # specify which folder to extract to
my_tar.close()


# It is possible to read all the **netcdf** files contained in a directory and to concatenate them in temporal order with the following instruction:

# In[27]:


data_test = xr.open_mfdataset('../work/*.nc', concat_dim='time', combine='nested')
data_test


# <span style="color:red">**Important:**</span> from this point onwards, the operations will use only the data retrieved with **option 1** which is easier and faster for the purpose of the short course.

# <hr>

# ## **Data extraction and analysis**

# ### **Demo tasks**
# 1) plot a map of cloud fraction (CFC) using a level-3 daily product;
# 2) plot a map of temporal averaged monthly mean cloud top pressure (CTP) over a specific region (Iberian Peninsula); make a time series of CTP over one geographical point (Madrid);
# 3) plot diurnal cycle of fraction of liquid clouds (CPH) for a specific geographical point (Madrid); plot CPH at midday over Europe.

# ### 1) Plot cloud fraction using daily mean data on 2020/07/15

# Exctract cfc (cloud fraction) variable from ds_l3_dm_cfc

# In[4]:


cfc_daily_mean = ds_l3_dm_cfc['cfc']


# Show content of the variable

# In[5]:


cfc_daily_mean


# Plot daily cloud fraction in geostationary projection:
# - specify figure size
# - specify colormap
# - choose projection to plot the data on
# - use *transform* parameter to reproject data
# - draw coasts and gridlines

# In[6]:


fig = plt.figure(figsize=(8,6))
colormap = 'Blues'

ax = plt.axes(projection= ccrs.Geostationary())

cfc_daily_mean.plot(transform=ccrs.PlateCarree(),cmap=colormap)

ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.8,
              color='black', alpha=0.5, linestyle='--', draw_labels=True)
ax.coastlines(color='black')
fig.tight_layout()
plt.show()


# ### 2) draw temporal-averaged map of monthly/seasonal CTP over the Iberian Peninsula and make time series of monthly CTP over a specified location

# #### 2-1) Plot temporal average of monthly mean CTP

# Extract 'ctp' variable

# In[12]:


ctp_mm = ds_l3_mm_ctp['ctp']    
ctp_mm


# Average cloud top pressure over the Iberian Peninsula and plot it in a PlateCarree projection (equidistant cylindrical projection)

# In[13]:


ctp_avg = ctp_mm.mean(dim='time')

fig=plt.figure(figsize=(10,6))
ax = plt.axes(projection=ccrs.PlateCarree())
map_eur = ctp_avg.plot(ax = ax, add_colorbar=False)
cbar = plt.colorbar(map_eur,shrink = 0.8, pad = 0.1, orientation='horizontal',label='cloud top pressure [hPa]')
ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.8,
              color='black', alpha=0.5, linestyle='--', draw_labels=True)
ax.coastlines(color='black')
ax.set_title('Cloud Top Pressure 2004 - 2023/03')
plt.show()


# #### 2-2) Compute seasonal averages, plot and compare

# Group data set by seasons and compute an average of each 3 months within a season

# In[14]:


ctp_avg_seasons = ctp_mm.groupby('time.season').mean(dim="time")
ctp_avg_seasons


# Plot seasonal averages as subplots

# In[15]:


fig, axes = plt.subplots(ncols=2,nrows=2, figsize=(12,8), subplot_kw={'projection': ccrs.PlateCarree()}, constrained_layout=True)
for i, season in enumerate(("DJF", "MAM", "JJA", "SON")):
    ax_i = axes.flatten()[i]
    ctp_avg_seasons.sel(season=season).plot(ax=ax_i, transform=ccrs.PlateCarree(), add_colorbar=True,
                                           vmin=300,vmax=700)
    ax_i.coastlines(color='black')
plt.show()


# #### 2-3) Plot time series of cloud top pressure over Madrid

# In[16]:


# Choose grid box with Madrid
location = 'Madrid'
latitude = 40.4189
longitude = -3.6919

ctp_mm_1gridbox = ctp_mm.sel(lat= latitude, lon= longitude, method='nearest')
ctp_mm_1gridbox


# In[17]:


fig = plt.figure(figsize=(10,6))
ax = plt.axes()
ctp_mm_1gridbox.plot(ax=ax, label='monthly CTP')
ax.grid()
ax.legend()
ax.set_title('Cloud top pressure over '+location+' Monthly means 2004 - 2023/03.')


# ### 3) plot diurnal cycle of fraction of liquid clouds (CPH) for a specifid location and plot CPH at midday over Europe

# #### 3-1) Plot long-term mean diurnal Cycle over Madrid

# Exctract 'cph_mmdc' variable from ds_l3_cph_md

# In[20]:


cph_md = ds_l3_cph_md['cph_mmdc']
cph_md


# Find a gridbox with Madrid

# In[21]:


location = 'Madrid'
latitude = 40.4189
longitude = -3.6919

cph_md_1gridbox = cph_md.sel(lat= latitude, lon= longitude, method='nearest')


# Compute a long-term average for each hour of the day

# In[22]:


cph_md_1gridbox_avg = cph_md_1gridbox.groupby('time.hour').mean()
cph_md_1gridbox_avg


# The prepared product has monthly averaged CPH values for each hour of the day, just plot it as a line.

# In[23]:


fig = plt.figure(figsize=(10,6))
ax = plt.axes()
cph_md_1gridbox_avg.plot(ax=ax, label='monthly CTP')
ax.grid()
ax.legend()
plt.title('Fraction of liquid cloud cover over '+location+'\nMonthly Mean Diurnal Cycle 2004-2023')
plt.show()


# #### 3-2) Plot a map of liquid cloud fraction at midday

# Extract monthly CPH at midday (12 UTC) 

# In[24]:


cph_md_12 = cph_md.sel(time=(cph_md.time.dt.hour == 12))


# Average all monthly mean values at midday

# In[25]:


cph_md_12_avg = cph_md_12.mean(dim='time')


# Plot a map of CPH

# In[26]:


fig=plt.figure(figsize=(10,6))
ax = plt.axes(projection=ccrs.PlateCarree())
colormap = 'Spectral'
map_eur = cph_md_12_avg.plot(ax = ax, transform=ccrs.PlateCarree(), cmap=colormap, add_colorbar=False)
cbar = plt.colorbar(map_eur,shrink = 0.8, pad = 0.1, orientation='horizontal',label='fraction of liquid clouds [%]')
ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.8,
              color='black', alpha=0.5, linestyle='--', draw_labels=True)
ax.coastlines(color='black')
plt.title('Fraction of liquid cloud cover\nMonthly Mean Diurnal Cycle at 12UTC 2004-2023')
plt.show()


# <hr>

# ### **Further tasks**

# #### This section contains further tasks for analysing CLAAS-3 data. The same products and the same variables can be used.

# ### 1.1) Plot cloud fraction over Europe or another region of your choice using daily mean data on 2020/07/15
# - extract variable to plot
# - cut region you interested in
# - plot on a map

# In[27]:


cfc_daily_mean = ds_l3_dm_cfc['cfc']


# In[28]:


lats_europe = slice(34, 72)
lons_europe = slice(-25, 45)
cfc_daily_mean_eur = cfc_daily_mean.sel(lat=lats_europe).sel(lon=lons_europe)


# In[29]:


fig = plt.figure(figsize=(8,6))
colormap = 'Blues'

ax = plt.axes(projection= ccrs.Geostationary())

cfc_daily_mean_eur.plot(transform=ccrs.PlateCarree(),cmap=colormap)

ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.8,
              color='black', alpha=0.5, linestyle='--', draw_labels=True)
ax.coastlines(color='black')
fig.tight_layout()
plt.show()


# ### 1.2) Plot cloud fraction over Europe or another region of your choice using different projection using daily mean data on 2020/07/15
# - extract variable to plot
# - cut region you interested in
# - define a projection for plot, e.g. Stereographic (list of projections: https://scitools.org.uk/cartopy/docs/v0.15/crs/projections.html)
# - plot on a map

# In[30]:


cfc_daily_mean = ds_l3_dm_cfc['cfc']


# In[31]:


lats_europe = slice(34, 72)
lons_europe = slice(-25, 45)
cfc_daily_mean_eur = cfc_daily_mean.sel(lat=lats_europe).sel(lon=lons_europe)


# In[32]:


fig = plt.figure(figsize=(8,6))
colormap = 'Blues'

ax = plt.axes(projection= ccrs.Stereographic())

cfc_daily_mean_eur.plot(transform=ccrs.PlateCarree(),cmap=colormap)

ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.8,
              color='black', alpha=0.5, linestyle='--', draw_labels=True)
ax.coastlines(color='black')
fig.tight_layout()
plt.show()


# ### 2.1) Draw temporal-averaged map of CTP over the Iberain Peninsula for each month
# - extract variable to plot
# - group data by month and compute an average of all monthly means on a specific month
# - make a plot of 12 subplots for each month: use the same limits for each subplot, use months' names to title subplots

# In[33]:


ctp_mm = ds_l3_mm_ctp['ctp']


# In[35]:


ctp_avg_months = ctp_mm.groupby('time.month').mean(dim="time")
ctp_avg_months


# In[36]:


fig, axes = plt.subplots(ncols=3,nrows=4, figsize=(12,12), subplot_kw={'projection': ccrs.PlateCarree()}, constrained_layout=True)
for i, mon in enumerate(range(1,13)):
    ax_i = axes.flatten()[i]
    ctp_avg_months.sel(month=mon).plot(ax=ax_i, transform=ccrs.PlateCarree(), add_colorbar=True,
                                           vmin=300,vmax=700)
    ax_i.gridlines(crs=ccrs.PlateCarree(), linewidth=0.8,
                  color='black', alpha=0.5, linestyle='--', draw_labels=False)
    ax_i.coastlines(color='black')
    ax_i.set_title(calendar.month_name[mon])


# ### 2.2) Plot time series of cloud top pressure over a selected point
# - define geografical coordinates of point of your choise on the Iberian Peninsula, here Valencia
# - find a CLAAS-3 grid box with defined coordinates 
# - plot time series of CTP

# In[37]:


location = 'Valencia'
latitude = 39.466667
longitude = -0.375

ctp_mm_1gridbox = ctp_mm.sel(lat= latitude, lon= longitude, method='nearest')


# In[38]:


fig = plt.figure(figsize=(10,6))
ax = plt.axes()
ctp_mm_1gridbox.plot(ax=ax, label='monthly CTP')
ax.set_xlabel('Year')
ax.grid()
ax.legend()
ax.set_title('Cloud top pressure over '+location+' 2004 - 2020.')


# ### 2.3) Compute annual mean and plot it together with the monthly-based time series of cloud top pressure over a selected point
# - choose only complete years
# - compute annual average
# - plot monthly and annualy time series of CTP in one plot

# In[39]:


ctp_mm_1gridbox = ctp_mm_1gridbox.sel(time=slice('2004','2022'))
ctp_mm_1gridbox


# In[40]:


ctp_mm_1gridbox_annual = ctp_mm_1gridbox.resample(time="Y").mean() # resample == groupby 
ctp_mm_1gridbox_annual


# In[41]:


fig = plt.figure(figsize=(10,6))
ax = plt.axes()
ctp_mm_1gridbox.plot(ax=ax, label='monthly CTP')
ctp_mm_1gridbox_annual[:-1].plot(ax=ax, label='annual mean CTP', color='tab:orange', marker='o')
ax.set_xlabel('Year')
ax.grid()
ax.legend()
ax.set_title('Cloud top pressure over '+location+' 2004 - 2022.')


# ### 2.4) Weighted annual average and linear trend over a selected location
# - define annual mean by a proper weighting operation of monthly values
#     - find length of each month in the time range 
#     - compute the weight to associate to each month
#     - compute annual mean applying weights
# - compute linear trend
#     - retrieve *year* variable
#     - find the slope and intercept using fuction *linregress* from Scipy library
#     - plot annual means and linear trend in one plot

# In[42]:


# ---- getting length of each month ----
month_length = ctp_mm_1gridbox.time.dt.days_in_month
month_length


# In[43]:


# ---- computing weights ----
wgts = month_length.groupby("time.year") / month_length.groupby("time.year").mean()
wgts


# In[44]:


# ---- computing weighted average ----
ctp_mean = (ctp_mm_1gridbox * wgts).resample(time="AS").mean(dim="time")
ctp_mean


# In[45]:


# ---- variable Year extraction----
dates = pd.DatetimeIndex(ctp_mean.time)
dates.year


# In[46]:


# ---- retrieving the linear trend ----
y = ctp_mean.values
x = np.array(dates.year)
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
ann_lin_trend = intercept + slope*x


# In[47]:


#----- data visualization ------
plt.figure(figsize=(10,6))
plt.plot(dates.year,ctp_mean.values, label = 'annual mean',color = 'tab:orange', marker='o')
plt.plot(dates.year,ann_lin_trend, label = 'linear trend', color = 'tab:red')
plt.title('Cloud top pressure annual means over '+location, fontsize = 12)
plt.xlabel('Year')
plt.ylabel('CTP [hPa]')
plt.grid()
plt.legend()
plt.show()


# ### 2.5) Monthly anomalies and linear trend over a selected location
# - group the available data by month
# - calculate anomalies as the subtraction : **monthly values** - **monthly climatology** 
# - remove not-valid values
# - calculate linear trend
# - plot anomalies and linear trend in one plot

# In[48]:


monthly_values = ctp_mm_1gridbox.groupby('time.month')
monthly_climatology = ctp_mm_1gridbox.groupby('time.month').mean()


# In[49]:


anomalies = monthly_values - monthly_climatology
anomalies


# In[50]:


# ---- retrieving the linear trend ----
y = anomalies.values[1:]
x_time = anomalies.time[1:]
x = np.arange(0,np.shape(y)[0],1,dtype = int)
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
lin_trend = intercept + slope*x


# In[51]:


#----- data visualization ------
plt.figure(figsize= (10,6))
plt.xlabel('year')
plt.ylabel('CTP [hPa]')
plt.plot(x_time,y,color='tab:blue', label = 'anomalies')
plt.plot(x_time,lin_trend,color='tab:red', label = 'linear trend')
plt.title('Cloud Top Pressure - Monthly anomalies over '+location, fontsize=12)
plt.grid()
plt.legend()
plt.show()


# ### 3.1) Compare long-term monthly mean diurnal cycle of CPH over specifies location in January and July 
# - extract variable to plot
# - specify location
# - choose values for July and for January
# - average values
# - plot both diurnal cycles in one plot

# In[52]:


cph_md = ds_l3_cph_md['cph_mmdc']


# In[53]:


location = 'Valencia'
latitude = 39.466667
longitude = -0.375
cph_md_1gridbox = cph_md.sel(lat= latitude, lon= longitude, method='nearest')


# In[54]:


cph_md_1gridbox_july = cph_md_1gridbox.sel(time=(cph_md.time.dt.month == 7))
cph_md_1gridbox_jan = cph_md_1gridbox.sel(time=(cph_md.time.dt.month == 1))


# In[55]:


cph_md_1gridbox_july_avg = cph_md_1gridbox_july.groupby('time.hour').mean()
cph_md_1gridbox_jan_avg = cph_md_1gridbox_jan.groupby('time.hour').mean()


# In[56]:


fig = plt.figure(figsize=(10,6))
ax = plt.axes()
cph_md_1gridbox_jan_avg.plot(ax=ax, label='January')
cph_md_1gridbox_july_avg.plot(ax=ax, label='July')
ax.grid()
ax.legend()
plt.title('Fraction of liquid cloud cover over '+location+'\nMonthly Mean Diurnal Cycle')
plt.show()


# ### 3.2) Plot a map of long-term mean liquid cloufd fraction at midnight (0 UTC)
# - extract monthly CPH at midnight (0 UTC) 
# - average all monthly mean values at midnight
# - plot a map of CPH

# In[57]:


cph_md_00 = cph_md.sel(time=(cph_md.time.dt.hour == 0))


# In[58]:


cph_md_00_avg = cph_md_00.mean(dim='time')


# In[59]:


fig=plt.figure(figsize=(10,6))
ax = plt.axes(projection=ccrs.PlateCarree())
colormap = 'Spectral'
map_eur = cph_md_00_avg.plot(ax = ax, transform=ccrs.PlateCarree(), cmap=colormap, add_colorbar=False)
cbar = plt.colorbar(map_eur,shrink = 0.8, pad = 0.1, orientation='horizontal',label='fraction of liquid clouds [%]')
ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.8,
              color='black', alpha=0.5, linestyle='--', draw_labels=True)
ax.coastlines(color='black')
plt.title('Fraction of liquid cloud cover\nMonthly Mean Diurnal Cycle at 00UTC 2004-2023')
plt.show()


# <hr>

# <img src='https://gitlab.eumetsat.int/eumetlab/oceans/ocean-training/tools/frameworks/-/raw/main/img/CMSAF_Name_Colour.png' align='left' width='25%'/>
# <img src='https://gitlab.eumetsat.int/eumetlab/oceans/ocean-training/tools/frameworks/-/raw/main/img/eumet_logo.png' align='right' width='30%'/>
