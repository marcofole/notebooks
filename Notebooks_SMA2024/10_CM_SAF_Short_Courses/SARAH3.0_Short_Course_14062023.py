#!/usr/bin/env python
# coding: utf-8

# <img src='https://gitlab.eumetsat.int/eumetlab/oceans/ocean-training/tools/frameworks/-/raw/main/img/CMSAF_Name_Colour.png' align='left' width='25%'/>
# <img src='https://gitlab.eumetsat.int/eumetlab/oceans/ocean-training/tools/frameworks/-/raw/main/img/eumet_logo.png' align='right' width='30%'/>

# # **SARAH - A Climate Data Record on Surface Solar Radiation**

# ## **14 June 2023 - ONLINE Short Course**

# <hr>

# ## <span style="color:blue">**Using the data - Jupyter Notebook** </span>.

# ### **Introduction**

# This notebook will introduce two products of <span style="color:red">**SARAH-3**</span>, the third release of the **SARAH** climate data record, showing how to download the data and how to extract and analyze the parameters contained in the products.

# ### **Outline**
# 
# * [**1. Loading SARAH-3 Products**](#load)
# * [**2. Reading SARAH-3 Products**](#read)
# * [**3. Data Extraction and Analysis**](#analysis)

# ### **Importing required libraries**

# In[1]:


import netCDF4 as nc
import glob
import tarfile
import xarray as xr
import tempfile
import warnings
import numpy as np
warnings.filterwarnings("ignore")
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.shapereader as shpr
from matplotlib import cm as cm
import matplotlib.pyplot as plt
import geopandas
import pandas as pd
import rioxarray as rio
from shapely.geometry import mapping
import nclcmaps
import seaborn as sns
import os


# ### <a id="load"></a>**1. Loading SARAH-3 products**

# Depending on the data that are needed by the user, there are different ways to load SARAH-3 products:

# <span style="color:blue">**Option 1** </span>: **_quick, small file size_**: load CLAAS-3 products using the corresponding **URL** (these are pre-modified files containing pre-specified region and one parameter per NetCDF file and they can accessed from a remote repository)
# The following example will load the daily Sunshine Duration (SDU) product for 1st May 2023, cut over Europe for the purpose of this short course.
# 

# In[2]:


daily_SDU_cut = xr.open_dataset('https://public.cmsaf.dwd.de/data/stkothe/SARAH_SC/SDUms202305010000004UD10001I1UD.nc'+'#mode=bytes')


# Now the product has been loaded into the variable named **daily_SDU** and it is possible to visualize its content just typing the name of the variable and running the cell:

# In[3]:


daily_SDU_cut


# The same method is applied when the products have been already uploaded directly in the remote Jupyter environment. Infact in this case the user has just to replace the **URL** with the **path** to the directory where the file is stored. For example the SDU product just loaded from an external URL is also located at the following internal path: 
# 
# `./eodata/SARA3.0_WS/SDUms202305010000004UD10001I1UD.nc`
# 
# So, we can use the same instruction shown above,but using the file path; the output is the same:

# In[8]:


daily_SDU_cut = xr.open_dataset('./eodata/SARA3.0_WS/SDUms202305010000004UD10001I1UD.nc')
daily_SDU_cut 


# <span style="color:blue">**Option 2.** </span> _**working with original products**_: if the user needs to work directly with the original products, the process to download the products directly into the local or remote jupyter environment can be described in 3 steps:
# 
# - select a product, choose the related options from the [**Web User Interface**](https://wui.cmsaf.eu/safira/action/viewHome) and submit an order
# 
# - the email address associated to the CM SAF account will receive an email with the download options  for each submitted order; each email will contain a **wget** instruction that can be used to download the product as a **.tar** file. The user can get the product running the following Python line in a cell of code:
# 
#      `os.system('<wget instruction in the email>')`
# 
#      the **.tar** file will be downloaded directly in the directory where the notebook has been created.
#      
# - unpack the **.tar** file to open the original products. 
# 
# The following example will show the whole process to get an original product ready to be used in the jupyter environment.

# #### <span style="color:orange">**Loading an original SARAH-3 product - example** </span>

# After submitting an order to get the SDU daily product (whole disk) for 1st of May 2023, an email confirming the product availability will provide the **wget** instruction for download, which will be similar to the following, valid for the example:
# 
# `wget -r -np -nH --cut-dirs=1 --reject="index.html" --user=routcm --password=4gVdHUdpq8UhHcIJIP https://cmsaf.dwd.de/data/ORD50063/`
# 
# and the Python code will be:

# In[11]:


os.system('wget -r -np -nH --cut-dirs=1 --reject="index.html" --user=routcm --password=4gVdHUdpq8UhHcIJIP https://cmsaf.dwd.de/data/ORD50063/')


# The user will find a subfolder named **ORD50063** in the same directory of the notebook. The subfolder will contain the TAR file that needs to be unpacked. 
# The code to unpack a TAR file is always the same. The user should change the path of the original product and the destination path, when necessary. In this example:

# In[15]:


my_tar = tarfile.open('ORD50063/ORD50063.tar')
my_tar.extractall() # a "work" directory is created for this product
my_tar.close()


# Now the unpacked product `SDUds2023061100000042310001I1MA.nc` is in the same folder where the notebook is located, ready to be opened as it has been shown for the NetCDF file used in the **option 1** section.

# <hr>

# ### <a id="read"></a>**2. Reading SARAH-3 products**

# #### 2.1 'Single product' case

# The previous section has described how to load and open a single SARAH-3 product that is contained in a NetCDF file. As already shown, in this case the file can be opened as follows:

# In[5]:


daily_SDU = xr.open_dataset('SDUds2023061100000042310001I1MA.nc')
daily_SDU


# The selected file is the one just downloaded in the last example, containing the daily SDU data for the whole SEVIRI disk.
# The structure of the file, besides showing the geospatial and temporal coordinates, contains five **Data Variables** which can be extracted from the file. In particular, the user should be able to extract and manipulate the <span style="color:red">**SDU** </span> variable. This is possible with the following operation:

# In[6]:


SDU_data = daily_SDU['SDU'] ##'SDU' is the name of the variable to extract
SDU_data


# The extracted variable is represented with a matrix of 2600 x 2600 pixels and has 4 attributes, instead of 39 as the original product. When the user wants to visualize a single variable for a single instantaneous/daily/monthly SARAH-3 product, the code should specify:
# 
# - the figure size
# - the colormap  (find list of available colormaps in the [**Matplotlib**](https://matplotlib.org/stable/gallery/color/colormap_reference.html) page)
# - the projection to plot the data on (find list of available projections in the [**Cartopy**](https://scitools.org.uk/cartopy/docs/v0.15/crs/projections.html) page) 
# - the *transform* parameter to reproject data
# - whether to draw **coasts** and **gridlines** on the map

# In order to get a quick and easy visualization, the user can always follow this structure:

# In[14]:


fig = plt.figure(figsize=(12,10))  #figure size
colormap = 'afmhot'                #colormap

#-----Defining the color of nodata pixels------
cmap = plt.get_cmap(colormap)
cmap.set_bad('grey')


#---------------Product visualization--------------------
ax = plt.axes(projection= ccrs.PlateCarree())   #projection 
SDU_data.plot(transform=ccrs.PlateCarree(),cmap=colormap)  #data visualization
ax.gridlines( crs=ccrs.PlateCarree(),   #gridlines options
              linewidth=0.8,
              color='black', 
              alpha=0.5, 
              linestyle='--', 
              draw_labels=True)

cmap = plt.get_cmap(colormap)
cmap.set_bad('grey')
ax.coastlines(color='black')

fig.tight_layout()

plt.savefig('SDU_DISK.png', dpi = 100, bbox_inches='tight') #saving the map as a .png file

plt.show()


# The user should have noticed that it is possible to save the resulting map as a png file: that means that each participant of the Jupyter session could produce the various maps described in the notebook and save them, so that they can be downloaded. 
# 
# It is worth mentioning that during and after the course the participants could share their maps and plots in the **CM SAF** [**PADLET**](https://padlet.com/CMSAF/the-cm-saf-padlet-azeujpu4vc9cbfcy) where results can be commented by CM SAF experts and others.

# #### 2.2 'Multiple products' case 

# Besides accessing and visualizing a single product, dealing with SARAH-3 products mainly means to manipulate a data record that can cover 40 years of data. For this course some pre-modified NetCDF files have been prepared for both SIS and SDU products, providing data only for a spatial subset of the SEVIRI full disk. This allows to skip the whole process that starts by downloading the original products. As an example, for this course the SIS full data record has been stored in the NetCDF file:
# 
# `SIS_1983-01-01-2022-12-01.nc`
# 
# which can be opened from the remote repository:

# In[37]:


SIS_month_rec = xr.open_dataset('https://public.cmsaf.dwd.de/data/stkothe/SARAH_SC/SIS_1983-01-01-2022-12-01.nc'+'#mode=bytes')
SIS_month_rec


# If compared with the structure of a single product, the relevant difference is given by the **time** coordinate that is now equal to 480,so the dataset is made of 480 stacked monthly products, covering from 1983 to 2022. The timestamps of all the products can be visible extracting the **time** coordinate:

# In[38]:


SIS_month_rec['time']

# SIS_month_rec['time'].values  to see the whole list


# Starting from the overall structure, if the user wanted to extract a single variable, it would be necessary first to select which timestamp to consider.:

# In[39]:


SIS_month = SIS_month_rec.sel(time = '2022-12-01')
SIS_month


# At this point, the SIS variable can be extracted:

# In[23]:


SIS_data = SIS_month['SIS']
SIS_data


# The user will recognize a structure that is similar to what has been obtained in the previous sections for the single SDU products, so it is possible to visualize the corresponding map as it has been done for the SDU variable:

# In[34]:


fig = plt.figure(figsize=(12,10))  #figure size
colormap = 'RdYlBu_r'                #colormap warning: the ramp of colors has been inverted 

#-----Defining the color of nodata pixels------
cmap = plt.get_cmap(colormap)
cmap.set_bad('grey')


#---------------Product visualization--------------------
ax = plt.axes(projection= ccrs.PlateCarree())   #projection 
SIS_data.plot(transform=ccrs.PlateCarree(),cmap=colormap)  #data visualization
ax.gridlines( crs=ccrs.PlateCarree(),   #gridlines options
              linewidth=0.8,
              color='black', 
              alpha=0.5, 
              linestyle='--', 
              draw_labels=True)

cmap = plt.get_cmap(colormap)
cmap.set_bad('grey')
ax.coastlines(color='black')

plt.savefig('SIS_cut_2022_12.png', dpi = 100, bbox_inches='tight') #saving the map as a .png file

plt.show()


# The choice of the colormap aims to distinguish the SDU variable from the SIS one, but in this case the selected colormap has been used after inverting the ramp of colors, in order to keep the map coherent with the values of the variable.
# From a visual point of view the size of the colorbar seems to be too large if compared with the map, so the default setting of the colorbar is not the best one for this layout. The user could modify the setting related to the colorbar by adding some instruction in the previous cell of code:

# In[35]:


fig = plt.figure(figsize=(12,10))  #figure size
colormap = 'RdYlBu_r'                #colormap warning: the ramp of colors has been inverted 

#-----Defining the color of nodata pixels------
cmap = plt.get_cmap(colormap)
cmap.set_bad('grey')


#---------------Product visualization--------------------
ax = plt.axes(projection= ccrs.PlateCarree())   #projection 


#--------modifying the colorbar settings--------------------
map_eur = SIS_data.plot(ax = ax, add_colorbar=False,cmap = cmap) #data visualization
cbar = plt.colorbar(map_eur, shrink = 0.8, pad = 0.1, orientation='horizontal',label='solar surface irradiance [W/m2]')



ax.gridlines( crs=ccrs.PlateCarree(),   #gridlines options
              linewidth=0.8,
              color='black', 
              alpha=0.5, 
              linestyle='--', 
              draw_labels=True)

cmap = plt.get_cmap(colormap)
cmap.set_bad('grey')
ax.coastlines(color='black')
ax.set_title('Monthly solar surface irradiance (SIS) 2022/12')

#fig.tight_layout()

plt.savefig('SIS_cut_2022_12.png', dpi = 100, bbox_inches='tight') #saving the map as a .png file

plt.show()


# #### 2.3 'Multiple products' case (**OPTIONAL**)

# When the submitted order returns a .TAR file that contains multiple NetCDF elements, in case of long term analysis that should produce anomaly maps or extract information from a time series, reading product by product is not the proper way to go. 
# The **xarray** library of Python provides the **open_mfdataset** function that can read multiple NetCDF files at the same time, if they are all in the same directory. The result is a new **xarray.Dataset** which includes the data from all the files stored in the source directory. 

# #### <span style="color:orange">**Loading multple original SARAH-3 products - the CDR and ICDR example** </span>

# If the user is working on the online Jupyter platform, inside the **eodata** directory, the path **SARA3.0_WS/SIS** folder contains the subfolders **ORD49903** and **ORD49903_2** providing the CDR and ICDR, respectively, of the **Solar Surface Irradiance** (SIS) product.
# 
# If the user is working on a local machine, he could download the two .TAR files related to the two subfolders from the repository mentioned at the beginning of the notebook. Then it will be necessary to unpack the two .TAR files as shown in **section 1** , for example storing all the NetCDF files in a new folder named **SIS**.
# 
# In both of the cases it will be possible to read all the files inside the subfolders with the following instruction:

# In[16]:


SIS_month_rec = xr.open_mfdataset('./eodata/SARA3.0_WS/SIS/**/*.nc', concat_dim='time', combine='nested')


# The string `./eodata/SARA3.0_WS/SIS/**/*.nc` includes all the files that are contained in subfolders of the SIS directory and with the NetCDF extension. The contents of the files will be then temporally sorted and with the coordinates properly aligned thanks to the **open_mfdataset**. And the new combined xarray.Dataset is structured as follows:

# In[17]:


SIS_month_rec


# The code to extract a single product and properly visualize it has been described in the previous section, so here it will be merely replicated just to confirm the same map in output:

# In[40]:


SIS_month = SIS_month_rec.sel(time = '2022-12-01')

SIS_data = SIS_month['SIS']

fig = plt.figure(figsize=(12,10))  #figure size
colormap = 'RdYlBu_r'                #colormap warning: the ramp of colors has been inverted 

#-----Defining the color of nodata pixels------
cmap = plt.get_cmap(colormap)
cmap.set_bad('grey')


#---------------Product visualization--------------------
ax = plt.axes(projection= ccrs.PlateCarree())   #projection 


#--------modifying the colorbar settings--------------------
map_eur = SIS_data.plot(ax = ax, add_colorbar=False,cmap = cmap) #data visualization
cbar = plt.colorbar(map_eur, shrink = 0.8, pad = 0.1, orientation='horizontal',label='solar surface irradiance [W/m2]')



ax.gridlines( crs=ccrs.PlateCarree(),   #gridlines options
              linewidth=0.8,
              color='black', 
              alpha=0.5, 
              linestyle='--', 
              draw_labels=True)

cmap = plt.get_cmap(colormap)
cmap.set_bad('grey')
ax.coastlines(color='black')
ax.set_title('Monthly solar surface irradiance (SIS) 2022/12')

#fig.tight_layout()

#plt.savefig('SIS_cut_2022_12.png', dpi = 100, bbox_inches='tight') #saving the map as a .png file

plt.show()


# <hr>

# ### <a id="analysis"></a>**3. Data Extraction and Analysis**

# Three different tasks are proposed to the users, to highlight which could be the processings that can be applied to the SARAH-3 products and which results can be obtained.

# #### **3.1 - Task 1: map of the long-term monthly mean of SIS for June, with reference period 1983-2020, on Europe cutout**
# 
# - extract the variable of interest
# - select the time range of interest
# - group data by month and compute the long term mean for each month 
# - extract the data for the required month
# - plot the map for the required month

# From the previous sections the whole SIS data record is already available, as an xarray.DataArray object named **SIS_month_rec** :

# In[44]:


SIS_month_rec


# **Step1**. The temporal coverage spans from 1983 to 2022, so it is necessary to extract the required temporal range **1983-2020**

# In[47]:


start_time = '1983-01-01'
end_time = '2020-12-31'

SIS_data = SIS_month_rec['SIS']
SIS_data_subset = SIS_data.sel(time = slice(start_time,end_time))  #extraction of the required temporal range
SIS_data_subset


# The temporal verage has ben set as required, now the long term monthly mean can be computed, in one line , for each of the 12 months.

# **Step 2.** Long term monthly mean calculation

# In[48]:


SIS_monthly_mean = SIS_data_subset.groupby('time.month').mean(dim='time')
SIS_monthly_mean


# Now the SIS data have been grouped by month and then the monthly mean has been calculated for each month across the temporal coverage, so the result provides the mean values for each month on the spatial coverage of the product. It is now required to select the month of June.

# **Step 3.** Single month extraction (June)

# In[49]:


SIS_mean_June = SIS_monthly_mean.sel(month = 6)
SIS_mean_June


# A single month has been extracted and now the corresponding map can be visualized.

# **Step 4.** Data visualization

# In[50]:


fig = plt.figure(figsize=(12,10))  #figure size
colormap = 'RdYlBu_r'                #colormap warning: the ramp of colors has been inverted 

#-----Defining the color of nodata pixels------
cmap = plt.get_cmap(colormap)
cmap.set_bad('grey')


#---------------Product visualization--------------------
ax = plt.axes(projection= ccrs.PlateCarree())   #projection 


#--------modifying the colorbar settings--------------------
map_eur = SIS_mean_June.plot(ax = ax, add_colorbar=False,cmap = cmap) #data visualization
cbar = plt.colorbar(map_eur, shrink = 0.8, pad = 0.1, orientation='horizontal',label='solar surface irradiance [W/m2]')



ax.gridlines( crs=ccrs.PlateCarree(),   #gridlines options
              linewidth=0.8,
              color='black', 
              alpha=0.5, 
              linestyle='--', 
              draw_labels=True)

cmap = plt.get_cmap(colormap)
cmap.set_bad('grey')
ax.coastlines(color='black')
ax.set_title('Long term SIS monthly mean for June  Ref.Period (1983-2020)' )

plt.savefig('SIS_June_mean.png', dpi = 100, bbox_inches='tight') #saving the map as a .png file

plt.show()


# #### **3.2 - Task 2: Plot time series of monthly SIS for a selected location in the period 1983-2020**

# - define geografical coordinates of the selected point
# - find a SARAH-3 grid box with defined coordinates 
# - plot time series of monthly SIS over the requested period

# **Step 1.** Define the geographical coordinates of the selected point

# In[56]:


Location = 'Offenbach'
Latitude = 50.1
Longitude = 8.75 


# **Step 2.** Extract the time series for the point location using a SARAH-3 grid box

# It is possible to use directly the SIS data record that has been already extracted for 1983-2020 in the previous task

# In[57]:


SIS_Offenbach = SIS_data_subset.sel(lat = Latitude, lon = Longitude, method = 'nearest')
SIS_Offenbach


# As the time series is referred to a single point, now the only coordinate that defines the variabilty of the parameter is the time.

# **Step 3.** Plot the time series for the requested period.

# In[61]:


fig = plt.figure(figsize=(14,8))
ax = plt.axes()
SIS_Offenbach.plot(ax=ax, label='monthly SIS')
ax.set_xlabel('Year')
ax.grid(alpha = 0.3)
ax.legend()
ax.set_title('Solar Surface Irradiance over '+Location+' (1983 - 2020).')
plt.savefig('SIS_Offenbach_TS.png', dpi = 100, bbox_inches='tight')
plt.show()


# #### **3.3 - Task 3: Monthly anomaly of SDU for May 2023 Ref.Period (1991-2020)**

# - Load and open the SDU data record
# - Extract the SDU data for the required reference period
# - Calculate the monthly long term mean for May
# - Extract the SDU data for May 2023
# - Calculate the monthly anomaly for May as **data of May 2023** - **May long term monthly mean**
# - Visualize the May monthly anomaly map of SDU

# **Step 1.** The SDU data are contained  in the predefined NetCDF file `SDU_1991-01-01-2020-12-01.nc` that can be directly loaded from the remote repository provided by CM SAF

# In[63]:


SDU_month_rec = xr.open_dataset('https://public.cmsaf.dwd.de/data/stkothe/SARAH_SC/SDU_1991-01-01-2020-12-01.nc'+'#mode=bytes')
SDU_month_rec


# **Step 2.** Extract the data for the reference period of interest (1991-2020)

# In[64]:


start_date = '1991-01-01'
end_date = '2020-12-31'

SDU_data = SDU_month_rec['SDU']   #variable extraction
SDU_data_subset = SDU_data.sel(time = slice(start_time,end_time))  #extraction of the required temporal range
SDU_data_subset


# **Step 3.** Calculate the long term monthly mean for May. It is the same operation already implemented for the SIS data in Task-1

# In[89]:


SDU_monthly_mean = SDU_data_subset.groupby('time.month').mean(skipna = False)
SDU_monthly_mean


# From the stack of the twelve monthly means, it is possible to extract May:

# In[93]:


SDU_mean_May = SDU_monthly_mean.sel(month = 5)
SDU_mean_May


# **Step 4.** Extract the SDU data for May 2023.

# It is necessary to read the data directly from the remote repository.

# In[109]:


SDU_May_2023 = xr.open_dataset('https://public.cmsaf.dwd.de/data/stkothe/SARAH_SC/SDUms202305010000004UD10001I1UD.nc'+'#mode=bytes')
SDU_May_2023
SDU_May_2023_data = SDU_May_2023['SDU']


# **Step 5.** Calculate the monthly anomaly 

# In[110]:


SDU_May23_anom = SDU_May_2023_data - SDU_mean_May
SDU_May23_anom


# **Step 6.** Visualize the monthly anomaly map

# In[112]:


fig = plt.figure(figsize=(12,10))  #figure size
colormap = 'coolwarm'                #colormap warning: the ramp of colors has been inverted 

#-----Defining the color of nodata pixels------
cmap = plt.get_cmap(colormap)
cmap.set_bad('grey')


#---------------Product visualization--------------------
ax = plt.axes(projection= ccrs.PlateCarree())   #projection 


#--------modifying the colorbar settings--------------------
map_eur = SDU_May23_anom.plot(ax = ax, add_colorbar=False,cmap = cmap) #data visualization
cbar = plt.colorbar(map_eur, shrink = 0.8, pad = 0.1, orientation='horizontal',label='SDU [h]')



ax.gridlines( crs=ccrs.PlateCarree(),   #gridlines options
              linewidth=0.8,
              color='black', 
              alpha=0.5, 
              linestyle='--', 
              draw_labels=True)

cmap = plt.get_cmap(colormap)
cmap.set_bad('grey')
ax.coastlines(color='black')
ax.set_title('SDU Monthly anomaly for May 2023' )

plt.savefig('SDU_May_anom.png', dpi = 100, bbox_inches='tight') #saving the map as a .png file

plt.show()


# <img src='https://gitlab.eumetsat.int/eumetlab/oceans/ocean-training/tools/frameworks/-/raw/main/img/CMSAF_Name_Colour.png' align='left' width='20%'/>
# <img src='https://gitlab.eumetsat.int/eumetlab/oceans/ocean-training/tools/frameworks/-/raw/main/img/eumet_logo.png' align='right' width='20%'/>
