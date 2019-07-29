# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 09:39:16 2019

@author: merks004
"""

from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import ticker
from os import listdir
from os.path import isfile, join
from matplotlib import rcParams
rcParams['font.family'] = 'Arial'

## make colormap with midpoint 0 white
class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


file_path = r'D:/UserData/merks004/Nexus/qtot/LPJmL/reldiff/ensmean'
file_list = listdir(file_path)
save_path = '//WURNET.NL/Homes/merks004/My Documents/WUR/Nexus/images'
#save_path2 = '//WURNET.NL/Homes/merks004/My Documents/WUR/C3S/Glorius/paper/figs'

nc = NetCDFFile(file_path +'//'+ file_list[0])

tas  = nc.variables['qtot']
lons = nc.variables['lon'][:]
lats = nc.variables['lat'][:]

tas = np.average(tas,axis=0)

fig, axes = plt.subplots(3,2, figsize=(9,7.5))
plt.subplots_adjust(top=0.93, bottom=0.01, left=0.02, right=0.92, hspace=0.10,
                    wspace=0.02)
# %%
ax1 = plt.subplot(3,2,1)

## load basemap
m = Basemap(projection='cyl',llcrnrlat=-60,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')


lon, lat = np.meshgrid(lons, lats)
x, y = m(lon, lat)

## draw contours
clevs = np.arange(-30.,80., 5)
cs = m.contourf(x,y,np.squeeze(tas),clevs,vmin=-31,norm=MidpointNormalize(midpoint=0),extend='both',cmap='RdBu')

## draw countries & coastlines
##m.drawcountries()
m.drawcoastlines()


plt.title('RCP2.6 2011-2040')

cbaxes = fig.add_axes([0.93, 0.68, 0.01, 0.22]) 
cb = plt.colorbar(cs, cax = cbaxes)
cb.set_label('[%]')

#%%
nc = NetCDFFile(file_path +'//'+ file_list[1])
tas  = nc.variables['qtot']
tas = np.average(tas,axis=0)

ax2 = plt.subplot(3,2,3)

## load basemap
m = Basemap(projection='cyl',llcrnrlat=-60,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')

lon, lat = np.meshgrid(lons, lats)
x, y = m(lon, lat)

## draw contours
clevs = np.arange(-30.,80., 5)
cs = m.contourf(x,y,np.squeeze(tas),clevs,vmin=-31,norm=MidpointNormalize(midpoint=0),extend='both',cmap='RdBu')

## draw countries & coastlines
##m.drawcountries()
m.drawcoastlines()


plt.title('RCP2.6 2041-2070')

## add colorbar future
#cbaxes2 = fig.add_axes([0.93, 0.15, 0.01, 0.35]) 
#cb2 = plt.colorbar(cs, cax = cbaxes2)
#cb2.set_label('[%] change')

#%%
nc = NetCDFFile(file_path +'//'+ file_list[2])
tas  = nc.variables['qtot']
tas = np.average(tas,axis=0)

ax6 = plt.subplot(3,2,5)

## load basemap
m = Basemap(projection='cyl',llcrnrlat=-60,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')

lon, lat = np.meshgrid(lons, lats)
x, y = m(lon, lat)

## draw contours
clevs = np.arange(-30.,80., 5)
cs = m.contourf(x,y,np.squeeze(tas),clevs,vmin=-31,norm=MidpointNormalize(midpoint=0),extend='both',cmap='RdBu')

## draw countries & coastlines
#m.drawcountries()
m.drawcoastlines()


plt.title('RCP2.6 2071-2099')

#%%
### HYPE

nc2 = NetCDFFile(file_path +'//'+ file_list[3])
tas2  = nc2.variables['qtot']
lons2 = nc2.variables['lon'][:]
lats2 = nc2.variables['lat'][:]
tas2 = np.average(tas2,axis=0)

ax3 = plt.subplot(3,2,2)


## load basemap
m = Basemap(projection='cyl',llcrnrlat=-60,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')


lon2, lat2 = np.meshgrid(lons2, lats2)
x, y = m(lon2, lat2)

## draw contours
clevs = np.arange(-30.,80., 5)
cs = m.contourf(x,y,np.squeeze(tas2),clevs,vmin=-31,norm=MidpointNormalize(midpoint=0),extend='both',cmap='RdBu')

## draw countries & coastlines
#m.drawcountries()
m.drawcoastlines()


# add title
plt.title('RCP6.0 2011-2040')


#%%
nc2 = NetCDFFile(file_path +'//'+ file_list[4])
tas2  = nc2.variables['qtot']
lons2 = nc2.variables['lon'][:]
lats2 = nc2.variables['lat'][:]
tas2 = np.average(tas2,axis=0)

ax4 = plt.subplot(3,2,4)


## load basemap
m = Basemap(projection='cyl',llcrnrlat=-60,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
#m = Basemap(projection='merc',llcrnrlat=-60,urcrnrlat=80,\
#            llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c' )
#m = Basemap(projection='cea',llcrnrlat=-90,urcrnrlat=90,\
#            llcrnrlon=-180,urcrnrlon=180,resolution='c')

lon2, lat2 = np.meshgrid(lons2, lats2)
x, y = m(lon2, lat2)


## draw contours
clevs = np.arange(-30.,80., 5)
cs = m.contourf(x,y,np.squeeze(tas2),clevs,vmin=-31,norm=MidpointNormalize(midpoint=0),extend='both',cmap='RdBu')

## draw countries & coastlines
#m.drawcountries()
m.drawcoastlines()

#parallels = np.arange(0.,90,15.)
#m.drawparallels(np.arange(-90., 91., 20.))
#
#meridians = np.arange(180.,360.,15.)
#m.drawmeridians(np.arange(-180., 180., 30.))

# add title
plt.title('RCP6.0 2041-2070')


#%%
## load data
nc2 = NetCDFFile(file_path +'//'+ file_list[5])
tas2  = nc2.variables['qtot']
lons2 = nc2.variables['lon'][:]
lats2 = nc2.variables['lat'][:]
tas2 = np.average(tas2,axis=0)

## make subplot
ax5 = plt.subplot(3,2,6)

## load basemap
m = Basemap(projection='cyl',llcrnrlat=-60,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')

lon2, lat2 = np.meshgrid(lons2, lats2)
x, y = m(lon2, lat2)


## draw contours
clevs = np.arange(-30.,80., 5)
cs = m.contourf(x,y,np.squeeze(tas2),clevs,vmin=-31,norm=MidpointNormalize(midpoint=0),extend='both',cmap='RdBu')

#cs = m.contourf(x,y,np.squeeze(tas2),clevs,vmin=-100.1, locator=ticker.LogLocator(),cmap='YlGnBu')

## draw countries & coastlines
#m.drawcountries()
m.drawcoastlines()


# add title
plt.title('RCP6.0 2071-2099')
#%%

## add colorbar hist
#cbaxes = fig.add_axes([0.93, 0.68, 0.01, 0.22]) 
#cb = plt.colorbar(cs, cax = cbaxes)
#cb.set_label('[mm/month]')

plt.suptitle('Relative difference runoff LPJmL', fontsize=16)


#plt.savefig(save_path +  '\\Runoff_C3S_ISIMIP_RCP8p5_complete')
plt.savefig(save_path +  '\\LPJmL_reldiff.jpeg', dpi=300, bbox_inches='tight')
