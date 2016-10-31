# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 14:22:58 2016

@author: xiaojian
"""
#http://www.ngdc.noaa.gov/mgg/coast/
# coast line data extractor

import numpy as np
#from pydap.client import open_url
import matplotlib.pyplot as plt
#from SeaHorseLib import *
#from datetime import *
#from scipy import interpolate
#import sys
#from SeaHorseTide import *
#import shutil
import matplotlib.mlab as mlab
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import  Axes3D
from matplotlib import  cm
from matplotlib.ticker import LinearLocator,FormatStrFormatter
import netCDF4
###################HARDCODE########################3
scale1=1.
scale2=4.
def sh_bindata(x, y, z, xbins, ybins):
    """
    Bin irregularly spaced data on a rectangular grid.

    """
    ix=np.digitize(x,xbins)
    iy=np.digitize(y,ybins)
    xb=0.5*(xbins[:-1]+xbins[1:]) # bin x centers
    yb=0.5*(ybins[:-1]+ybins[1:]) # bin y centers
    zb_mean=np.empty((len(xbins)-1,len(ybins)-1),dtype=z.dtype)
    zb_median=np.empty((len(xbins)-1,len(ybins)-1),dtype=z.dtype)
    zb_std=np.empty((len(xbins)-1,len(ybins)-1),dtype=z.dtype)
    zb_num=np.zeros((len(xbins)-1,len(ybins)-1),dtype=int)    
    for iix in range(1,len(xbins)):
        for iiy in range(1,len(ybins)):
#            k=np.where((ix==iix) and (iy==iiy)) # wrong syntax
            k,=np.where((ix==iix) & (iy==iiy))
            zb_mean[iix-1,iiy-1]=np.mean(z[k])
            zb_median[iix-1,iiy-1]=np.median(z[k])
            zb_std[iix-1,iiy-1]=np.std(z[k])
            zb_num[iix-1,iiy-1]=len(z[k])
            
    return xb,yb,zb_mean,zb_median,zb_std,zb_num
def get_nc_data(url, *args):
    '''
    get specific dataset from url

    *args: dataset name, composed by strings
    ----------------------------------------
    example:
        url = 'http://www.nefsc.noaa.gov/drifter/drift_tcs_2013_1.dat'
        data = get_url_data(url, 'u', 'v')
    '''
    nc = netCDF4.Dataset(url)
    data = {}
    for arg in args:
        try:
            data[arg] = nc.variables[arg]
        except (IndexError, NameError, KeyError):
            print 'Dataset {0} is not found'.format(arg)
    #print data
    return data
"""
from netCDF4 import Dataset

# read in etopo5 topography/bathymetry.
url = 'http://ferret.pmel.noaa.gov/thredds/dodsC/data/PMEL/etopo5.nc'
etopodata = Dataset(url)

topoin = etopodata.variables['ROSE'][:]
lons = etopodata.variables['ETOPO05_X'][:]
lats = etopodata.variables['ETOPO05_Y'][:]
# shift data so lons go from -180 to 180 instead of 20 to 380.
topoin,lons = shiftgrid(180.,topoin,lons,start=False)
"""



"""
BATHY=np.genfromtxt('necscoast_noaa.dat',dtype=None,names=['coast_lon', 'coast_lat'])
coast_lon=BATHY['coast_lon']
coast_lat=BATHY['coast_lat']
"""

#BATHY=np.genfromtxt('coastlineNE.dat',names=['coast_lon', 'coast_lat'],dtype=None,comments='>')
#coast_lon=BATHY['coast_lon']
#coast_lat=BATHY['coast_lat']


# www.ngdc.noaa.gov
# world vector shoreline ascii
FNCL='necscoast_worldvec.dat'
# lon lat pairs
# segments separated by nans
"""
nan nan
-77.953942	34.000067
-77.953949	34.000000
nan nan
-77.941035	34.000067
-77.939568	34.001241
-77.939275	34.002121
-77.938688	34.003001
-77.938688	34.003881
"""
CL=np.genfromtxt(FNCL,names=['lon','lat'])


FN='binned.npz'
#FN='binned_model.npz'
Z=np.load(FN) 
xb=Z['xb']
yb=Z['yb']
ub_mean=Z['ub_mean']
ub_median=Z['ub_median']
ub_std=Z['ub_std']
ub_num=Z['ub_num']
vb_mean=Z['vb_mean']
vb_median=Z['vb_median']
vb_std=Z['vb_std']
vb_num=Z['vb_num']
Z.close()
for a in np.arange(len(ub_num)):
    for b in np.arange(len(ub_num[0])):
        if ub_num[a][b]<10:
            ub_num[a][b]=0
            ub_mean[a][b]=np.nan
            vb_mean[a][b]=np.nan
duf=ub_mean
dvf=vb_mean
aa=np.empty((len(ub_mean),len(ub_mean[0])),dtype=ub_mean.dtype)
for a in np.arange(len(ub_mean)):
    for b in np.arange(len(ub_mean[0])):
        aa[a,b]=np.sqrt(ub_mean[a,b]*ub_mean[a,b]+vb_mean[a,b]*vb_mean[a,b])
        
dv=aa
#cmap = matplotlib.cm.jet
#cmap.set_bad('w',1.)
xxb,yyb = np.meshgrid(xb, yb)
cc=np.arange(-1.5,1.500001,0.05)
#cc=np.array([-1., -.75, -.5, -.25, -0.2, -.15, -.1, -0.05, 0., 0.05, .1, .15, .2, .25, .5, .75, 1.])

plt.figure()
ub = np.ma.array(ub_mean, mask=np.isnan(ub_mean))
vb = np.ma.array(vb_mean, mask=np.isnan(vb_mean))
Q=plt.quiver(xxb,yyb,ub.T,vb.T,scale=scale1)
#Q=plt.arrow(xxb,yyb,ub.T,vb.T,head_width=.05)
qk=plt.quiverkey(Q,0.8,0.1,0.1, r'$10cm/s$',labelpos='E',linewidth=0.005, fontproperties={'weight': 'bold','size':'larger','family':'monospace'})#0.8 0.1 0.1

plt.title(FN[:-4]+', mean')
#plt.title('binned_drifter_mean')
#plt.plot(coast_lon,coast_lat,'b.')
plt.plot(CL['lon'],CL['lat'])
plt.axis([-70.75,-70,41.63,42.12])
plt.savefig('binned_drifter_mean')
plt.show()

plt.figure()
ub = np.ma.array(ub_std, mask=np.isnan(ub_mean))
vb = np.ma.array(vb_std, mask=np.isnan(vb_mean))
Q=plt.quiver(xxb,yyb,ub.T,vb.T,scale=scale2)
qk=plt.quiverkey(Q,0.8,0.1,0.4, r'$40cm/s$',labelpos='W', linewidth=0.005,fontproperties={'weight': 'bold','size':'larger','family':'monospace'})

plt.title(FN[:-4]+', std')

#plt.plot(coast_lon,coast_lat,'b.')
plt.plot(CL['lon'],CL['lat'])
plt.axis([-70.75,-70,41.63,42.12])
plt.savefig('binned_drifter_std')
plt.show()

ubn = np.ma.array(ub_num, mask=np.isnan(ub_num))
vbn = np.ma.array(vb_num, mask=np.isnan(vb_num))
plt.figure()
#plt.plot(xxb+0.015,yyb+0.015,'r-')
#plt.plot((xxb+0.015).T,(yyb+0.015).T,'r-')
#Q=plt.quiver(xxb,yyb,ub.T,vb.T,scale=5.)
#qk=plt.quiverkey(Q,0.8,0.1,0.5, r'$50cm/s$', fontproperties={'weight': 'bold'})
for a in np.arange(len(xxb[0])):
    for b in np.arange(len(yyb)):
        plt.text(xxb[0][a],yyb[b][0],ub_num[a][b])
plt.plot(CL['lon'],CL['lat'])
plt.axis([-70.75,-70,41.63,42.12])
plt.title('binned_drifter_num')
plt.savefig('binned_drifter_num')
plt.show()
  
plt.figure()
plt.title('sea_depth')
   
data = np.genfromtxt('sea.csv',dtype=None,names=['x','y','h'],delimiter=',')    
x=[]
y=[]
h=[]
x=data['x']
y=data['y']
h=data['h']
xi = np.arange(-70.75,-70.00,0.03)
yi = np.arange(41.63,42.12,0.03)
xb,yb,hb_mean,hb_median,hb_std,hb_num = sh_bindata(x, y, h, xi, yi)
xxxb,yyyb = np.meshgrid(xb, yb)
CS=plt.contour(xxxb, yyyb, -abs(hb_mean.T))
plt.clabel(CS, inline=1, fontsize=10)
plt.plot(CL['lon'],CL['lat'])
plt.axis([-70.75,-70,41.63,42.12])
plt.savefig('sea_depth')   
plt.show()
