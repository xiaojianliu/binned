# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 11:41:04 2016
s1_mk1h_hour.py
compare drifter dataset velocity with fvcom hincast
following the drifter trajectory

step 1
fvcom data from local hourly files

To do: 
1. spacial interpolation of fvcom velocity data
currently, the nearest neighbor
fixed to polygonal baricentric coordinate interpolation

before running
reset kdhour.txt file

@author: xiaojian
"""

import numpy as np
from SeaHorseLib import *
from datetime import *
import sys
import matplotlib.pyplot as plt
from SeaHorseTide import *
def RataDie1(yr,mo,da,hr,mi,se):
    """

RD = RataDie(yr,mo=1,da=1,hr=0,mi=0,se=0)
RD = RataDie(yr,mo=1,da=1)

returns the serial day number in the (proleptic) Gregorian calendar
or elapsed time in days since 0001-01-00.

Vitalii Sheremet, SeaHorse Project, 2008-2013.
"""
#
#    yr+=(mo-1)//12;mo=(mo-1)%12+1; # this extends mo values beyond the formal range 1-12
    RD=367*yr-(7*(yr+((mo+9)//12))//4)-(3*(((yr+(mo-9)//7)//100)+1)//4)+(275*mo//9)+da-396;
    RD=RD+(((se/60.+mi)/60.)+hr)/24.
    return RD    
def RataDie(yr,mo,da,hr):
    """

RD = RataDie(yr,mo=1,da=1,hr=0,mi=0,se=0)
RD = RataDie(yr,mo=1,da=1)

returns the serial day number in the (proleptic) Gregorian calendar
or elapsed time in days since 0001-01-00.

Vitalii Sheremet, SeaHorse Project, 2008-2013.
"""
#
#    yr+=(mo-1)//12;mo=(mo-1)%12+1; # this extends mo values beyond the formal range 1-12
    RD=367*yr-(7*(yr+((mo+9)//12))//4)-(3*(((yr+(mo-9)//7)//100)+1)//4)+(275*mo//9)+da-396;
    RD=RD*24+hr
    return RD    
FList = np.genfromtxt('data1/FList.csv',dtype=None,names=['FNs'],delimiter=',',skip_header=1)
FNs=list(FList['FNs'])
print FNs
kdr=0
try:
    P = np.genfromtxt('kdrhour.txt',dtype=None,names=['INDEX'],delimiter=',')
    kdr=P['INDEX']
#    kdr = np.load('kdr.npy')
    #print 'continue with drifter '+str(kdr)
except IOError:
    kdr=0
    #print 'start with drifter '+str(kdr)
filename='FList.csv'
FN3='data2/'+filename
csvfile1 = file(FN3, 'wb')

while kdr in range(len(FNs)):
     
    FN='data1/'+FNs[kdr]
    """
    # ID - drifter ID
    # TIME - RataDie, serial day number, days since 0001-00-00 00:00:00, GMT
    # TIME_GMT - timestamp, yyyy-mm-dd, GMT
    # YRDAY0_GMT - yearday, zero based indexing, days since Jan-01 00:00:00, GMT
    # LON_DD - longitude, degrees east
    # LAT_DD - latitude, degrees north
    # TEMP - temperature, degrees Celsius
    # DEPTH_I - depth of instrument drogue, meters
    #Vars: ID, TimeRD, TIME_GMT, YRDAY0_GMT, LON_DD, LAT_DD, TEMP, DEPTH_I
    19954471,728387.5687,1995-04-04,93.5687,-67.39,44.167,9.99,40.0
    19954471,728387.9056,1995-04-04,93.9056,-67.429,44.213,9.99,40.0
    19954471,728387.975,1995-04-04,93.975,-67.428,44.205,9.99,40.0
    19954471,728388.0451,1995-04-05,94.0451,-67.431,44.202,9.99,40.0
    19954471,728388.2576,1995-04-05,94.2576,-67.44,44.224,9.99,40.0
    """
    d1=np.genfromtxt(FN,dtype=None,names=['ID','TIME','LON','LAT','DEPTH'],delimiter=',')
    if len(d1['ID'])<3:
        kdr=kdr+1
        f=open('kdrhour.txt','w')
        f.write(str(kdr))
        f.close()
        continue
    D = np.genfromtxt(FN,dtype=None,names=['ID','TIME','LON','LAT','DEPTH'],delimiter=',',skip_header=1)
    td=[]
    td1=[]
    print 11111111111
    
  
    for a in np.arange(len(D['TIME'])):
            #print 2222222222
        td1.append(datetime.strptime(D['TIME'][a],'%Y-%m-%dT%H:%M:%SZ'))
    for a in np.arange(len(D['TIME'])):
        td.append(RataDie1(int(D['TIME'][a][0:4]),int(D['TIME'][a][5:7]),int(D['TIME'][a][8:10]),int(D['TIME'][a][11:13]),int(D['TIME'][a][14:16]),int(D['TIME'][a][17:19])))
    latd=np.array(D['LAT'])
    lond=np.array(D['LON'])
    depd=np.median(np.array(D['DEPTH']))
    '''
    start_time=td1[0]
    if start_time.minute==0 and start_time.second==0:
        a=start_time.year,start_time.month,start_time.day,start_time.hour
    else:
        a=start_time.year,start_time.month,start_time.day,start_time.hour+1
    a=str(a[0])+'-'+str(a[1])+'-'+str(a[2])+'T'+str(a[3])+':'+'00'+':'+'00'+'Z'
    a=datetime.strptime(a,'%Y-%m-%dT%H:%M:%SZ')
    end_time=td1[-1]
    t=(end_time-a)
    h=t.days*24+t.seconds/3600
    start_time=a
    latd=np.array(D['LAT'])
    lond=np.array(D['LON'])
    depd=np.median(np.array(D['DEPTH']))
    print 'depd',depd
    #start and end times close to a whole hour
    t1=RataDie(start_time.year,start_time.month,start_time.day,start_time.hour)
    t2=t1+h
    tdh=np.arange(t1,t2,1./24.)
    '''    
        #start and end times close to a whole hour
    t1=np.ceil(np.min(td)*24.)/24.
    t2=np.floor(np.max(td)*24.)/24.
    
    # available range of FVCOM GOM3 30yr hindcast
    # record must be long enogh for interpolation to work   
        #if (t1>RataDie(1978,1,1)) and len(td)>3:
        
            # interpolate drifter trajectory to hourly 
            #tdh=np.arange(t1,t2,1./24.)
    tdh=np.arange(t1,t2,1./24.)
            #latdh=np.interp(tdh,td,latd)
            #londh=np.interp(tdh,td,lond)
    
            # optionally cubic interpolation: scipy interpolate interp1d
            #fip1=interpolate.interp1d(td,latd,kind='cubic',bounds_error=False) #returns function not array
            #latdh=fip1(tdh)
            #fip2=interpolate.interp1d(td,lond,kind='cubic',bounds_error=False) #returns function not array
            #londh=fip2(tdh)         
            # the above algorithm gets very slow for long arrays
            # and fails for many long tracks because of lack of memory
            # it must have used matrix inversion to calculate spline over all data points.
            
            # I use cubic Hermite polynomial spline, SeaHorseTide.sh_interp3.
            # It is a local algorithm only using neighbor data
    latdh=sh_interp3(tdh,td,latd)
    londh=sh_interp3(tdh,td,lond)
            
            
            
            
    plt.figure()
    plt.plot(D['LON'],D['LAT'],'b.-')
    plt.plot(londh,latdh,'r.')
    plt.show()
            
            
    ND=len(tdh)
    print 'points in track', ND 
    '''
    kvdh=np.zeros(ND,dtype=int)
    u0=np.zeros(ND,dtype=float)   
    v0=np.zeros(ND,dtype=float)   
    uws=np.zeros(ND,dtype=float)   
    vws=np.zeros(ND,dtype=float)   
    '''
            
    udh=np.zeros(ND,dtype=float)   
    vdh=np.zeros(ND,dtype=float)
    Coef=111111./86400. # deg/day -> m/s
            #print Coef
    for i in range(1,ND-1):
        udh[i]=(londh[i+1]-londh[i-1])/(tdh[i+1]-tdh[i-1])*Coef*np.cos(latdh[i]*np.pi/180.)
        vdh[i]=(latdh[i+1]-latdh[i-1])/(tdh[i+1]-tdh[i-1])*Coef
            
            
    FN1=FNs[kdr]
    FN2=FN1[0:-4]
    FN2='data2/'+FN2+'.npz'
    print FN2
    np.savez(FN2,tdh=tdh,londh=londh,latdh=latdh,udh=udh,vdh=vdh)
            
    kdr=kdr+1 # next drifter
    #    np.save('kdr.npy',kdr)
    f=open('kdrhour.txt','w')
    f.write(str(kdr))
    f.close()
    
    writer1 = csv.writer(csvfile1)
    a=FN2[6:]
    b=[]
    b.append(a)
    writer1.writerow(b)
 #csvfile1.close()