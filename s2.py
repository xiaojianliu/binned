# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 13:12:09 2016

@author: xiaojian
"""
import csv
import numpy as np
import matplotlib.pyplot as plt
from SeaHorseLib import *
from datetime import *
from scipy import interpolate
import sys
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

SOURCEDIR='data2/'
DESTINDIR='data3/'

FList = np.genfromtxt(SOURCEDIR+'FList.csv',dtype=None,names=['FNs'],delimiter=',')
FNs=list(FList['FNs'])
kdr=0
filename='FList.csv'
FN3='data3/'+filename
csvfile1 = file(FN3, 'wb')

while kdr in range(len(FNs)):
#while kdr in range(0,284):
#while kdr in range(9,10):

    FN=FNs[kdr]
    FN1=SOURCEDIR+FN
    Z=np.load(FN1)
    tdh=Z['tdh'];londh=Z['londh'];latdh=Z['latdh'];
    udh=Z['udh'];vdh=Z['vdh'];
    Z.close()

    print kdr, FN

# need the raw drifter tracks to find gaps    
    FND=FN[0:-4]
    FND='data1/'+FND+'.csv'
#driftfvcom_data/ID_19965381_dm.npz
    D = np.genfromtxt(FND,dtype=None,names=['ID','TIME','LON','LAT','DEPTH'],delimiter=',',skip_header=1)
    td=[]
    for a in np.arange(len(D['TIME'])):
        td.append(RataDie1(int(D['TIME'][a][0:4]),int(D['TIME'][a][5:7]),int(D['TIME'][a][8:10]),int(D['TIME'][a][11:13]),int(D['TIME'][a][14:16]),int(D['TIME'][a][17:19])))    
    
    latd=np.array(D['LAT'])
    lond=np.array(D['LON'])
    
    #start and end times close to a whole hour
    t1=np.ceil(np.min(td)*24.)/24.
    t2=np.floor(np.max(td)*24.)/24.

# remove tidal signal
    udm=sh_rmtide(udh,ends=np.NaN)
    vdm=sh_rmtide(vdh,ends=np.NaN)
    udti=udh-udm
    vdti=vdh-vdm
    '''
    umom=sh_rmtide(umoh,ends=np.NaN)
    vmom=sh_rmtide(vmoh,ends=np.NaN)
    umoti=umoh-umom
    vmoti=vmoh-vmom

# it does not make much sense to remove tides from the wind stress
# but it might be helpful if comparing with the rmtide filtered velocities
    uwsm=sh_rmtide(uwsh,ends=np.NaN)
    vwsm=sh_rmtide(vwsh,ends=np.NaN)
    '''

# find gaps
    tgap=tdh*0.+1./24.
    
    for i in range(len(tdh)):
        ti=tdh[i]
        i1=max(np.argwhere(td<=ti))
        i2=min(np.argwhere(td>ti))
       
        tgap[i]=td[i2]-td[i1]
        
    flag_gap=tdh*0.+1.
    for i in range(len(tdh)):
        if (tgap[i] > 9./24.): 
            flag_gap[i]=np.NaN
            
    flag=udm*0.+1. # 1 inside GOM3 NaN outside
    flag=flag*flag_gap

    plt.figure(1)
    plt.plot(londh*flag,latdh*flag,'r.-')
    
 
    """   
    plt.figure()
    plt.plot(tdhgap,udh,'b.-',tdhgap,umoh,'r.-',tdhgap,udm,'k-',tdhgap,umom,'g-')
    plt.title('U '+FN0)
    plt.legend(['Drift','FVCOM','Drift rmtide','FVCOM rmtide'])
    plt.show()

    plt.figure()
    plt.plot(tdhgap,vdh,'b.-',tdhgap,vmoh,'r.-',tdhgap,vdm,'k-',tdhgap,vmom,'g-')
    plt.title('V '+FN0)
    plt.legend(['Drift','FVCOM','Drift rmtide','FVCOM rmtide'])
    plt.show()
    """    
    FN2=DESTINDIR+FN
    np.savez(FN2,tdh=tdh,tgap=tgap,flag=flag,londh=londh,latdh=latdh,udh=udh,vdh=vdh,udm=udm,vdm=vdm,udti=udti,vdti=vdti)

    kdr=kdr+1
    
    writer1=csv.writer(csvfile1)
    a=FN2[6:]
    print 'a',a
    b=[]
    b.append(a)
    print 'b',b
    writer1.writerow(b)
    
plt.show()
