# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 13:12:09 2016

@author: xiaojian
"""
# -*- coding: utf-8 -*-
import numpy as np
#from pydap.client import open_url
import matplotlib.pyplot as plt
from SeaHorseLib import *
from datetime import *
from scipy import interpolate
import sys
from SeaHorseTide import *
import shutil

SOURCEDIR='data3/'
DESTINDIR='data4/'

FList = np.genfromtxt(SOURCEDIR+'FList.csv',dtype=None,names=['FNs'],delimiter=',')
FNs=list(FList['FNs'])

filename='FList.csv'
FN3='data4/'+filename
csvfile1 = file(FN3, 'wb')

k=0
while k in range(len(FNs)):
    try:
        P = np.genfromtxt('kpic.txt',dtype=None,names=['INDEX'],delimiter=',')
        kpic=P['INDEX']
        k=kpic
    except IOError:
        k=0
    
    
    FN=FNs[k]
    #ID_19965381.npz
    FN1=SOURCEDIR+FN
    print k, FN1
    Z=np.load(FN1) 
    tdh=Z['tdh'];londh=Z['londh'];latdh=Z['latdh'];
    udh=Z['udh'];vdh=Z['vdh'];
    udm=Z['udm'];vdm=Z['vdm'];
    tgap=Z['tgap'];flag=Z['flag'];
    udti=Z['udti'];vdti=Z['vdti'];
    
    Z.close()
    
            
    flag=udh*0.+1. # 1 inside GOM3 NaN outside
    
    flag_gap=tdh*0.+1.
    i=np.argwhere(tgap>12./24.) # exceeds 12h
    flag_gap[i]=np.NaN
    
    flag=flag*flag_gap
    
    
    v=np.sqrt(udh*udh+vdh*vdh)    
    i=np.argwhere(v>2.7) # 2 m/s
    flag_v=tdh*0.+1.
    flag_v[i]=np.NaN
    
    plt.close('all')
    plt.figure(1)
    plt.plot(londh,latdh,'r.-')
    plt.plot(londh*flag,latdh*flag,'m.-')
    flag=flag*flag_v
    plt.plot(londh*flag,latdh*flag,'b.-')
    plt.title(FN1)
    plt.show()
    
    i=np.argwhere(flag == 1).flatten()
    if len(i)>0:
        FN2=DESTINDIR+FN
        print FN2
        #        np.savez(FN2,tdh=tdh,londh=londh,latdh=latdh,udh=udh,vdh=vdh,umoh=umoh,vmoh=vmoh)
        shutil.copyfile(FN1,FN2)
            
        writer1=csv.writer(csvfile1)
        a=FN2[6:]
        print 'a',a
        b=[]
        b.append(a)
        print 'b',b
        writer1.writerow(b)
    
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
    
    k=k+1    
    f=open('kpic.txt','w')
    f.write(str(k))
    f.close()    
    
