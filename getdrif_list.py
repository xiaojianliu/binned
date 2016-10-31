# code to get drifter ids given gbox and time
# JiM and Xiaojian in Oct 2016
from datetime import datetime as dt
import pytz
import pandas as pd
import numpy as np
import csv
input_time=[dt(1988,10,11,0,0,0,0,pytz.UTC),dt(2017,1,1,0,0,0,0,pytz.UTC)] # start time and end time
gbox=[-70.,-70.75, 42.12, 41.63]
depth=1.0


def getobs_drift_byrange(gbox,input_time,depth):
    """
   Function written by Huanxin and used in "getdrifter_erddap.py"
   get data from url, return id, latitude,longitude, and times
   gbox includes 4 values, maxlon, minlon,maxlat,minlat, like:  [-69.0,-73.0,41.0,40.82]
   input_time can either contain two values: start_time & end_time OR one  value:interval_days
   and they should be timezone aware
   example: input_time=[dt(2012,1,1,0,0,0,0,pytz.UTC),dt(2012,2,1,0,0,0,0,pytz.UTC)]
   """
    lon_max=gbox[0];lon_min=gbox[1];lat_max=gbox[2];lat_min=gbox[3]
    mintime=input_time[0].strftime('%Y-%m-%d'+'T'+'%H:%M:%S'+'Z')  # change time format
    maxtime=input_time[1].strftime('%Y-%m-%d'+'T'+'%H:%M:%S'+'Z')
    # open url to get data
    url='http://comet.nefsc.noaa.gov:8080/erddap/tabledap/drifters.csv?id,depth&time>='\
    +str(mintime)+'&time<='+str(maxtime)+'&latitude>='\
    +str(lat_min)+'&latitude<='+str(lat_max)+'&longitude>='+str(lon_min)+'&longitude<='+str(lon_max)+'&depth<='+str(depth)+'&orderBy("id")'
    df=pd.read_csv(url,skiprows=[1])
    ids=df.id.values
    deep=df.depth.values
    idss=[]
    for a in np.arange(len(deep)):        
        if abs(deep[a])==1:
            idss.append(ids[a])
    news_ids = []
    for id in idss:
        if id not in news_ids:
            news_ids.append(id)
    return news_ids
ids=getobs_drift_byrange(gbox,input_time,depth)
filename='FList.csv'
csvfile = file(filename, 'wb')
writer = csv.writer(csvfile)
writer.writerow(['ids'])
for a in np.arange(len(ids)):
    writer.writerow([str(ids[a])+'.csv'])
csvfile.close()

lon_max=gbox[0];lon_min=gbox[1];lat_max=gbox[2];lat_min=gbox[3]
mintime=input_time[0].strftime('%Y-%m-%d'+'T'+'%H:%M:%S'+'Z')  # change time format
maxtime=input_time[1].strftime('%Y-%m-%d'+'T'+'%H:%M:%S'+'Z') 
url='http://comet.nefsc.noaa.gov:8080/erddap/tabledap/drifters.csv?id,time,longitude,latitude,depth&time>='\
    +str(mintime)+'&time<='+str(maxtime)+'&latitude>='\
    +str(lat_min)+'&latitude<='+str(lat_max)+'&longitude>='+str(lon_min)+'&longitude<='+str(lon_max)+'&depth<='+str(depth)+'&orderBy("id,time")'
df=pd.read_csv(url,skiprows=[1])
print 1111111111111111
b=0
while b<len(ids):
    filename=str(ids[b])+'.csv'
    print 2222222222
    csvfile = file(filename, 'wb')
    writer = csv.writer(csvfile)
    print 333333333
    writer.writerow(['id', 'time', 'lon','lat','depth'])
    #writer.writerow(['id','time','lon','lat','depth'])
    print 444444444
    for a in np.arange(len(df['id'])):
        
        if df['id'][a]==ids[b]:
            
            
            print 'aaaaaaaa'
            writer.writerow([df['id'][a],df['time'][a],df['longitude'][a],df['latitude'][a],-abs(df['depth'][a])])
        else:
            continue
    csvfile.close()
    b=b+1





















