#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import urllib
from urllib import request
import requests 
import pandas as pd
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np
from numpy import loadtxt
import json
import pandas as df
from matplotlib import dates


# In[ ]:


def check_url(url):
    ret = requests.get(url)
    if ret.status_code == 200:
        return True
    else:
        return False


# In[ ]:


def plotEPS(year,month,pfuthreshold,instrument): 
    f = open(str(year)+'-'+str(month)+'--eps.json')
    data = json.load(f)

    # loading timestamps, pflux, satellite data in a usable manner 
    timestamps = data['time_tag']
    timestamps = np.array([datetime.strptime(date, '%Y-%m-%d %H:%M:%S') for date in timestamps])
    pflux = np.array(data['p3_flux_ic'], dtype=float)
    satellite_ids = np.array(data['satellite_id'])
    available_ids = []
    for sid in satellite_ids:
        if (sid not in available_ids):
            available_ids.append(sid)
    print("Available IDs:", available_ids)
    indexes = np.where( (satellite_ids == instrument) & (pflux > 0.0) )

    
    plt.figure(figsize = (15,10))
    pflux=np.log10(pflux[indexes])
    timestamps = timestamps[indexes]
    wherey=[]# indices of pflux where > pfuthreshold for event records 
    wherex=[]# indices of timestamp where > pfuthreshold for event records 
    counter=-1 

    for i in pflux:
        counter+=1
        if i*pfuthreshold > np.log10(pfuthreshold): 
            wherex.append(timestamps[counter])
            wherey.append(pflux[counter])
            
    start=[] # noting start of an event 
    end=[] # noting start of an event 
    sflux = [] # noting flux of start of an event to align with time 
    eflux = [ ]# noting flux of end of an event to align with time 
    event = False

    for i in range(0,len(wherey)):
        
        if (wherey[i] >= np.log10(pfuthreshold) # spike conditions, do nothing 
            and wherey[i-1] < np.log10(pfuthreshold) 
            and wherey[i+1] < np.log10(pfuthreshold)
            and wherey[i-2] < np.log10(pfuthreshold)
            and wherey[i+2] < np.log10(pfuthreshold)
            and wherex[i-1] - wherex[i] < timedelta(minutes=10) # spikes considered <10 minute differences in time
            and event == False):
            continue
            
        if (wherey[i] >= np.log10(pfuthreshold) # beginning of event, taking out spike possibilites 
            and wherey[i+1] > np.log10(pfuthreshold)
            and wherey[i-1] < np.log10(pfuthreshold)
            and wherey[i+2] > np.log10(pfuthreshold)
            and wherex[i+1] - wherex[i] < timedelta(minutes=30) # events with time difference at least 30 minutes from next time instance 
            and event == False):
            start.append(wherex[i])
            sflux.append(wherey[i])
            event = True
            
        if (wherey[i] >= np.log10(pfuthreshold) # end of event, saves after event over 
            and wherey[i+1] < np.log10(pfuthreshold)
            and event==True):
            end.append(wherex[i+1])
            eflux.append(wherey[i+1])
            event = False 

    # plotting event duration
    for i in range(0,len(start)):
        if len(end) != 0: 
            for i in range(0,len(start)):
                plt.axvline(start[i],color='green')
                plt.axvline(end[i],color='red')
        else:
            continue            
                    
    # table for event dates-peak-peakfluxval-fluence       
    with open(str(year)+'-'+str(month)+str(instrument)+"EPS.txt", "w") as savehere:
        for a,b,c,d in zip (start,sflux,end,eflux): 
            #start date, start flux, end date, end flux 
            try:
                ### aligning times with event flux instances as new arrays 
                beg = np.where(timestamps==a)[0][0] 
                fin = np.where(timestamps==c)[0][0]
                eventtimes=timestamps[beg:fin]  
                
                fstart = np.where(pflux==b)
                tstart = np.where(timestamps==a)
                
                startindex = np.intersect1d(fstart, tstart)[0]
                
                fend = np.where(pflux==d)
                tend = np.where(timestamps==c)
                      
                endindex = np.intersect1d(fend, tend)[0]
                
                eventtimes=timestamps[startindex:endindex]
                eventflux = pflux[startindex:endindex]
                peakval = np.max(eventflux)
                fluence = np.sum(eventflux) 
                ###
                
                # finding peak flux timing within new arrays 
                cc = -1
                for i,j in zip(pflux,timestamps):
                    cc+=1 
                    if i == peakval and i in eventflux and j in eventtimes:
                        peaktime=j
                
                if fluence > 0: 
#                     print('\nEvent Start: ', a)
#                     print('Event End: ',c)
#                     print('Peak Flux Time : ', peaktime)
#                     print('Peak Flux Value: ', str(peakval)[:4])
#                     print('Event Fluence: ', str(fluence)[:6])
                    savehere.write(str(a)+'\t')
                    savehere.write(str(c)+'\t')
                    savehere.write(str(peaktime)+'\t')
                    savehere.write(str(peakval)+'\t')
                    savehere.write(str(fluence)+'\n')

            except:
                print('Error')
                pass

    plt.plot(timestamps,pflux)
    hline = np.log10(pfuthreshold)
    plt.axhline(hline, label='>10 MeV')
    plt.plot([], [], color='green',label="Beginning of SPE")
    plt.plot([], [], color='red', label="Ending of SPE")
    plt.xticks(rotation=50)
    plt.ylabel('Logarithmic proton flux (pfu)')
    plt.ylim(-2,5)
    plt.title(str(year)+'-'+str(month)+' SPE Events >' + str(pfuthreshold) + ' MeV from ' + str(instrument))
    plt.legend(loc=0,frameon=True)
#     plt.savefig(str(year)+'-'+str(month)+'g6EPS.png')
    plt.show()
    
plotEPS(1990,8,10,'GOES-06')


# In[ ]:


def plotEPEAD(year,month,pfuthreshold,instrument): 
    f = open(str(year)+'-'+str(month)+'--epead.json')
    data = json.load(f)
    timestamps = data['time_tag']
    timestamps = np.array([datetime.strptime(date, '%Y-%m-%d %H:%M:%S') for date in timestamps])

    pflux = (np.array(data['ZPGT10E'], dtype=float) + np.array(data['ZPGT10W'], dtype=float)) / 2.0
    
    pflux_qe = np.array(data['ZPGT10E_QUAL_FLAG'], dtype=int)
    pflux_qw = np.array(data['ZPGT10W_QUAL_FLAG'], dtype=int)
    satellite_ids = np.array(data['satellite_id'])
    available_ids = []
    for sid in satellite_ids:
        if (sid not in available_ids):
            available_ids.append(sid)
    print("Available IDs:", available_ids)
    indexes= np.where( (satellite_ids == instrument) & (pflux_qe == 0) & (pflux_qw == 0) ) 

    plt.figure(figsize = (15,10))
    pflux=np.log10(pflux[indexes])
    timestamps = timestamps[indexes]
    wherey=[]# indices of pflux where > pfuthreshold 
    wherex=[]# indices of timestamp where > pfuthreshold 
    counter=-1 

    for i in pflux:
        counter+=1
        if i*pfuthreshold > np.log10(pfuthreshold): 
            wherex.append(timestamps[counter])
            wherey.append(pflux[counter])
            
    start=[]
    end=[]
    sflux = []
    eflux = []
    event = False

    for i in range(0,len(wherey)):
        
        if (wherey[i] >= np.log10(pfuthreshold) # spike conditions, do nothing 
            and wherey[i-1] < np.log10(pfuthreshold) 
            and wherey[i+1] < np.log10(pfuthreshold)
            and wherey[i-2] < np.log10(pfuthreshold)
            and wherey[i+2] < np.log10(pfuthreshold)
            and wherex[i-1] - wherex[i] < timedelta(minutes=10)
            and event == False):
            continue
            
        if (wherey[i] >= np.log10(pfuthreshold) # beginning of event 
            and wherey[i+1] > np.log10(pfuthreshold)
            and wherey[i-1] < np.log10(pfuthreshold)
            and wherey[i+2] > np.log10(pfuthreshold)
            and wherex[i+1] - wherex[i] < timedelta(minutes=30)
            and event == False):
            start.append(wherex[i])
            sflux.append(wherey[i])
            event = True
            
        if (wherey[i] >= np.log10(pfuthreshold) # end of event 
            and wherey[i+1] < np.log10(pfuthreshold)
            and event==True):
            end.append(wherex[i+1])
            eflux.append(wherey[i+1])
            event = False 
            
    for i in range(0,len(start)):
        if len(end) != 0: 
            for i in range(0,len(start)):
                plt.axvline(start[i],color='green')
                plt.axvline(end[i],color='red')
        else:
            continue 
            
    ### table for event dates-peak-peakfluxval-fluence       
    with open(str(year)+'-'+str(month)+str(instrument)+"EPEAD.txt", "w") as savehere:
        for a,b,c,d in zip (start,sflux,end,eflux): 
            #start date, start flux, end date, end flux 
            try:
                beg = np.where(timestamps==a)[0][0]
                fin = np.where(timestamps==c)[0][0]
                eventtimes=timestamps[beg:fin]  
                
                fstart = np.where(pflux==b)
                tstart = np.where(timestamps==a)
                
                startindex = np.intersect1d(fstart, tstart)[0]
                
                fend = np.where(pflux==d)
                tend = np.where(timestamps==c)
                      
                endindex = np.intersect1d(fend, tend)[0]
                
                eventtimes=timestamps[startindex:endindex]
                eventflux = pflux[startindex:endindex]
                peakval = np.max(eventflux)
                fluence = np.sum(eventflux)        
                
                cc = -1
                for i,j in zip(pflux,timestamps):
                    cc+=1 
                    if i == peakval and i in eventflux and j in eventtimes:
                        peaktime=j
                
                if fluence > 0: 
                    print('\nEvent Start: ', a)
                    print('Event End: ',c)
                    print('Peak Flux Time : ', peaktime)
                    print('Peak Flux Value: ', str(peakval)[:4])
                    print('Event Fluence: ', str(fluence)[:6])
                    savehere.write(str(a)+'\t')
                    savehere.write(str(c)+'\t')
                    savehere.write(str(peaktime)+'\t')
                    savehere.write(str(peakval)+'\t')
                    savehere.write(str(fluence)+'\n')
            except:
                print('Error')
                pass

    plt.plot(timestamps,pflux)
    hline = np.log10(pfuthreshold)
    plt.axhline(hline, label='>10 MeV')
    plt.plot([], [], color='green',label="Beginning of SPE")
    plt.plot([], [], color='red', label="Ending of SPE")
    plt.xticks(rotation=50)
    plt.ylabel('Logarithmic proton flux (pfu)')
    plt.ylim(-2,5)
    plt.title(str(year)+'-'+str(month)+' SPE Events >' + str(pfuthreshold) + ' MeV from ' + instrument)
    plt.legend(loc=0,frameon=True)
#     plt.savefig(str(year)+'-'+str(month)+'g14EPEAD.png')
    plt.show()
    
plotEPEAD(2017,9,10,'GOES-13')


# In[ ]:


# Automating Generation : 

# command line list empty txt files : 
# find . -name "*.txt" -size 0k // add prefix of different GOES satellites to find 'empty dates'  and save terminal output as *empty.txt to use with emptydatadates.ipynb

# command line delete empty txt files : 
# find . -type f -size 0 -exec rm -f '{}' +

# command line delete empty plot files : 
# find . -name "*.png" -type 'f' -size -19k -delete


# In[ ]:


# looping through all GOES instruments 
ins = ['GOES-06','GOES-07','GOES-08','GOES-09','GOES-10','GOES-11','GOES-12']


# EPS : 

# In[ ]:


for year in range(1986,2012):
    for month in range (1,13):
        for i in ins: 
            try: 
                plotEPS(year,month,10,i)
            except:
                print(str(year) + '-' + str(month) + ' Data not available')


# EPEAD:

# In[ ]:


# looping through all GOES instruments 
ins = ['GOES-13','GOES-14','GOES-15']


# In[ ]:


for year in range(2010,2021):
    for month in range (1,13):
        for i in ins: 
            try: 
                plotEPEAD(year,month,10,i)
            except:
                print(str(year) + '-' + str(month) + ' Data not available')

