#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[2]:


def check_url(url):
    ret = requests.get(url)
    if ret.status_code == 200:
        return True
    else:
        return False


# In[3]:


# EPS:
# need to run only the first time: 

# base_url = ' https://sun.njit.edu/devSEP_Slava/SEP_API/apiget_obs_EPS.php?start_time='
# end_url='&end_time='

# for i in range(1986,2012): #years 1986-2011 
#     for j in range (1,13):#allmonths
#         if j < 10:
#             for k in range(2,32):
#                 try:
#                     temp = base_url+str(i)+'-0'+str(j)+'-0'+str(1)+'T00:00:00'+end_url+str(i)+'-0'+str(j)+'-'+f"{k:02d}"+'T00:00:00'
#                     print(check_url(temp))
#                     if check_url(temp) == True: 
#                         site = urllib.request.urlopen(temp).read()
#                         if len(site) != 2: 
#                             print('here',k,temp)
#                             r = requests.get(url = temp)
#                             data = r.json()
#                             with open(str(i)+'-'+str(j)+'--eps.json', 'w') as f:
#                                 json.dump(data, f)
#                 except:
#                     print('Server Error')
#                     continue 
                
#         if j >=10:
#             for k in range(2,32):
#                 try:
#                     temp = base_url+str(i)+'-'+str(j)+'-0'+str(1)+'T00:00:00'+end_url+str(i)+'-'+str(j)+'-'+f"{k:02d}"+'T00:00:00'
#                     print(temp)
#                     if check_url(temp) == True: 
#                         site = urllib.request.urlopen(temp).read()
#                         if len(site) != 2: 
#                             print('here',k,temp)
#                             r = requests.get(url = temp)
#                             data = r.json()
#                             with open(str(i)+'-'+str(j)+'--eps.json', 'w') as f:
#                                 json.dump(data, f)
#                 except:
#                     print('Server Error')
#                     continue 


# In[ ]:


# EPEAD : 
# ONLY NEED TO RUN ONCE 

# base_url = 'https://sun.njit.edu/devSEP_Slava/SEP_API/apiget_obs_EPEAD.php?start_time='
# end_url='&end_time='

# for i in range(2010,2021): #years 1986-2020
#     for j in range (1,13):#allmonths
#         if j < 10:
#             for k in range(2,32):
#                 try:
#                     temp = base_url+str(i)+'-0'+str(j)+'-0'+str(1)+'T00:00:00'+end_url+str(i)+'-0'+str(j)+'-'+f"{k:02d}"+'T00:00:00'
#                     if check_url(temp) == True: 
#                         site = urllib.request.urlopen(temp).read()
#                         if len(site) != 2: 
#                             print('here',k,temp)
#                             r = requests.get(url = temp)
#                             data = r.json()
#                             with open(str(i)+'-'+str(j)+'--epead.json', 'w') as f:
#                                 json.dump(data, f)
#                 except:
#                     print('Server Error')
#                     continue 
#         if j >= 10: 
#             for k in range(2,32):
#                 try:
#                     temp = base_url+str(i)+'-'+str(j)+'-0'+str(1)+'T00:00:00'+end_url+str(i)+'-'+str(j)+'-'+f"{k:02d}"+'T00:00:00'
#                     print(temp)
#                     if check_url(temp) == True: 
#                         site = urllib.request.urlopen(temp).read()
#                         if len(site) != 2: 
#                             print('here',k,temp)
#                             r = requests.get(url = temp)
#                             data = r.json()
#                             with open(str(i)+'-'+str(j)+'--epead.json', 'w') as f:
#                                 json.dump(data, f)
#                 except:
#                     print('Server Error')
#                     continue 


# In[5]:


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
#     pflux=fixer(pflux,m=0) # reducing spikes 
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
    eflux = []# noting flux of end of an event to align with time 
    event = False # starting with no event happening
    #     ind = np.where(pflux==np.max(pflux))[0][0]
#     print(pflux[ind],pflux[ind-1],pflux[ind+1])
#     print(timestamps[ind])
#     print(pflux[ind]>=np.log10(pfuthreshold)
#              ,pflux[ind-1] <= np.log10(pfuthreshold)
#              ,pflux[ind+1] <= np.log10(pfuthreshold))
    for i in range(0,len(pflux)):
            
        if (pflux[i] >= np.log10(pfuthreshold) # beginning of event, taking out spike possibilites 
            and pflux[i+1] >= np.log10(pfuthreshold)
            and pflux[i-1] <= np.log10(pfuthreshold)
            and pflux[i+2] >= np.log10(pfuthreshold)
            and pflux[i+3] >= np.log10(pfuthreshold)
            and timestamps[i+1] - timestamps[i] < timedelta(minutes=30) # events with time difference at least 30 minutes from next time instance 
            and event == False):
            print('START: i-1,i,i+1,i+2\n',pflux[i-1],pflux[i],pflux[i+1])
            print('at',timestamps[i])
            start.append(timestamps[i])
            sflux.append(pflux[i])
            event = True

        if (pflux[i] >= np.log10(pfuthreshold) # end of event, saves after event over 
            and pflux[i+1] <= np.log10(pfuthreshold)
            and event==True):
            end.append(timestamps[i+1])
            eflux.append(pflux[i+1])
            event = False 
            
        if (event == False # spike condition, replace flux with previous flux
            and pflux[i]>=np.log10(pfuthreshold)
            and pflux[i-1] <= np.log10(pfuthreshold)
            and pflux[i+1] <= np.log10(pfuthreshold)):
            pflux[i]=pflux[i-1]
#             plt.axvline(timestamps[i],color='pink')    
                    
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
                    
                    plt.axvline(a,color='green')
                    plt.axvline(c,color='red')

            except:
                print('Error')
                pass
    ##
    ind = np.where(pflux==np.max(pflux))[0][0]
    print('Spikes fluxes: (i-1,i,i+1)',pflux[ind-1],pflux[ind],pflux[ind+1])
    print('at date',timestamps[ind])  
    ##
    plt.plot(timestamps,pflux)
    hline = np.log10(pfuthreshold)
    plt.hlines(hline,timestamps[0],timestamps[-1], label='>10 MeV',color='steelblue')
    plt.plot([], [], color='green',label="Beginning of SPE")
    plt.plot([], [], color='red', label="Ending of SPE")
    plt.xticks(rotation=50)
    plt.ylabel('Logarithmic proton flux (pfu)')
    plt.ylim(-2,5)
    plt.title(str(year)+'-'+str(month)+' SPE Events >' + str(pfuthreshold) + ' MeV from ' + str(instrument))
    plt.legend(loc=0,frameon=True)
#     plt.savefig(str(year)+'-'+str(month)+str(instrument)+'EPS.png')
    plt.show()
    
plotEPS(1989,7,10,'GOES-06')
# plotEPS(1987,3,10,'GOES-07')
# event at : plotEPS(2006,12,10,'GOES-11')


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

    for i in range(0,len(pflux)):
            
        if (pflux[i] >= np.log10(pfuthreshold) # beginning of event, taking out spike possibilites 
            and pflux[i+1] >= np.log10(pfuthreshold)
            and pflux[i-1] <= np.log10(pfuthreshold)
            and pflux[i+2] >= np.log10(pfuthreshold)
            and pflux[i+3] >= np.log10(pfuthreshold)
            and timestamps[i+1] - timestamps[i] < timedelta(minutes=30) # events with time difference at least 30 minutes from next time instance 
            and event == False):
            print('START: i-1,i,i+1,i+2\n',pflux[i-1],pflux[i],pflux[i+1])
            print('at',timestamps[i])
            start.append(timestamps[i])
            sflux.append(pflux[i])
            event = True

        if (pflux[i] >= np.log10(pfuthreshold) # end of event, saves after event over 
            and pflux[i+1] <= np.log10(pfuthreshold)
            and event==True):
            end.append(timestamps[i+1])
            eflux.append(pflux[i+1])
            event = False 
            
        if (event == False # spike condition, replace flux with previous flux
            and pflux[i]>=np.log10(pfuthreshold)
            and pflux[i-1] <= np.log10(pfuthreshold)
            and pflux[i+1] <= np.log10(pfuthreshold)):
            pflux[i]=pflux[i-1]
#             plt.axvline(timestamps[i],color='pink')    
                    
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
                    plt.axvline(a,color='green')
                    plt.axvline(c,color='red')
            except:
                print('Error')
                pass
    ##
    ind = np.where(pflux==np.max(pflux))[0][0]
    print('Spikes fluxes: (i-1,i,i+1)',pflux[ind-1],pflux[ind],pflux[ind+1])
    print('at date',timestamps[ind])  
    ##
    plt.plot(timestamps,pflux)
    hline = np.log10(pfuthreshold)
    plt.hlines(hline,timestamps[0],timestamps[-1], label='>10 MeV',color='steelblue')
    plt.plot([], [], color='green',label="Beginning of SPE")
    plt.plot([], [], color='red', label="Ending of SPE")
    plt.xticks(rotation=50)
    plt.ylabel('Logarithmic proton flux (pfu)')
    plt.ylim(-2,5)
    plt.title(str(year)+'-'+str(month)+' SPE Events >' + str(pfuthreshold) + ' MeV from ' + instrument)
    plt.legend(loc=0,frameon=True)
#     plt.savefig(str(year)+'-'+str(month)+str(instrument)+'EPEAD.png')
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


# In[ ]:





# In[ ]:




