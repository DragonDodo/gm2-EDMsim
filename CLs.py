#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
#script to implement the CLs method for setting limits
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

#make plots in g-2 style
plt.style.use('gm2.mplstyle')
plt.ion()

def pVal(mean,sigma,x):
    #functional shorthand to calc p-value
    p = norm.cdf(x,mean,sigma)
    return p

def findCLsLimit(xrange,amp,err,sig):
    #function to calc the CLs defined limit
    limit = 0
    for t in xrange:
        pb = pVal(0,err,t)
        psb = pVal(amp,err,t)
        CLs = psb/pb
        #print(t, CLs)
        
        if CLs <= sig: #iterate over bins to find limit CLs = sig
            limit = t
            
    
    con_limit = limit + (xrange[1]-xrange[0]) #add extra bin to overcover
    return con_limit
            
xrange = np.linspace(0,100e-19,1000) #this range is what you scan over: so setting it from 0 limits the edm > 0. 
amp = 0.3e-19
err = 2e-19

amps = np.linspace(0e-19,10e-19,40) #input du values (not really amplitudes any more)

clist = []
olist = []

for i in amps:
    sol = findCLsLimit(xrange,i,err,0.95) # CLs   
    pb = norm.interval(0.95,loc=i,scale=err)[1] #old but two-tailed method
    clist.append(sol)
    olist.append(pb)
    
    
plt.figure(1)  
plt.scatter(amps,clist,label='CLs method')  
plt.scatter(amps,olist,label='Naive frequentist')
plt.xlim(min(amps),max(amps))
plt.ylim(min(clist),max(olist))
plt.xlabel(r"$d_{u}\ from\ amplitude\ [e\ cm]$",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0, labelpad=40)
plt.ylabel("Limit set [e cm]",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0,labelpad=20)
    
    
    
