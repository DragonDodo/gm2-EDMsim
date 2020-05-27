#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
#script to implement the CLs method for setting limits
'''

from scipy.stats import norm

def findCLsLimit(xrange,s,b,sig,err):
    #function to calc the CLs defined limit
    limit = 0
    for t in xrange:
        pb = norm.cdf(s,b,err)
        psb = norm.cdf(s,t+b,err)
        CLs = psb/pb
        #print(t, CLs, pb, psb)
        
        if 1-CLs <= sig: #iterate over bins to find limit CLs = sig
            limit = t
            
    
    con_limit = limit + (xrange[1]-xrange[0]) #add extra bin to overcover
    return con_limit
 
if __name__ == "__main__":  
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    #make plots in g-2 style
    plt.style.use('gm2.mplstyle')
    plt.ion() 
    
    

           
    xrange = np.linspace(0,10,100) #this range is what you scan over: so setting it from 0 limits the edm > 0. 
    amp = 0
    err = 1
    
    signal = 0.4
    background = 0
    
    sol = findCLsLimit(xrange,signal,background,0.9,err)
    print(sol)
    
    
    '''
    amps = np.linspace(-2,8,40) #input du values (not really amplitudes any more)
    
    clist = []
    olist = []
    
    for i in amps:
        sol = findCLsLimit(xrange,i,background,0.9,err)  
        
        clist.append(sol)
        
        
        
    plt.figure(1)  
    plt.plot(amps,clist,label='CLs method')  
    plt.xlim(min(amps),max(amps))
    plt.ylim(0,max(clist))
    plt.xlabel(r"$d_{u}\ from\ amplitude\ [e\ cm]$",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0, labelpad=40)
    plt.ylabel("Limit set [e cm]",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0,labelpad=20)
    '''
        
    
