#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 11:03:02 2020

@author: dodo
"""
#FC implementation

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

#plt.style.use('gm2.mplstyle')
#plt.ion()


#start with a test mu value to find the limit for 
#find the R ratio of that gaussian/the mu best gaussian for each x
#order them based on R and add points to the band until the prob in the band hits desired sig
#get min/max from that for each mu: that defines the band 

#this code only tries to calculate one mu atm 

def FC(x_bins,mean,width,sig,val):
    
    def Rscan(bins):
        xwidth = (max(bins)-min(bins))/len(bins)
        Rs = []          
    
        for xs in bins:
            if xs >= 0:
                mu_max = xs
            elif xs < 0:
                mu_max = 0

                
            '''
            #this is the same calculation done numerically, also doesn't work
            if xs >= 0:
                R = np.exp(-0.5*(xs-val)*(xs-val))
            elif xs < 0:
                R = np.exp(0.5*(xs*val-val*val))
            '''
                       

            R = norm(val,width).pdf(xs)/norm(mu_max,width).pdf(xs)            
            Rs.append([xs,R])      
    
            
        Rarr = sorted(Rs, key=lambda x: x[1],reverse=True)
        #print(Rarr)
        
        psum = 0 
        xlimvals = []



        for p in Rarr:                       
            
            pv = norm.pdf(p[0],val,width)        
            #ab0 = norm.cdf(0,val+mean,width) 
            prob = pv*xwidth #have to scale the pdf prob as it's not normalised right
            psum += prob
            xlimvals.append(p[0])
            #print(psum)
               
            if psum > sig:
                break
        
        print(xlimvals)
        print(psum)
        up = max(xlimvals)   
        down = min(xlimvals)
        
        
        if up == max(x_bins):
            print('Warning: upper limit setting not complete, check mu range!')
         
                  
        return down,up+xwidth,xwidth
    
    firstscan = Rscan(x_bins) 
    
   
    print(firstscan)

    return firstscan



mean = 0
width = 1


x_bins=np.linspace(-10,10,1000)
rough = FC(x_bins,mean,width,0.90,3)


#the values it should be getting are in FCpaper and downs for the up/down limits
'''
FCpaper = [0.4,0.81,1.64, 2.14, 2.64, 2.94, 3.14,3.34,3.54,3.64, 4.64]
downs = [0,0,0,0,0,0.02,0.22,0.38,0.51,0.58,1.37]

testmus = np.linspace(-2,3,10)
fct = [-2,-1,0,0.5,1,1.3,1.5,1.7,1.9,2,3]
ls =[]
ds = []

for t in testmus:
    lim = FC(x_bins,mean,width,0.90,t)
    ls.append(lim[1])
    ds.append(lim[0])
    

plt.figure(1)  
plt.plot(testmus,ls,color='C0', label='This implementation')
plt.plot(testmus,ds,color='C0')
plt.plot(fct,FCpaper,color='C1', label='FC paper values')
plt.plot(fct,downs,color='C1')

#plt.ylim(bottom=0)
#plt.xlim(0,4)
plt.xlabel(r"$x/\sigma$",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
plt.ylabel("$\mu/\sigma$",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0,labelpad=20)
'''
