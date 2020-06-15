#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FC implementation

from numpy import linspace,exp,sqrt,pi
from scipy.stats import norm
import CLs


def folded_norm(x,bin_width,mu,sigma):
    upbranch = exp(-1*(x-mu)*(x-mu)/(2*sigma*sigma))
    downbranch = exp(-1*(x+mu)*(x+mu)/(2*sigma*sigma))
    
    return (upbranch+downbranch)*1/(sigma*sqrt(2*pi))*bin_width


def FC_band(x_bins,mu_bins,s_width,significance):
    #Function to calculate the x-values of the feldman-cousins interval
    #Normalised by sigma for speed! All scaling is done by helper functions.
    
    #x_bins - range of x to scan over
    #width - standard deviation of the gaussian
    #sig - desired significance level 
    #val - measured x to set the limit using
    print('Calculating band...')
    
    def Rscan(bins,val):
        #function to scan for a single mu value for all x
        xwidth = (max(bins)-min(bins))/len(bins)
        Rs = []          
    
    
        #define mu_max for all x - mu_max = max(0,x)
        for xs in bins:
            if xs >= 0:
                mu_max = xs
            elif xs < 0:
                mu_max = 0
    
            #R ratio P(x|mu)/P(x|mu_max)
            R = folded_norm(xs,xwidth,val,1)/folded_norm(xs,xwidth,mu_max,1)
            #R = norm(val,1).pdf(xs)/norm(mu_max,1).pdf(xs)            
            Rs.append([xs,R])      
    
        #sort by R value, largest first    
        Rarr = sorted(Rs, key=lambda x: x[1],reverse=True)
        
        psum = 0 
        xlimvals = []

        #numerically integrate out adding intervals of highest R first
        for p in Rarr:                       
            
            #pv = norm.pdf(p[0],val,1)   
            pv = folded_norm(p[0],xwidth,val,s_width)
            
            #prob = pv*xwidth #have to scale the pdf prob as it's not normalised right
            prob = pv
            psum += prob
            xlimvals.append(p[0])

               
            if psum > significance:
                break
        
        up = max(xlimvals)
        down = min(xlimvals)
        #print(psum)
        
        if up == max(x_bins):
            print('Warning: upper limit setting not complete, check mu range!')         
        #print(down,up)         
        return down,up
    
    band = []
    for mu in mu_bins:    
        scan = Rscan(x_bins,mu) 
        band.append(scan)

    return band


def find_FC_limits(band,x_bins,mu_bins,testval): 
    #Function to convert the band into limits on x
    print('Setting limits...')
    
    def limscan(xs,mu_bins,testx):
        for i in range(len(xs)):
            lim = mu_bins[i]
            if xs[i] >= testx:
                break
            
        return lim      
    
    upx = []
    downx = []
    
    for i,j in zip(band,mu_bins):
        downx.append(i[1])
        upx.append(i[0])
        
    testx = testval
        
    
    downlim = limscan(downx,mu_bins,testx)
    uplim = limscan(upx,mu_bins,testx) 

        
    if downlim < 0:
        downlim = 0.0
        
    return downlim,uplim,upx,downx


def FC_limit(x,s_width,bgmean,x_bins=linspace(-10,10,400),mu_bins=linspace(0,8,100),alpha=0.9):
    #Helper function to collect it all together and scale the output 
    
    testx = abs(x-bgmean)
    
    binwidth = (max(x_bins)-min(x_bins))/len(x_bins)
    
    if testx ==0:
        testx = binwidth #smallest limit that can be calcuated with the step size
        
    print(testx)

    
    
    if testx > max(mu_bins):
        raise ValueError('x/sigma not in band, scan for wider mean value!')
        
    rough = FC_band(x_bins,mu_bins,s_width,alpha)
    lims = find_FC_limits(rough,x_bins,mu_bins,testx) 
    
    print('Limits (x10^(n)):',lims[0]*s_width,lims[1]*s_width)
    
    return lims[0]*s_width,lims[1]*s_width


def one_sided_limit(mean,sigma,x_bins,alpha=0.9):
    
    xwidth = (max(x_bins)-min(x_bins))/len(x_bins)
    
    psum = 0
    uplim = 0
    for i in x_bins:
        pv = folded_norm(i,xwidth,mean,sigma)
        #pv = norm.pdf(i,mean,sigma)*xwidth
        psum += pv
        
        if psum >= alpha:
            uplim = i
            break
        
    return uplim
            

if __name__ == "__main__":  

    import matplotlib.pyplot as plt
    import numpy as np


    plt.style.use('gm2.mplstyle')
    plt.ion()
      
    
    x_bins=linspace(0,20,1000)
    mu_bins=linspace(0,10,1000)
    xwidth = (max(x_bins)-min(x_bins))/len(x_bins)
    
    x = 10
    sigma = 4
    #sigma = 1
    
    lim2 = FC_limit(x,sigma,0,x_bins,mu_bins,0.9)

    '''
    plt.figure(1)
    plt.plot(x_bins,folded_norm(x_bins,xwidth,4,2),label=r'$\mu=4, \sigma=2$')
    plt.xlabel("x",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
    plt.ylabel("PDF of x",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0,labelpad=20)
    plt.legend(loc='upper right')
    '''

    
    #halfgaus = folded_norm(x_bins,xwidth,2,2)
    #normed = halfgaus*xwidth


    '''

    
    

    testmus = np.linspace(0,8,100)
    
    ls =[]
    ds = []
    ls2 =[]
    ds2 = []
    
    gup=[]
    gdown=[]
    
    sid = []    
    cid = []


    x = 1
    sigma = 1
    

    
       
    
        
    
    
    lim = FC_band(x_bins,testmus,sigma,0.9)
    lim2 = FC_limit(x,sigma,0,x_bins,mu_bins,0.9)
    onesid = one_sided_limit(x,sigma,x_bins,0.95)
    interval = norm.interval(0.95,loc=x,scale=sigma)
    print(lim2)
    
    for i in lim:
        ls.append(i[0])
        ds.append(i[1])
        
        
  
    plt.figure(1)  
    plt.plot(testmus,ls2)
    plt.plot(testmus,ds2,label=r'$\mu_{bkg}=0.0$')
    
    

    
    plt.xlabel(r"$Measured\ value\ x$",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
    plt.ylabel("$Parameter\ value\ \mu$",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0,labelpad=20)
    plt.legend(loc='upper left')
    
       
    for i in testmus:
        interval = norm.interval(0.9,loc=i,scale=1) 
        gup.append(interval[1])
        gdown.append(interval[0])
        onesid = one_sided_limit(i,1,x_bins)
        sid.append(onesid)
        #sol = CLs.findCLsLimit(x_bins,i,0,0.9,sigma) 
        #cid.append(sol)





    fct = [-2,-1,0,0.5,1,1.3,1.5,1.7,1.9,2,3] #test x values
    FCpaper = [0.4,0.81,1.64, 2.14, 2.64, 2.94, 3.14,3.34,3.54,3.64, 4.64] #upper limits
    downs = [0,0,0,0,0,0.02,0.22,0.38,0.51,0.58,1.37] #lower limits
    
      
    plt.figure(1)
    #plt.plot(fct,FCpaper,color='C1', label='FC paper values')
    #plt.plot(fct,downs,color='C1')
    plt.plot(ls,testmus,color='r',label = 'Feldman-Cousins method')
    plt.plot(ds,testmus,color='r')
    plt.plot(testmus,gup,color='black',linestyle='--', label='Two-sided gaussian')
    plt.plot(testmus,gdown,color='black',linestyle='--')
    plt.plot(testmus,sid,color='blue', label = 'One-sided folded gaussian')
    plt.plot(testmus,[0 for i in testmus],color='blue')
    #plt.plot(testmus,cid,color='green', label = 'CLs method')

    
    plt.xlim(0,2.5)
    plt.ylim(-3,4)
    plt.legend(loc='lower right')
    plt.xlabel(r"$x/\sigma$",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
    plt.ylabel("$\mu/\sigma$",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0,labelpad=20)

    
    plt.axvline(0.4,color='#d3d3d3',linestyle='--')
    

    
    plt.axhline(l[0]/sigma,color='#d3d3d3',linestyle='--')
    plt.axhline(l[1]/sigma,color='#d3d3d3',linestyle='--')
    
    
    

    
    
    #the values it should be getting are in FCpaper and downs for the up/down limits
    
    fct = [-2,-1,0,0.5,1,1.3,1.5,1.7,1.9,2,3] #test x values
    FCpaper = [0.4,0.81,1.64, 2.14, 2.64, 2.94, 3.14,3.34,3.54,3.64, 4.64] #upper limits
    downs = [0,0,0,0,0,0.02,0.22,0.38,0.51,0.58,1.37] #lower limits
    

    
    plt.figure(1)  
    plt.plot(testmus,ls,color='C0', label='This implementation')
    plt.plot(testmus,ds,color='C0')
    
    #plt.ylim(bottom=0)
    #plt.xlim(0,4)
    
    '''
