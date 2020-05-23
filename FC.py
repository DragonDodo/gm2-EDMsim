#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FC implementation

from numpy import linspace
from scipy.stats import norm

def FC_band(x_bins,mu_bins,width,sig):
    #Function to calculate the x-values of the feldman-cousins interval
    
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
            R = norm(val,1).pdf(xs)/norm(mu_max,1).pdf(xs)            
            Rs.append([xs,R])      
    
        #sort by R value, largest first    
        Rarr = sorted(Rs, key=lambda x: x[1],reverse=True)
        
        psum = 0 
        xlimvals = []

        #numerically integrate out adding intervals of highest R first
        for p in Rarr:                       
            
            pv = norm.pdf(p[0],val,1)        
            prob = pv*xwidth #have to scale the pdf prob as it's not normalised right
            psum += prob
            xlimvals.append(p[0])

               
            if psum > sig:
                break
        
        up = max(xlimvals)   
        down = min(xlimvals)        
        
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
        
    return downlim,uplim


def FC_limit(x,sigma,x_bins=linspace(-10,10,400),mu_bins=linspace(0,8,100),alpha=0.9):
    #Helper function to collect it all together and scale the output 
    testx = x/sigma
    
    if testx > max(mu_bins):
        raise ValueError('x/sigma not in band, scan for wider mean value!')
        
    rough = FC_band(x_bins,mu_bins,sigma,alpha)
    lims = find_FC_limits(rough,x_bins,mu_bins,testx) 
    
    print('Limits per sigma:',lims[0],lims[1])
    print('Scaled limits:',lims[0]*sigma,lims[1]*sigma)
    
    return lims[0]*sigma,lims[1]*sigma


if __name__ == "__main__":  

    import matplotlib.pyplot as plt


    plt.style.use('gm2.mplstyle')
    plt.ion()
      
    x = 0
    sigma = 2
    
    l = FC_limit(x,sigma)
    print(l)



    '''
    fct = [-2,-1,0,0.5,1,1.3,1.5,1.7,1.9,2,3] #test x values
    FCpaper = [0.4,0.81,1.64, 2.14, 2.64, 2.94, 3.14,3.34,3.54,3.64, 4.64] #upper limits
    downs = [0,0,0,0,0,0.02,0.22,0.38,0.51,0.58,1.37] #lower limits
    
      
    plt.figure(1)
    plt.plot(lims[2],mu_bins,color='C0',label = 'This implementation')
    plt.plot(lims[3],mu_bins,color='C0')
    #plt.axvline(testx,color='#d3d3d3',linestyle='--')
    #plt.axhline(lims[0],color='#d3d3d3',linestyle='--')
    #plt.axhline(lims[1],color='#d3d3d3',linestyle='--')
    
    
    plt.xlim(0,7)
    plt.ylim(0,5)
    plt.legend()
    plt.xlabel(r"$x/\sigma$",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
    plt.ylabel("$\mu/\sigma$",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0,labelpad=20)
    
    #the values it should be getting are in FCpaper and downs for the up/down limits
    
    fct = [-2,-1,0,0.5,1,1.3,1.5,1.7,1.9,2,3] #test x values
    FCpaper = [0.4,0.81,1.64, 2.14, 2.64, 2.94, 3.14,3.34,3.54,3.64, 4.64] #upper limits
    downs = [0,0,0,0,0,0.02,0.22,0.38,0.51,0.58,1.37] #lower limits
    
    testmus = np.linspace(-2,3,10)
    
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
    
    '''
