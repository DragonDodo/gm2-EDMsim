#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#who needs ROOT anyway?

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import mplhep as hep

#make plots in g-2 style
plt.style.use(hep.style.ATLAS)
plt.style.use({"axes.labelsize":'20'})
plt.rcParams['axes.unicode_minus'] = False #hacky fix for fonts not having minus signs
plt.ion()

#data here
#GM2 files have no errors included, EDM files do!
tgm2, gm2counts = list(np.loadtxt("GM2.txt",unpack=True,delimiter=','))
t,counts,errors = list(np.loadtxt("EDM.txt",unpack=True))



T= 2*np.pi/1.5e-3
modt = np.array([i%T for i in t])


nbins = 30
n, _ = np.histogram(t, bins=nbins)
sy, _ = np.histogram(t, bins=nbins, weights=counts)
sy2, _ = np.histogram(t, bins=nbins, weights=counts*counts)
binned_counts = sy / n
std = np.sqrt((sy2/n) - binned_counts**2)
bin_edges = (_[1:] + _[:-1])/2





# Fit the wiggle with a shifted sine 
fitfunc = lambda p, x: p[0]*np.sin(p[1]*x+p[2]) + p[3] # Target function
errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
p0 = [1000,1.5e-3,1,2000]

p1g2, successgm2 = optimize.leastsq(errfunc, p0[:], args=(tgm2, gm2counts))

#T = 2*np.pi/p1g2[1]
#print(T)



p0edm = [5.8e-3,1.45e-3, 0,0] # Initial guess for the parameters
p1, success = optimize.leastsq(errfunc, p0edm[:], args=(t, counts))

    
fit = fitfunc(p1,t)


ndof = len(t) - len(p0edm)   

    
residuals = []
for i,j,k in zip(counts,fit,errors):
    res = (i-j)/(k)
    residuals.append(res)
        
    
      
        
chisq = sum([r**2 for r in residuals])/(ndof)
print(chisq)




fig, ax = plt.subplots()
'''
plt.figure(1)
time = np.linspace(tgm2.min(), tgm2.max(), 1000)
plt.scatter(tgm2,gm2counts,marker='.')
plt.plot(time,fitfunc(p1g2, time))
ax.set_xlabel("Time [ns]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
ax.set_ylabel("Number of positrons with E>2000 MeV",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0);
'''
plt.figure(1)
time = np.linspace(t.min(), t.max(), 1000)
plt.errorbar(t, counts, errors,ls='none') # Plot of the data and the fit

#plt.errorbar(bin_edges, mean, yerr=std, fmt='bo')
plt.plot(time,fitfunc(p1, time),linewidth=2)
#plt.fill_between(time,fitfunc(pd, time),fitfunc(pu, time),color='orange',alpha=0.7)


#plt.plot(time,fitfunc(p0,time))


plt.xlim(0,max(t))
ax.set_xlabel("Time % g-2 frequency [ns]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
ax.set_ylabel("Average angle [rad]",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0);
#plt.text(2100,0.0070,'Amplitude from fit: %f' %p1[0],horizontalalignment='left')
#plt.text(2100,0.0060,'Chi-squared/ndf: %f' %chisq,horizontalalignment='left')
#plt.legend(('Wiggle', 'Fit'))



input('Any key to continue')

plt.close()



plt.show()



