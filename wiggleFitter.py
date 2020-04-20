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
tgm2, gm2counts, gm2errors = list(np.loadtxt("GM2.txt",unpack=True,delimiter=','))
t,counts,errors = list(np.loadtxt("EDM.txt",unpack=True,delimiter=','))



T= 2*np.pi/1.5e-3
modt = np.array([i%T for i in t])

'''
nbins = 30
n, _ = np.histogram(t, bins=nbins)
sy, _ = np.histogram(t, bins=nbins, weights=counts)
sy2, _ = np.histogram(t, bins=nbins, weights=counts*counts)
binned_counts = sy / n
std = np.sqrt((sy2/n) - binned_counts**2)
bin_edges = (_[1:] + _[:-1])/2
'''




# Fit the wiggle with a shifted sine 
fitfunc = lambda p, x: p[0]*np.sin(p[1]*x+p[2]) + p[3] # Target function
errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
p0gm2 = [1000,1.5e-3,1,2000]

p1g2, successgm2 = optimize.leastsq(errfunc, p0gm2[:], args=(tgm2, gm2counts))

#T = 2*np.pi/p1g2[1]
#print(T)



p0edm = [5.8e-3,1.45e-3, 0,0] # Initial guess for the parameters
p1edm, success = optimize.leastsq(errfunc, p0edm[:], args=(t, counts))


def chiSquare(time_array,data_array,fit_params,error_array):
    ndof = len(time_array) - len(fit_params) 
    fit = fitfunc(fit_params,time_array)
    
    residuals = []
    for i,j,k in zip(data_array,fit,error_array):
        res = (i-j)/(k)
        residuals.append(res)      
        
    chisq = sum([r**2 for r in residuals])/(ndof)
    return chisq

chilist = [0.9064453676110634,1.0632147174507445,0.9605416617436612,0.9134226996919688,1.2105460796003507,1.1871313386931648,1.0212660986090978,1.1772792026978454,0.8902910376299037, 1.1634836078839093,0.8002650149994733, 0.9054801071526796, 1.0720157244884008, 0.9207173192895274
,0.9208960760747018,1.043713058619897,1.023182588121364,1.0937785104251003,0.9321739898677243,1.0180883591527838,1.027715890058236
,0.904279760630539,0.9811697305892801,0.9583329576274016,1.0656585083170707,1.1088759752509622,0.9752644626152667,0.9362351888583991
,1.1024195116300794,0.926048765829955]

#print(np.mean(chilist))


fig, ax = plt.subplots()

#plt.figure(1)
#plt.hist(chilist,bins=6,histtype='step')
'''
atts = np.linspace(1,30,30)
plt.scatter(atts,chilist,marker='o',color='k',label='Fit results')
mean = np.mean(chilist)
plt.axhline(mean, linewidth = 2, color = 'b',label='Mean value')

ax.set_xlabel("Chi-squared value",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
ax.set_ylabel("Count/bin",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0)
ax.legend()
'''

plt.figure(1)
time = np.linspace(tgm2.min(), tgm2.max(), 1000)
plt.fill_between(tgm2,gm2counts-gm2errors,gm2counts+gm2errors,alpha=0.7,color='#ff7f0e')
plt.scatter(tgm2,gm2counts,marker='.', color='k')
plt.plot(time,fitfunc(p1g2, time),color='r',linewidth=2)
plt.xlabel("Time [ns]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
plt.ylabel("Number of positrons with E>2000 MeV",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0);


chigm2 = chiSquare(tgm2,gm2counts,p1g2,gm2errors)
chiedm = chiSquare(t,counts,p1edm,errors)


print(chiedm)


plt.figure(2)
time = np.linspace(t.min(), t.max(), 1000)
#plt.errorbar(t, counts, errors,ls='none') # Plot of the data and the fit
plt.fill_between(t,counts-errors,counts+errors,alpha=0.8,color='#942ed2')
plt.scatter(t,counts,marker='.', color='k')

#plt.errorbar(bin_edges, mean, yerr=std, fmt='bo')
plt.plot(time,fitfunc(p1edm, time),linewidth=2,color='r')
#plt.fill_between(time,fitfunc(pd, time),fitfunc(pu, time),color='orange',alpha=0.7)


#plt.plot(time,fitfunc(p0,time))


plt.xlim(0,max(t))
plt.xlabel("Time [ns]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
plt.ylabel("Average vertical angle [rad]",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0);
#plt.text(2100,0.0070,'Amplitude from fit: %f' %p1[0],horizontalalignment='left')
#plt.text(2100,0.0060,'Chi-squared/ndf: %f' %chisq,horizontalalignment='left')
#plt.legend(('Wiggle', 'Fit'))



input('Any key to continue')

plt.close()



plt.show()



