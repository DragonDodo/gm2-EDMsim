#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#who needs ROOT anyway?

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import stats

#data here
#GM2 files have no errors included, EDM files do!
t,counts,errors = np.loadtxt("EDM.txt",unpack=True)

# Fit the wiggle with a shifted sine 
fitfunc = lambda p, x: p[0]*np.sin(p[1]*x+p[2]) + p[3] # Target function
errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
p0 = [150, 1.4e-3, 9,2000] # Initial guess for the parameters
p1, success = optimize.leastsq(errfunc, p0[:], args=(t, counts))

#chi-square calculations
ndof = len(t) - len(p1)

fit = fitfunc(p1,t)

residuals = []
for i,j,k in zip(counts,fit,errors):
    res = (i-j)/k
    residuals.append(res**2)
  
    
chisq = sum(residuals)/ndof
print(chisq)

time = np.linspace(t.min(), t.max(), 1000)
plt.errorbar(t, counts, errors) # Plot of the data and the fit
plt.plot(time,fitfunc(p1, time))
plt.xlabel('Time [ns]')
plt.ylabel('Average angle')
plt.text(7500,0.0016,'Amplitude from fit: %f' %p1[0])
plt.text(7500,0.00145,'Chi-squared/ndf: %f' %chisq)
#plt.legend(('Wiggle', 'Fit'))


plt.show()