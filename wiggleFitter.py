
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#who needs ROOT anyway?

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


gamma = 29.3 #'magic' gamma value for g-2

#make plots in g-2 style
plt.style.use('gm2.mplstyle')
plt.rcParams['axes.unicode_minus'] = False #hacky fix for fonts not having minus signs
plt.ion()


def chiSquare(function,datax,datay,fit_params,error_array):
    #function to calculate a general chi-squared value from a fit
    ndof = len(datax) - len(fit_params) 
    fit = function(fit_params,datax) 
    
    residuals = []
    for i,j,k in zip(datay,fit,error_array):
        res = (i-j)/k
        residuals.append(res)

        
    chisq = sum([r**2 for r in residuals])/(ndof)
    return chisq, residuals

def fit_leastsq(p0, datax, datay, function):
    #function to fit using least squares but also calculate parameter errors from the jacobian 
    errfunc = lambda p, x, y: function(p,x) - y #distance to target function

    pfit, pcov, infodict, errmsg, success = optimize.leastsq(errfunc, p0, args=(datax, datay), full_output=1, epsfcn=0.0001) #perform the fit

    if (len(datay) > len(p0)) and pcov is not None:
        s_sq = (errfunc(pfit, datax, datay)**2).sum()/(len(datay)-len(p0))
        pcov = pcov * s_sq #jacobian is the fractional cov matrix so x by residual variance
    else:
        pcov = np.inf

    error = [] 
    for i in range(len(pfit)):
        try:
          error.append(np.absolute(pcov[i][i])**0.5)
        except:
          error.append( 0.00 ) #to stop /0 issues
    pfit_leastsq = pfit
    perr_leastsq = np.array(error) 
    return pfit_leastsq, perr_leastsq 


#data here
tgm2, gm2counts, gm2errors = list(np.loadtxt("GM2.txt",unpack=True,delimiter=','))
t,counts,errors = list(np.loadtxt("EDM.txt",unpack=True,delimiter=','))


# Five-parameter g-2 fit
fitfuncgm2 = lambda p, x: p[0]*np.exp(-x/(gamma*p[1]))*(1-p[2]*np.sin(p[3]*x+p[4])) # 5-parameter wiggle function
p0gm2 = [100,1e-6*1e9,100,1.5e-3,1]
p1gm2, errgm2 = fit_leastsq(p0gm2, tgm2, gm2counts, fitfuncgm2) 
chigm2 = chiSquare(fitfuncgm2,tgm2,gm2counts,p1gm2,gm2errors)[0]


# EDM wiggle fit with a shifted sine
fitfuncedm = lambda p, x: p[0]*np.sin(p[1]*x+p[2]) + p[3] 
p0edm = [5.8e-3,1.45e-3, 0,0] 
p1edm, erredm = fit_leastsq(p0edm,t,counts,fitfuncedm)
chiedm = chiSquare(fitfuncedm,t,counts,p1edm,errors)[0]


plt.figure(1)
time = np.linspace(tgm2.min(), tgm2.max(), 1000)
plt.fill_between(tgm2,gm2counts-gm2errors,gm2counts+gm2errors,alpha=0.7,color='#ff7f0e')
plt.scatter(tgm2,gm2counts,marker='.', color='k')
plt.plot(time,fitfuncgm2(p1gm2, time),color='r')
plt.xlabel("Time [ns]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
plt.ylabel("Count with E>2000 MeV/300ns",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0)
plt.text(30000,3800,r'$\chi^{2}/dof: %.2f$' %chigm2,horizontalalignment='right')
plt.text(30000,3500,r'$N_{0}: (%.2f \pm %.2f) \times10^{3}$' %(p1gm2[0]*1e-3,errgm2[0]*1e-3),horizontalalignment='right')
plt.text(30000,3200,r'$\tau: (%.2f \pm %.2f) \times10^{2} ns$' %(p1gm2[1]*1e-2,errgm2[1]*1e-2),horizontalalignment='right')
plt.text(30000,2900,r'$A: %.2f \pm %.3f$' %(p1gm2[2],errgm2[2]),horizontalalignment='right')
plt.text(30000,2600,r'$\omega: (%.2f \pm %.3f) \times10^{-3} ns^{-1}$' %(p1gm2[3]*1e3,errgm2[3]*1e3),horizontalalignment='right')
plt.text(30000,2300,r'$\phi: %.2f \pm %.2f$' %(p1gm2[4],errgm2[4]),horizontalalignment='right')


plt.figure(2)
time = np.linspace(t.min(), t.max(), 1000)
plt.fill_between(t,counts-errors,counts+errors,alpha=0.8,color='#942ed2')
plt.scatter(t,counts,marker='.', color='k')
plt.plot(time,fitfuncedm(p1edm, time),color='r')
plt.xlim(0,max(t))
plt.xlabel("Time % g-2 period [ns]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
plt.ylabel("Average vertical angle [rad]",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0);
plt.text(4000,0.0034,r'$\chi^{2}/dof: %.2f$' %chiedm,horizontalalignment='right')
plt.text(4000,0.004,r'$Amplitude: (%.2f \pm %.2f) \times10^{-3} rad$' %(p1edm[0]*1e3,erredm[0]*1e3),horizontalalignment='right')


plt.show()



