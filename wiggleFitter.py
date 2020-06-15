
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#who needs ROOT anyway?

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.stats import norm
import FC

gamma = 29.3 #'magic' gamma value for g-2

#make plots in g-2 style
plt.style.use('gm2.mplstyle')
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

def EDMamplitude_to_tilt(a):
    consts = 9.158e15       
    d_u = np.arctan(gamma*np.tan(a/consts))
    return d_u

def limitCalc(a,err,CL):
    power = np.floor(np.log10(err))
    amp = a*10**(-1*power)
    err = err*10**(-1*power)
    
    l = FC.FC_limit(amp,err,0,x_bins=x_bins,mu_bins=mu_bins,alpha=CL)
    up = l[1]*10**(power)
    down = l[0]*10**(power)
    
    return down,up
    



#data here
tgm2, gm2counts, gm2errors = list(np.loadtxt("gm2.txt",unpack=True,delimiter=','))
t,counts,errors = list(np.loadtxt("hadd.txt",unpack=True,delimiter=','))


# Five-parameter g-2 fit
fitfuncgm2 = lambda p, x: p[0]*np.exp(-x/(gamma*p[1]))*(1-p[2]*np.sin(p[3]*x+p[4])) # 5-parameter wiggle function
p0gm2 = [100,1e-6*1e9,100,1.5e-3,1]
p1gm2, errgm2 = fit_leastsq(p0gm2, tgm2, gm2counts, fitfuncgm2) 
chigm2 = chiSquare(fitfuncgm2,tgm2,gm2counts,p1gm2,gm2errors)[0]

w = p1gm2[3]
phi = p1gm2[4]-np.pi/2

# EDM wiggle fit with a shifted sine
fitfuncedm = lambda p, x: p[0]*np.sin(w*x+phi) + p[1]
p0edm = [0,0] 
p1edm, erredm = fit_leastsq(p0edm,t,counts,fitfuncedm)
chiedm = chiSquare(fitfuncedm,t,counts,p1edm,errors)[0]

#limit setting stuff 
fitTilt = EDMamplitude_to_tilt(p1edm[0])
fitTilterr = EDMamplitude_to_tilt(erredm[0])

xrange = np.linspace(-1e-18,1e-18,1000)
limitGaus = norm.pdf(xrange,fitTilt,fitTilterr)

interval = norm.interval(0.90,loc=fitTilt,scale=fitTilterr)
#limit = interval[1]
#print(limit)

#this can be a little tempramental - if you get errors, adjust the ranges until you don't or limits will be too low
x_bins=np.linspace(0,40,1000)
mu_bins=np.linspace(0,20,1000)




lim = limitCalc(fitTilt,fitTilterr,0.9)
print(lim)

#----------------------






#---------------------





plt.figure(1)
time = np.linspace(tgm2.min(), tgm2.max(), 1000)
plt.fill_between(tgm2,gm2counts-gm2errors,gm2counts+gm2errors,alpha=0.7,color='#ff7f0e')
plt.scatter(tgm2,gm2counts,marker='.', color='k')
plt.plot(time,fitfuncgm2(p1gm2, time),color='r')
plt.xlabel("Time [ns]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
plt.ylabel("Count with E>2000 MeV /300ns",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0,labelpad=20)
plt.text(30000,1200,r'$\chi^{2}/dof: %.2f$' %chigm2,horizontalalignment='right')
plt.text(30000,1100,r'$N_{0}: (%.2f \pm %.2f) \times10^{3}$' %(p1gm2[0]*1e-3,errgm2[0]*1e-3),horizontalalignment='right')
plt.text(30000,1000,r'$\tau: %.2f \pm %.2f \times10^{3} ns$' %(p1gm2[1]*1e-3,errgm2[1]*1e-3),horizontalalignment='right')
plt.text(30000,900,r'$A: %.3f \pm %.3f$' %(p1gm2[2],errgm2[2]),horizontalalignment='right')
plt.text(30000,800,r'$\omega: (%.3f \pm %.3f) \times10^{-3} ns^{-1}$' %(p1gm2[3]*1e3,errgm2[3]*1e3),horizontalalignment='right')
plt.text(30000,700,r'$\phi: %.2f \pm %.2f$' %(p1gm2[4],errgm2[4]),horizontalalignment='right')


plt.figure(2)
time = np.linspace(t.min(), t.max(), 1000)
plt.fill_between(t,counts-errors,counts+errors,alpha=0.8,color='#942ed2')
plt.scatter(t,counts,marker='.', color='k')
plt.plot(time,fitfuncedm(p1edm, time),color='r')
plt.xlim(0,max(t))
#plt.ylim(top=0.8e-3)
plt.gca().ticklabel_format(style='sci', scilimits=(0,1), axis='y')
plt.xlabel("Time % g-2 period [ns]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
plt.ylabel("Average vertical angle [rad]",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0,labelpad=20);
plt.text(4000,0.7e-3,r'$\chi^{2}/dof: %.2f$' %(chiedm),horizontalalignment='right')
plt.text(4000,0.6e-3,r'$Amplitude: (%.2f \pm %.2f) \times10^{-5} rad$' %(p1edm[0]*1e5,erredm[0]*1e5),horizontalalignment='right')

xlow = np.linspace(0,1.68*2*fitTilterr,100)
band = norm.pdf(xlow,fitTilt,fitTilterr)




plt.figure(3)

plt.plot(xrange,limitGaus,label = 'Probability gaussian')
'''

plt.xlabel("EDM value [e cm]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0,labelpad = 40)
plt.ylabel("Counts/bin [arbitrary units]",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=1.0)
plt.axvline(lim[0],color='#d3d3d3',linestyle='--', label='90% FC limit')
plt.axvline(lim[1],color='#d3d3d3',linestyle='--')
plt.axvline(fitTilt,color='C0',linestyle='--')
#plt.text(max(xrange),max(limitGaus)*0.95,r'$|d_{u}|\ < %.2f \times 10^{-19} e\ cm$' %(interval[1]*1e19),horizontalalignment='right')
plt.legend(loc='upper left')
'''






