#!/usr/bin/env
# -*- coding: utf-8 -*-

'''
TODO:
    -Test EDM oscialltion
    -Test fitting and chi-squared calcs
    - Check impact on w_tot by EDM
'''

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spint
from scipy.stats import norm
from scipy.optimize import curve_fit



#Crude 'mode' selection: 
#Options:
#   'EDM' to make vertical angle oscialltion plots
#   'GM2' to make usual wiggle plot
#   Blank string for single angle stationary snapshot


option = "EDM"
npoints = 100 #number of points to plot for wiggles


#MeV units
n_muons = 10000
muonM = 105.66
eM = 0.51 
Emax = muonM/2

class muon:
    #class to define a muon decay object
    def __init__(self):
        eE = sample_muE(1)
        self.E = eE
        self.angle=sample_angle(1,eE)
        self.lifetime = 2.1969811e-6 #seconds
        self.P = np.sqrt(self.E**2-eM**2)
        self.vangle = sample_vertical_angle(1)
        
        #unchanged self variables for debugging
        self.oE = self.E
        self.oA = self.angle
        self.oP = self.P
        self.oV = self.vangle
    
    def boost(self):
        
        #method to boost the muon into the detector frame  
        muonM = 105.66
        gamma = 29.3 #'magic' gamma value for g-2
        mu_momentum = 3.094e3 #3.094 GeV/c
        mu_v = -mu_momentum/(gamma*muonM) #-tive because boosting back into lab frame
        
                
        
        Px = self.P*np.sin(self.angle)*np.cos(self.vangle)
        Py = self.P*np.sin(self.vangle)
        Pz = self.P*np.cos(self.angle)*np.cos(self.vangle)
        
        boostedPz = gamma*(Pz-mu_v*self.E)
        boostedE = gamma*(self.E-mu_v*Pz) 
        
        self.E = boostedE
        self.angle = np.arctan(Px/boostedPz)
        self.vangle = np.arctan(Py/boostedPz)
        self.lifetime = self.lifetime*gamma
        self.P = np.sqrt(boostedPz**2+Px**2)
        
    def gm2Shift(self,t):
    #method to add in the g-2 oscillation by rotating spin vector by an angle
    #IMPORTANT: must do this before boosting!
        w = 1.5e-3
        shift = w*t
        self.angle=self.angle+shift
        self.oA = self.angle 
        
    def EDMShift(self,t):
    #method to add in an EDM oscillation, ie a vertical angle oscillation
    #Do before boosting!
        
        wEDM = 1.5e-3 # ~ g-2 freq   
        #add in a direct calculation here!
        
        d_u = 10*1.9e-19 #100BNL
        consts = 8.9e15
        
        A = d_u*consts
        phi = np.pi/2# EDM is pi/2 out of phase with g-2 oscillation        
        self.vangle = self.vangle + A*np.cos(wEDM*t-phi)   
        #self.vangle = self.vangle + t

#-------------------------------------------------------
#Define the energy distribution of muon decay and sample

def energyDecaySpectrum(E):
    #muon decay energy spectrum, assuming massless positrons     
    
    prefactor = 1e-5 #scale it close to 1
    bracket = (E**2)*np.float128((1-(4*E)/float(3*muonM)))
   
    return prefactor*bracket


def sample_muE(n): #using rejection sampling to get a sample of size n
    samples = []

    while len(samples) < n:
        z = np.random.uniform(eM, Emax)
        u = np.random.uniform(0, 0.05)

        if u <= energyDecaySpectrum(z):
            samples.append(z)

    if len(samples) ==1:
        return samples[0]
    else:
        return samples

#-------------------------------------------------------
#Define angular distribution for positrons and sample
        
def angleSpectrum(E,theta,Pu=1):
    
    x = E/Emax

    a = (Pu*np.cos(theta))*(1-2*x)
    b = 3-2*x    
    # 1 + 1/3cos(theta)
    return x*x*(b-a)

def sample_angle(n,E):
    
    samples = []

    while len(samples) < n:
        z = np.random.uniform(-np.pi, np.pi)
        u = np.random.uniform(0, 5.)

        if u <= angleSpectrum(E,z,1):
            samples.append(z)

    if len(samples) ==1:
        return samples[0]
    else:
        return samples
   

     
#--------------------------------------------------------    
#vertical distribution of positrons
'''
def vertical_angle(n):
    #vA = np.random.normal(0,0.01)
    #vA = np.random.uniform(-1*np.pi/2,np.pi/2)
    cos =  cosine.rvs(1)
    #vA = np.pi*2*cos -1
    return cos
'''

def verticalAngleSpectrum(theta):    
    

    a = np.cos(theta)
      
    # 1 + 1/3cos(theta)
    return a

def sample_vertical_angle(n):
    
    samples = []

    while len(samples) < n:
        z = np.random.uniform(-np.pi/2, np.pi/2)
        u = np.random.uniform(0, 20.)

        if u <= verticalAngleSpectrum(z):
            samples.append(z)

    if len(samples) ==1:
        return samples[0]
    else:
        return samples

#--------------------------------------------------------


def decayMuons(n,t,e):
    muonlist = []
    
    for i in range(n+1):
        m = muon()
        m.gm2Shift(t) #apply g-2 precession
        m.EDMShift(e) #apply EDM precession
        m.boost() #boost into lab frame  
        #m.boostback()
        muonlist.append(m)  
    return muonlist
   
#-------------------------------------------------------

if option == "GM2":   
    
   
    times = list(np.linspace(0., 30e3, npoints))
        
    counts = []
    for t in times:
        counter = 0
        genmuons = decayMuons(n_muons,t,0)
        print('Generating muons for time', t)
        
        for i in genmuons:
            #energy cut for wiggle plot
            if i.E> 2000:
                counter += 1
                
            
        
        counts.append(counter)   
    
    
    plt.figure(1)
    plt.scatter(times,counts,marker='.',label='')    
    plt.xlabel('Time [ns]')
    plt.ylabel('Number of positrons with E>2000 MeV')
    plt.legend()
    
    f=open('precession.txt','w')
    for time,count in zip(times,counts):
        f.write(str(time)+','+str(count)+'\n')
    f.close()
    

    

elif option == "EDM":
    
    times = list(np.linspace(0., 30e3, npoints))
    
    avAngle = []
    spreadAngle = []
    counts = []

    
    for t in times:
        counter = 0
        edmlist = []
        genmuons = decayMuons(n_muons,t,t)
        
        #gm2
        for i in genmuons:
            edmlist.append(i.vangle)
            
            if i.E> 2000:
                counter += 1
        
        counts.append(counter) 
        #plt.figure(3)
        
        b = np.linspace(-0.5,0.5,200)
        
        n, bins = np.histogram(edmlist,b,normed=True)
        
        #plt.hist(edmlist,b,normed=True,histtype='step')
        centers = (0.5*(b[1:]+b[:-1]))
        pars, cov = curve_fit(lambda x, mu, sig : norm.pdf(x, loc=mu, scale=sig), centers, n, p0=[1,1])
        
        
        
        #p = n/100
        #yerr = np.sqrt(np.sqrt(100)*p*(1-p)*1.95)
        #plt.errorbar(centers, n, yerr, fmt='C0_')
        
        
        #xs = np.linspace(-0.5,0.5,200)
        
        
        
        #fit = [norm.pdf(i,pars[0],pars[1]) for i in xs]
        #plt.plot(xs,fit, linewidth = 2,label = r'$\mu = 1.45x10^{-5}, \sigma = 0.010$')
        #plt.xlabel('Vertical angle [rad]')
        #plt.ylabel('Count [arbitrary units]')
        #plt.legend()
        
               
        m = pars[0]
        e = np.sqrt(cov[0,0])
        #sig = pars[1]
        
        
        #m = np.mean(edmlist)
        #s = np.std(edmlist)/np.sqrt(n_muons)

        
        avAngle.append(m)
        spreadAngle.append(e)     
                 
        print("Plotting EDM for time", t)
     
    plt.figure(1)
    plt.scatter(times,counts,marker='.',label='')    
    plt.xlabel('Time [ns]')
    plt.ylabel('Number of positrons with E>2000 MeV')
    #plt.legend()
    
    f=open('GM2.txt','w')
    for time,count in zip(times,counts):
        f.write(str(time)+','+str(count)+'\n')
    f.close()    
        
        
    plt.figure(2)
    plt.errorbar(times,avAngle,spreadAngle, marker='.',label='')

    
    plt.xlabel('Time [ns]')
    plt.ylabel('Average vertical angle')
    
    
    
    
    
    f=open('EDM.txt','w')
    for time,angle,err in zip(times,avAngle,spreadAngle):
        f.write(str(time)+'   '+str(angle)+'   '+ str(err)+'\n')
    f.close()
    
else:

    energies = []
    pangles = []
    momenta = []
    vangles = []
    olda = []
    othervan = []
    ov = []
    
    GM2_fixedangle = 0
    EDM_fixedangle = 1
    
    genmuons = decayMuons(n_muons,GM2_fixedangle,EDM_fixedangle)
            
    
    #get information out of all muon decays
    for i in genmuons:
        
        vangles.append(i.vangle)
        ov.append(i.oV)
        
        
        
        if i.E > 0:     
    
            energies.append(i.E)
            pangles.append(i.angle)
            momenta.append(i.P)
            olda.append(i.oA)
            vangles.append(i.vangle)
            ov.append(i.oV)
            print ('Muon', genmuons.index(i), ',positron energy ', i.E, 'MeV, momentum', i.P, 'and angle', i.angle)  
            

        
        
    
                
            
    #These plots can be used to check the samples follow the required distributions
    #Not needed unless you've tampered with those
    '''   
    #define checking plots for distributions
    es = range(0,54,1)
    espec = [energyDecaySpectrum(i) for i in es]
    especN = [i*1/float(spint.simps(espec,dx=1)) for i in espec]
    
    angles = np.arange(-np.pi,np.pi,0.1)
    dist = [1 + (1/3)*np.cos(theta) for theta in angles]
    distN = [i*1/float(spint.simps(dist,angles)) for i in dist]
    '''
    
    '''
    #plot positron energy spectrum
    plt.figure(1)
    plt.hist(energies,bins=50,histtype='step',normed=True)
    #plt.scatter(es,especN,marker='.')
    plt.xlabel('Positron energy [MeV]')
    plt.ylabel('Count [arbitrary units]')
    #plt.ylim(0,0.05)
    
    #plot positron angle spectrum
    plt.figure(2)
    #plt.hist(pangles,bins=50,histtype='step',label='Boosted')
    plt.hist((olda),bins=50,histtype='step',label='Unboosted')
    #plt.scatter(angles,distN,marker='.')
    plt.xlabel('cos(angle), pre boost')
    plt.ylabel('Count [arbitrary units]')
    #plt.legend()
    #plt.ylim(0,0.25)
    
    #plot momentum spectrum
    plt.figure(3)
    plt.hist(momenta,bins=50,histtype='step',normed=True)
    plt.xlabel('Positron momentum [MeV/c]')
    plt.ylabel('Count [arbitrary units]')
    
    
    #plot 2d angle vs energy 
    plt.figure(4)
    plt.hist2d(pangles, energies, (50, 50), cmap=plt.cm.jet)
    plt.xlabel('Angle from momentum direction')
    plt.ylabel('Positron energy [MeV]')
    plt.colorbar()
    #plt.xlim(0.925,1.)
    '''

    #plot vertical angle distribution
    plt.figure(5)
    bins = np.linspace(-0.5,0.5,200)
    
    plt.hist(vangles,bins=bins,histtype='step',normed=True)
    
    #plt.xlim(-np.pi/2,np.pi/2)
    
    n, bins = np.histogram(vangles,bins,normed=True)
    
    centers = (0.5*(bins[1:]+bins[:-1]))
    pars, cov = curve_fit(lambda x, mu, sig : norm.pdf(x, loc=mu, scale=sig), centers, n, p0=[1,1])
    
    print(pars[0],pars[1])
    
    x = np.linspace(-0.5,0.5,200)
    
    plt.plot(x,norm.pdf(x,pars[0],pars[1]))

    
    
    '''
    plt.hist(othervan, bins=600,histtype='step',label=r'$\theta_{V} \leq \pi/4$')
    plt.hist(vangles, bins=50,histtype='step',label=r'$\theta_{V} > \pi/4$')
    
    #plt.title(r'$lalala > \beta$')
    
    tsum = vangles+othervan
    plt.hist(tsum,bins=600, histtype='step',label='Sum')
    '''
    '''
    xt = plt.xticks()[0]  
    xmin, xmax = min(xt), max(xt)  
    lnspc = np.linspace(xmin, xmax, len(vangles))

    
    m, s = stats.norm.fit(vangles) 
    pdf_g = stats.norm.pdf(lnspc, m, s)   
    plt.plot(lnspc, pdf_g, label="Norm") 
    
    #plt.plot(lnspc,stats.norm.pdf(lnspc,np.mean(vangles), np.std(vangles)))
    
    #plt.hist(vangles,bins=50,histtype='step')
    '''
    plt.xlim(-0.5,0.5)
    plt.xlabel('Vertical angle [rad]')
    plt.ylabel('Count [arbitrary units]')
    #plt.legend()
    '''
    xs = np.linspace(-1.5,1.5,1000)
    g1 = stats.norm.pdf(xs,-0.00012,0.03221)*0.21
    g2 = stats.norm.pdf(xs,0.00165,0.72951)
    plt.plot(xs,g1,label=r'$\mu = -0.00012, \sigma =0.03221$')
    plt.plot(xs,g2,label=r'$\mu = 0.00165, \sigma =0.72951$')
    plt.legend()
    '''
    #angles = np.linspace(-np.pi/2,np.pi/2,100)
    #checks = [0.5*verticalAngleSpectrum(theta) for theta in angles]
    #checksN = [i*1/float(spint.simps(checks,dx=1)) for i in checks]
   # u = np.random.normal(0,0.8,n_muons)
    #plt.hist(u, normed = True, histtype='step')
    
    
    #plt.plot(angles,checks)
    
    
   
    #plt.plot(x,z)
    
    
    #plt.show()


