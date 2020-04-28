#!/usr/bin/env
# -*- coding: utf-8 -*-

'''
TODO:
    -...
'''

import numpy as np
import matplotlib.pyplot as plt
#import scipy.integrate as spint
from scipy.stats import norm
from scipy.optimize import curve_fit

n_events = 1000000
t_start = 0
t_end = 30000

#Crude 'mode' selection: 
#Options:
#   'EDM' to make vertical angle oscialltion plots
#   'GM2' to make usual wiggle plot
#   'test' for anything else. 
option = "EDM"

#MeV units
muonM = 105.66
eM = 0.51 
Emax = muonM/2
gamma = 29.3 #'magic' gamma value for g-2

class muon:
    #class to define a muon decay object
    def __init__(self):
        eE = sample_muE(1)
        self.E = eE
        self.angle=sample_angle(1,eE)
        self.lifetime = 2.1969811e-6 #seconds
        self.P = np.sqrt(self.E**2-eM**2)
        self.vangle = sample_vertical_angle(1)        
        self.decay_time = np.random.exponential(1/(self.lifetime*gamma))
        #self.decay_time = np.random.uniform(t_start, t_end)
        
        #unchanged self variables for debugging
        self.oE = self.E
        self.oA = self.angle
        self.oP = self.P
        self.oV = self.vangle
    
    def boost(self):
        
        #method to boost the muon into the detector frame  
        muonM = 105.66        
        mu_momentum = 3.094e3 #3.094 GeV/c
        mu_v = -mu_momentum/(gamma*muonM) #-tive because boosting back into lab frame
        
        Px = self.P*np.sin(self.angle)*np.cos(self.vangle)
        Py = self.P*np.sin(self.vangle)
        Pz = self.P*np.cos(self.angle)*np.cos(self.vangle)
        
        boostedPz = gamma*(Pz-mu_v*self.E)
        boostedE = gamma*(self.E-mu_v*Pz) 
        
        self.P = np.sqrt(boostedPz**2+Px**2+Py**2)
        self.E = boostedE
        self.angle = np.arctan(Px/boostedPz)
        self.vangle=np.arctan(np.tan(self.vangle)/gamma)
        self.lifetime = self.lifetime*gamma
        
        
    def gm2Shift(self,t):
    #method to add in the g-2 oscillation by rotating spin vector by an angle
    #IMPORTANT: must do this before boosting!
        w = 1.5e-3
        shift = w*t
        #shift = t
        self.angle=self.angle+shift
        self.oA = self.angle 
        
    def EDMShift(self,t):
    #method to add in an EDM oscillation, ie a vertical angle oscillation
    #Do before boosting!
        
        wEDM = 1.5e-3 # g-2 freq
        
        d_u = 100*1.9e-19 #100BNL
        consts = 8.9e15
        
        A = d_u*consts
        #A = 1
        phase = np.pi/2# EDM is pi/2 out of phase with g-2 oscillation        
        dphi = A*np.cos(wEDM*t-phase) 
        

        theta  = self.angle
        phi  = self.vangle

        tangle = np.sin(theta)/(np.cos(theta)*np.cos(dphi)-np.tan(phi)*np.sin(dphi))
        svangle = np.cos(theta)*np.cos(phi)*np.sin(dphi)+np.sin(phi)*np.cos(dphi)

        self.angle = np.arctan(tangle)
        self.vangle = np.arcsin(svangle)
        
    
        self.oV = self.vangle
        
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

def verticalAngleSpectrum(theta):   
    a = np.cos(theta)
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

def decayMuons(n):
    muonlist = []
    
    for i in range(n+1):
        m = muon()
        m.EDMShift(m.decay_time) #apply EDM precession
        m.gm2Shift(m.decay_time) #apply g-2 precession        
        m.boost() #boost into lab frame  
        
        muonlist.append(m)  
        #print('Generating decay',i,'of',n,',',int(i/n*100),'% done...', end='')
    return muonlist
   
#-------------------------------------------------------
    
genmuons = decayMuons(n_events)


T = 2*np.pi/1.5e-3

nbins = 100



if option == "GM2":  
    bins = np.linspace(t_start,t_end,nbins)
    
    high_E = []
    
    for i in genmuons:
        if i.E > 2000:
            high_E.append(i.decay_time)
            
    
    counts, edges = np.histogram(high_E,bins)    
    uncert = np.sqrt(counts)
    
    plt.figure(1)    
    plt.fill_between(bins[:-1],counts-uncert,counts+uncert,alpha=0.7,color='#ff7f0e')
    plt.scatter(bins[:-1],counts,marker='.',label='',color='k')    
    
    plt.xlabel('Time [ns]')
    plt.ylabel('Number of positrons with E>2000 MeV')
    plt.legend()
    
    f=open('GM2.txt','w')
    for time,count,err in zip(bins[:-1],counts,uncert):
        f.write(str(time)+','+str(count)+','+str(err)+'\n')
    f.close()
    

elif option == "EDM":
    bins = np.linspace(0,T,nbins)

    vangles = []
    times = []
    avAngle= []
    spreadAngle = []
    fit1 =[]
    fit2 = []   

    
    for i in genmuons:
        modtime = i.decay_time%T
        times.append(modtime)
        vangles.append(i.vangle)   
    
    digitized = np.digitize(times, bins)
    vangles = np.array(vangles)   

    for i in range(1,len(bins)):
        bin_contents = vangles[digitized==i]

        
        stmean = np.mean(bin_contents)
        std = np.std(bin_contents)
        
        b = np.linspace(stmean-1*std,stmean+1*std,500)   
        #b = np.linspace(-1.5,1.5,2000) 
        n, fitbins = np.histogram(bin_contents,b,density=True)      
        
        centers = (0.5*(b[1:]+b[:-1]))
        pars, cov = curve_fit(lambda x, mu, sig : norm.pdf(x, loc=mu, scale=sig), centers, n, p0=[1,1])
               
        m = pars[0]
        e = np.sqrt(cov[0,0])    
        
        stmean = np.mean(bin_contents)
        stderr = np.std(bin_contents)/np.sqrt(len(bin_contents))
        
        
        avAngle.append(stmean)
        spreadAngle.append(stderr)
        fit1.append(m)
        fit2.append(e)
        '''        
        plt.hist(bin_contents,b,histtype='step',normed=True)
        plt.xlabel('Phi angle post-boost')
        plt.plot(b,norm.pdf(b,m,pars[1]),label = 'fit')
        plt.legend()
        '''

    
    plt.figure(2)
    plt.errorbar(bins[:-1],avAngle,spreadAngle, marker='.')     
    plt.xlabel('Time [ns]')
    plt.ylabel('Average vertical angle')   
    
    
    f=open('EDM.txt','w')
    for time,angle,err in zip(bins[:-1],avAngle,spreadAngle):
        f.write(str(time)+','+str(angle)+','+ str(err)+'\n')
    f.close()
    

if option == "test":

    energies = []
    pangles = []
    momenta = []
    vangles = []
    olda = []
    times = []

    
    T = 2*np.pi/1.5e-3
    

            
    
    #get information out of all muon decays
    for i in genmuons:
        

        vangles.append(i.vangle)
        times.append(i.decay_time)
        

        
        '''
        #simple cuts on sim data
        if i.E > 0:     
    
            energies.append(i.E)
            pangles.append(i.angle)
            momenta.append(i.P)
            olda.append(i.oA)
            vangles.append(i.vangle)
            ov.append(i.oV)
            print ('Muon', genmuons.index(i), ',positron energy ', i.E, 'MeV, momentum', i.P, 'and angle', i.angle)  
            
        '''
        
                 
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
    minx = -0.5
    maxx = 0.5
    bins = np.linspace(minx,maxx,100)
    #bins = 50

      
    n, bins = np.histogram(vangles,bins,density=True)
    
    centers = (0.5*(bins[1:]+bins[:-1]))
    pars, cov = curve_fit(lambda x, mu, sig : norm.pdf(x, loc=mu, scale=sig), centers, n, p0=[1,1]) #can use cauchy
    
    print(pars[0],pars[1])
    
    err = np.sqrt(cov[0,0])
    
    x = np.linspace(minx,maxx,200)
    
    plt.hist(vangles,bins=bins,histtype='step',density=True,label="Using tan and gamma")
    plt.legend()
    

    #plt.plot(x,norm.pdf(x,pars[0],pars[1]))

    #plt.xlim(-0.5,0.5)
    plt.xlabel('Vertical angle [rad]')
    plt.ylabel('Count [arbitrary units]')
    #plt.legend()
    
    plt.figure(6)
    plt.hist(times,bins=50,histtype='step',density=True)

    plt.show()