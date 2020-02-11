#!/usr/bin/env
# -*- coding: utf-8 -*-

'''
TODO:
    - add a detector and time tracking
    - look for wiggles
'''

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spint

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
        self.oE = self.E
        self.oA = self.angle
    
    def boost(self):
    #method to boost the muon into the detector frame    
        gamma = 29.3 #'magic' gamma value for g-2
        mu_momentum = 3.094e3 #3.094 GeV/c
        mu_v = -mu_momentum/(gamma*muonM)
        
        vPz = mu_v*self.P*np.cos(self.angle)
        
        boostedE = gamma*(self.E-vPz)       
       
        
        boostedPz = gamma*(self.P*np.cos(self.angle)-mu_v*self.E)
        Px = self.P*np.sin(self.angle)
        
        boostedPTot = np.sqrt(boostedPz**2+Px**2)
        
        boostedAngle = np.arctan(Px/boostedPz)
        
        self.E = boostedE
        self.angle = boostedAngle
        self.lifetime = self.lifetime*gamma
        self.P = boostedPTot
        
    def gm2Shift(self,t):
        w = 0.1
        shift = w*t
        #update this to work with a t point and w
        self.angle=self.angle+shift

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
        z = np.random.uniform(0, Emax)
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


    


''' 
def gm2Oscillation(angle):
    #!!! NON-REL ATM, CORRECT THIS!!!
    lifetime = 2.1969811e-6 #s
    omega = 1.5e-3 #radians per second    
    spinDirection = lifetime*omega    
    return angle + spinDirection
'''       
#--------------------------------------------------------    
#lifetime dist

#--------------------------------------------------------


def decayMuons(n,t):
    muonlist = []
    
    for i in range(n+1):
        m = muon()
        m.gm2Shift(t)
        m.boost()        
        muonlist.append(m)
        print ('Muon', i, ', positron energy ', m.E, 'MeV, momentum', m.P, 'and angle', m.angle)
    
    return muonlist
   

genmuons = decayMuons(n_muons,0.)

#--------------------
energies = []
pangles = []
momenta = []



#get information out of all muon decays
for i in genmuons:
    if i.E > 0:
        energies.append(i.E)
        pangles.append(i.angle)
        momenta.append(i.P)
        


#define checking plots for distributions
es = range(0,54,1)
espec = [energyDecaySpectrum(i) for i in es]
especN = [i*1/float(spint.simps(espec,dx=1)) for i in espec]

angles = np.arange(-np.pi,np.pi,0.1)
dist = [1 + (1/3)*np.cos(theta) for theta in angles]
distN = [i*1/float(spint.simps(dist,angles)) for i in dist]



#plot positron energy spectrum
plt.figure(1)
plt.hist(energies,bins=50,histtype='step',normed=True,label='')
#plt.scatter(es,especN,marker='.')
plt.xlabel('Positron energy [MeV]')
plt.ylabel('Count [arbitrary units]')
plt.legend()
#plt.ylim(0,0.05)

#plot positron angle spectrum
plt.figure(2)
plt.hist(pangles,bins=50,histtype='step')
#plt.scatter(angles,distN,marker='.')
plt.xlabel('Angle from momentum direction')
plt.ylabel('Count [arbitrary units]')
plt.legend()
#plt.ylim(0,0.25)

plt.figure(3)
plt.hist(momenta,bins=50,histtype='step',normed=True)
plt.xlabel('Positron momentum [MeV/c]')
plt.ylabel('Count [arbitrary units]')


#plot 2d angle vs energy 
plt.figure(3)

plt.hist2d(pangles, energies, (50, 50), cmap=plt.cm.jet)
plt.xlabel('Angle from momentum direction')
plt.ylabel('Positron energy [MeV]')



plt.show()


        