#!/usr/bin/env
# -*- coding: utf-8 -*-

'''
TODO:
    -Test test test!
    -Figure out the energy bump at low E
    -Add in EDM oscillation and check it's visible
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
        self.vangle = vertical_angle()
        
        #unchanged self variables for debugging
        self.oE = self.E
        self.oA = self.angle
        self.oP = self.P
    
    def boost(self):
    #method to boost the muon into the detector frame    
        gamma = 29.3 #'magic' gamma value for g-2
        mu_momentum = 3.094e3 #3.094 GeV/c
        mu_v = -mu_momentum/(gamma*muonM)  
        
        
            
          
            
        #if self.angle < 0.:
            #A = A + 2*np.pi
        A = self.angle  
        #var = abs(self.angle)%(2*np.pi)
        eP = self.P
        
        #if var > np.pi:
            #var = var - 2*np.pi
            
            
            
        if abs(A) > np.pi/2:
            eP = -eP
            A = np.pi - self.angle  
            
            #print('backward,', eP*np.cos(A))
            boostedE = gamma*(self.E-mu_v*eP*np.cos(A))        
            boostedPz = gamma*(eP*np.cos(A)-mu_v*self.E)
            Px = eP*np.sin(A)
        
            boostedPTot = np.sqrt(boostedPz**2+Px**2)
        
            boostedAngle = np.arctan(Px/boostedPz) 
        elif abs(A) <= np.pi/2:
            eP = eP
            A = self.angle  
            
            #print('forward,', eP*np.cos(A))
            boostedE = gamma*(self.E-mu_v*eP*np.cos(A))        
            boostedPz = gamma*(eP*np.cos(A)-mu_v*self.E)
            Px = eP*np.sin(A)
        
            boostedPTot = np.sqrt(boostedPz**2+Px**2)
        
            boostedAngle = np.arctan(Px/boostedPz)
        else:
            print('o no')
            
            
        #if A < np.pi/2 or A > 3*np.pi/2:
            #eP=-eP

               
        
        self.E = boostedE
        self.angle = boostedAngle
        self.lifetime = self.lifetime*gamma
        self.P = boostedPTot
        
    def gm2Shift(self,t):
    #method to add in the g-2 oscillation by rotating spin vector by an angle
    #IMPORTANT: must do this before boosting!
        #w = 1.5e-3
        #shift = w*t
        #self.angle=self.angle+shift
        self.angle = self.angle+t
        

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
   

     
#--------------------------------------------------------    
#vertical dists
        
def vertical_angle():
    vA = np.random.normal(0,0.1)
    return vA    
    

#--------------------------------------------------------


def decayMuons(n,t):
    muonlist = []
    
    for i in range(n+1):
        m = muon()
        m.gm2Shift(t)
        m.boost()        
        muonlist.append(m)
        #print ('Muon', i,', positron energy ', m.E, 'MeV, momentum', m.P, 'and angle', m.angle)
        #if m.E < 50:
            #print('E before:', m.oE, 'angle before', m.oA, 'momentum before', m.oP)

    
    return muonlist
   
#--------------------
energies = []
pangles = []
momenta = []
vangles = []
'''

npoints = 100
times = list(np.linspace(0., 30e3, npoints))
#counts = [2815, 2805, 2794, 2611, 2439, 2373, 2212, 2030, 1861, 1631, 1517, 1354, 1197, 1094, 1048, 962, 889, 970, 1029, 1175, 1291, 1425, 1577, 1850, 1953, 2115, 2312, 2406, 2525, 2629, 2774, 2782, 2804, 2738, 2703, 2507, 2430, 2217, 2097, 1928, 1754, 1602, 1375, 1294, 1164, 1040, 963, 975, 936, 1027, 1068, 1186, 1324, 1550, 1698, 1890, 2001, 2242, 2380, 2458, 2665, 2764, 2782, 2823, 2766, 2715, 2554, 2465, 2269, 2102, 1979, 1783, 1586, 1531, 1394, 1233, 1053, 1033, 995, 957, 977, 1071, 1186, 1314, 1440, 1570, 1799, 2028, 2128, 2240, 2469, 2542, 2667, 2794, 2724, 2802, 2670, 2613, 2551, 2372, 2245, 2004, 1834, 1661, 1535, 1341, 1299, 1106, 1045, 1050, 936, 939, 1039, 1048, 1291, 1383, 1532, 1776, 1935, 2037, 2273, 2395, 2545, 2706, 2674, 2812, 2725, 2745, 2635, 2625, 2540, 2326, 2202, 1968, 1813, 1635, 1433, 1327, 1215, 1080, 956, 982, 1013, 995, 1073, 1186, 1261, 1474, 1676, 1816]

counts = []
for t in times:
    counter = 0
    genmuons = decayMuons(n_muons,t)
    print('Generating muons for time', t)
    
    for i in genmuons:
        if i.E> 2000:
            counter += 1
    
    counts.append(counter)
    print(counter)
 



plt.figure(1)
plt.scatter(times,counts,marker='.',label='')

plt.xlabel('Time [ns]')
plt.ylabel('Number of positrons with E>2000 MeV')
#plt.legend()


'''
olda = []

testangle = np.pi

genmuons = decayMuons(n_muons,testangle)

#get information out of all muon decays
for i in genmuons:
    
    if i.E > 0:
        

        energies.append(i.E)
        pangles.append(i.angle)
        momenta.append(i.P)
        olda.append(i.oA)
        vangles.append(i.vangle)
        #print('Energy,angle,momentum:',i.oE,i.oA,i.oP,'//',i.E,i.angle,i.P)
            
        
    
    
    

        
        


#define checking plots for distributions
es = range(0,54,1)
espec = [energyDecaySpectrum(i) for i in es]
especN = [i*1/float(spint.simps(espec,dx=1)) for i in espec]

angles = np.arange(-np.pi,np.pi,0.1)
dist = [1 + (1/3)*np.cos(theta) for theta in angles]
distN = [i*1/float(spint.simps(dist,angles)) for i in dist]



#plot positron energy spectrum
plt.figure(1)
plt.hist(energies,bins=50,histtype='step',normed=True,label='Forward positrons, aparr')
#plt.scatter(es,especN,marker='.')
plt.xlabel('Positron energy [MeV]')
plt.ylabel('Count [arbitrary units]')
plt.legend()
#plt.ylim(0,0.05)
'''
#plot positron angle spectrum
plt.figure(2)
plt.hist(pangles,bins=50,histtype='step',label='Boosted')
plt.hist(olda,bins=50,histtype='step',label='Unboosted')
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
plt.figure(4)
plt.hist2d(np.cos(pangles), energies, (50, 50), cmap=plt.cm.jet)
plt.xlabel('Angle from momentum direction')
plt.ylabel('Positron energy [MeV]')

plt.figure(5)
plt.hist(vangles,bins=50,histtype='step')
'''


plt.show()


        