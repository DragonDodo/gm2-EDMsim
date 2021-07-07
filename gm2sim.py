#!/usr/bin/env
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
#import random
#import scipy.integrate as spint
#from scipy.stats import norm, skewnorm
#from scipy.optimize import curve_fit
import time as timestamp
#import ring_decay
from scipy.interpolate import interp1d

#make plots in g-2 style
plt.style.use('gm2.mplstyle')
plt.ion()

#--------------------------
#Options area

n_events = 100000#number of events to generate
t_start = 0 #start time (ns)
t_end = 30000 #end time (ns)
nbins =  100 #number of time bins to use

#Crude 'mode' selection: 
#Options:
#   'MCgen' to generate muon decays across the time range and plot oscillations
#   'test' to check overall distribtuions, not really needed any more, just add plots to MCgen
option = "MCgen"
#--------------------------

#MeV units
muonM = 105.66
eM = 0.51 
Emax = muonM/2
gamma = 29.3 #'magic' gamma value for g-2
ringRadius = 7.1 #m

class muon:
    #class to define a muon decay 
    #It's more definining the positron attributes really
    
    def __init__(self):
        #Generate muon attributes in the MRF
        eE = sample_muE(1)
        self.E = eE
        self.angle=sample_angle(1,eE)
        self.lifetime = 2.1969811e-6*1e9 #ns
        self.P = np.sqrt(self.E**2-eM**2)
        self.vangle = sample_vertical_angle(1)        
        self.decay_time = np.random.exponential(self.lifetime*gamma)
        #self.decay_time = np.random.uniform(t_start, t_end)
        
        ########### beam distribtuons ############
        
        self.beamx = sample_beam_dist_x(1)
        self.beamy = sample_beam_dist_y(1)
        
        self.pathlength = np.NaN #need boosted P to calc this
        
        #these attributes are for the ring_decay MC module 
        self.ringAngle = np.random.uniform(0,np.pi*2)
        self.x = np.cos(self.ringAngle)*ringRadius
        self.y = np.sin(self.ringAngle)*ringRadius 
        self.xhit = 0
        self.yhit = 0
        self.hitTracker = False
        self.nhits = 0            
        self.hits = []
        
        #########################################
        
        #unchanged self variables for debugging
        self.oE = self.E
        self.oA = self.angle
        self.oP = self.P
        self.oV = self.vangle
        self.oy = self.beamy
        self.ox = self.beamx
        
    
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
        #self.oA = np.arctan(Px/Pz)
        
        self.vangle = np.arctan(Py/boostedPz)
        #print("----")
        #print(self.P)
        #print(self.vangle)
        #print(Py/boostedPz)
        #print(np.arcsin(Py/self.P))
        #print("----")
        self.lifetime = self.lifetime*gamma
        
        self.pathlength = sample_pathlength(self.P)
        
        
        
    def gm2Shift(self,t):
    #method to add in the g-2 oscillation by rotating spin vector by an angle
    #IMPORTANT: must do this before boosting!
    #Can also lock the polarisation direction here in the plane 
    
        w = 1.5e-3
        shift = w*t
        self.angle=self.angle+shift
        #self.angle=self.angle+np.pi/2 #lock the 'polarisation'
        #self.oA = self.angle
        
    def EDMShift(self,t):
    #method to add in an EDM oscillation, ie a vertical angle oscillation from a tilted plane
    #Do before boosting!
        
        wEDM = 1.5e-3 # g-2 freq
        
        d_u = 100*1.9e-19 #100BNL
        #d_u = 1.76e-19
        consts = 9.158e15
        
        A = consts*d_u
        #A = 0
        #A = 10000/1e6
        #A = 0.17
        #A = np.arctan(10000/1e6)
        #A = 1000
        phase = np.pi/2# EDM is pi/2 out of phase with g-2 oscillation        
        dphi = A*np.cos(wEDM*t-phase) 
        

        theta  = self.angle
        phi  = self.vangle

        tangle = (np.sin(theta)/(np.cos(theta)*np.cos(dphi)-np.tan(phi)*np.sin(dphi)))
        svangle = np.cos(np.arctan(tangle))*np.cos(phi)*np.sin(dphi)+np.sin(phi)*np.cos(dphi)
        
        
        self.vangle = np.arcsin(svangle)
        
    
        self.oV = self.vangle
        

        
        
    def smear(self,uncert):
        #method to simulate the tracker resolution by smearing with a gaussian
        #it's probably too simple
        noise = np.random.normal(0,uncert) #assume gaussian uncert
        self.vangle = self.vangle + noise
        
        
    #def radial_field:
    #radial field produces a fake EDM by tilting the plane towards the center? of the ring 
    #need to fake the circular geometry - use the momentum direction? 
    
    def radial_field(self,t,Br):
        #method to add in the impact of a radial field 
        
        #do need to add in a shift here for the error
        
        w = 1.5e-3 # g-2 freq    

        
        #Br_sample = abs(np.random.normal(0,Br))
        
        A = Br/(1e6)
        #A = 0.1
        
        
       # A = 0

        phase = np.pi/2
        dphi = A*np.cos(w*t-phase) 
        

        theta  = self.angle
        phi  = -self.vangle

        tangle = np.sin(theta)/(np.cos(theta)*np.cos(dphi)-np.tan(phi)*np.sin(dphi))
        svangle = np.cos(np.arctan(tangle))*np.cos(phi)*np.sin(dphi)+np.sin(phi)*np.cos(dphi)

        self.vangle = np.arcsin(svangle)
        #self.vangle=np.arctan(np.tan(self.vangle)/gamma)
    
        self.oV = self.vangle

    def longitudinal_field(self,t):
        #method to add in the longitudinal field, ie a vertical angle oscillation IN PHASE with wa
    
        
        wEDM = 1.5e-3 # g-2 freq
        
        #A = 0.017 #1700ppm
        A = 0.17e-3
        #A = 10
        #A = 1.16e-2
        phase = 0# in phase with g-2 oscillation        
        dphi = A*np.cos(wEDM*t-phase) 
        

        theta  = self.angle
        phi  = -self.vangle

        tangle = (np.sin(theta)/(np.cos(theta)*np.cos(dphi)-np.tan(phi)*np.sin(dphi)))
        svangle = np.cos(np.arctan(tangle))*np.cos(phi)*np.sin(dphi)+np.sin(phi)*np.cos(dphi)
        

        self.vangle = np.arcsin(svangle)
        
    
        self.oV = self.vangle
        
    def shift_beam(self):
        #General method for shifting beam relative to trackers 
        #Have used this for beam dist, misaligment
        
        phi = self.vangle
        length = self.pathlength
        
        #rc = ring_decay.r_c(self.P)
        
        yphi = length*np.tan(phi)*1000 #convert to mm
        
       # xtheta = 1000*4*rc*np.sin(length/(2*rc))*np.sin(self.angle/2)
        internal_align_shift = 0.2 #mm
        
        #tracker = random.randint(1, 8)
        
        #if tracker < 7:
         #   self.beamy =  self.beamy + yphi + tracker*internal_align_shift
        #elif tracker == 8:
        #    self.beamy =  self.beamy + yphi - internal_align_shift
        
       # if tracker < 6:
        #    self.beamy =  self.beamy + yphi - internal_align_shift
        #elif tracker > 5:
        #    self.beamy =  self.beamy + yphi + internal_align_shift
            
        
        
        
        self.beamy =  self.beamy + yphi
        
        #self.beamx = xtheta + self.beamx



        
        
        
#-------------------------------------------------------
#Define the energy distribution of muon decay and sample

def energyDecaySpectrum(E):
    #muon decay energy spectrum, assuming massless positrons  
    #E = muon energy
    
    prefactor = 1e-5 #scale it close to 1
    bracket = (E**2)*np.float128((1-(4*E)/float(3*muonM)))
   
    return prefactor*bracket


def sample_muE(n): 
    #using rejection sampling to get a sample of size n
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
    #E = positron energy, as angle is correlated with energy
    #theta = precession angle in MRF
    #Pu is the polarisation, set to 1 as default 
    
    x = E/Emax

    a = (Pu*np.cos(theta))*(1-2*x)
    b = 3-2*x
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
        u = np.random.uniform(0, 12.)

        if u <= verticalAngleSpectrum(z):
            samples.append(z)

    if len(samples) ==1:
        return samples[0]
    else:
        return samples

#--------------------------------------------------------
        
def sample_beam_dist_x(x):
    #something that should be a function of momentum...
 
    xpos = skewnorm.rvs(a=-6,loc = 30.7, scale=25)
    ran = np.random.uniform(0,1)
    
    if ran > 0.85:
        xpos = np.random.normal(-18,10)
    return xpos
    
def sample_beam_dist_y(n):
    ypos = np.random.normal(0,12.9)
    return ypos

def sample_pathlength(p):
    #put the momentum dependence here in bins 
    #this is 2D MC
    #cx = [0,500,1000,1500,2000,2500,2700,2850,3200]
    #cy = [0.7,1.1,1.5,2.1,3,4.5,5.7,8.5,30]
    
    #this is data p-arc 
    cxd = [0,500,1000,1500,2000,2500,3000,3200]
    cyd = [1.,1.3,1.7,2.1,2.8,3.8,6.1,8]
    
    fit = interp1d(cxd, cyd,kind='linear')
    
   
    #fit = interp1d(cx, cy,kind='cubic')
    
    #cofs = [1.15217966e-18,-8.46978380e-15,2.35601047e-11,-3.01736447e-08,1.73910224e-05,-1.42723959e-03,1.39080843e-03]
    
    #fit = np.poly1d(cofs)
    
    
    
    mean = float(fit(p))
    
    return abs(np.random.normal(loc=mean,scale=0.5))
    #return 0.8
  
#--------------------------------------------------------

def decayMuons(n):
    #main muon generating function
    #all functions act on the muon object
    #n = number of muons to generte
    muonlist = []
    
    
    for i in range(n+1):
        m = muon()
        m.EDMShift(m.decay_time) #apply EDM precession
        m.gm2Shift(m.decay_time) #apply g-2 precession   
        #m.longitudinal_field(m.decay_time)
        m.boost() #boost into lab frame 
        #m.longitudinal_field(m.decay_time)
        
        #m.radial_field(m.decay_time,10000)  
        #m.shift_beam()
        
        #m.smear(0.05) #uniform tracker resolution
        muonlist.append(m)  
        pc = i/n*100
        if pc%10 == 0:
            print('Generating decay',i,'of',n,',',int(i/n*100),'% done...','\n', end='')
    return muonlist
   
#-------------------------------------------------------

def detectorAcceptanceCut(data,val_low,val_high):
    #simple vertical angle cut to simulate detector acceptance between val_low and val_high
    output = []
    
    if val_low == None and val_high == None:
        return data
    
    else:    
        for i in data:
            #if abs(i.beamy) <= val:
            #ran = np.random.uniform(300,1000)
            if i.beamy <= val_high and i.beamy >= val_low:
                output.append(i)
            
        #return mom_acc_cut(output)
        return output
    
def mom_acc_cut(data):
    output = acc_filter(data,momAcceptanceFunction,np.linspace(0,3089,50))
    
    return output
    
    
def momAcceptanceFunction(x):
    #gas gun
    binz = np.linspace(0,3100,10)
    binz = np.insert(binz,1,300)
    corr = np.array([ 0.        ,  0, 0.14418263,  0.6452906 ,  0.92512741,  0.99520434,
        0.97328187,  0.880436  ,  0.72599633,  0.62138323,  0.        ]) 
    
    #fit = interp1d(binz, corr,kind='linear')
    
    nhitbins = np. array([    0.        ,   110.34482759,   220.68965517,   331.03448276,
         441.37931034,   551.72413793,   662.06896552,   772.4137931 ,
         882.75862069,   993.10344828,  1103.44827586,  1213.79310345,
        1324.13793103,  1434.48275862,  1544.82758621,  1655.17241379,
        1765.51724138,  1875.86206897,  1986.20689655,  2096.55172414,
        2206.89655172,  2317.24137931,  2427.5862069 ,  2537.93103448,
        2648.27586207,  2758.62068966,  2868.96551724,  2979.31034483,
        3109.65517241])
    
    nhits = np.array([ 0.        ,  0.        ,  0.        ,  0.06416584,  0.29821261,
        0.42457588,  0.48974247,  0.55030891,  0.60542291,  0.61363636,
        0.6358164 ,  0.64919663,  0.64183486,  0.67966775,  0.69102875,
        0.6816218 ,  0.67947861,  0.68691296,  0.70479928,  0.70959378,
        0.74560922,  0.73001066,  0.73916274,  0.77342382,  0.80515527,
        0.81412201,  0.80924349,  0.79277267,  0.        ])
    
    
    nhitshighp = np.array([ 0.        ,  0.        ,  0.        ,  0.01030928,  0.05365854,
        0.18      ,  0.32110092,  0.34710744,  0.44104803,  0.49145299,
        0.6300813 ,  0.58736059,  0.63846154,  0.5974026 ,  0.63859649,
        0.64052288,  0.6722973 ,  0.61309524,  0.59593023,  0.59638554,
        0.56648936,  0.47804878,  0.4354067 ,  0.42380952,  0.34285714,
        0.25297114,  0.17615894,  0.04780362,  0.        ])
    
    fit = interp1d(nhitbins, nhitshighp,kind='linear')
    
    
    
    #MC nhits
     
    
    return float(fit(x))

def acc_filter(data,function,xs):
    #function to perform cuts that are a function of their variable
    
    scan = [function(i) for i in xs]    
    maxfn = max(scan)   
    #print(maxfn)
    
    samples = []
    
    for d in data:
        #print(d)
        u = np.random.uniform(0, maxfn)
    
        if u <= momAcceptanceFunction(d.P):
            samples.append(d)
    return samples

timestr = timestamp.strftime("%Y%m%d-%H%M%S")

#generate n muon decay events with previously defined time distribution
genmuons = decayMuons(n_events)


if option == "MCgen":

    bins = np.linspace(t_start,t_end,nbins)
    T = 2*np.pi/1.5e-3
    modbins = np.linspace(0,T,nbins)
    

    #define tracker boundaries, if using the ring MC
    x_max = 6.9
    x_min = 6.7
    y_min = -1
    y_max = 0

    ### Apply cuts to the data 
    data1 = genmuons
    data3 = detectorAcceptanceCut(genmuons,-150,150) 
    #data2 = mom_acc_cut(data1)
    #data1 = ring_decay.hit_tracker(genmuons,x_min,x_max,y_min,y_max)[0]


    datalist = [data1]
    
    cnt = 0
    
    minx = -50
    maxx = 50
    bins22 = np.linspace(minx,maxx,30)
    

    #print(eff)
    
    for dset in datalist: 

        high_E = []
        vangles = []
        times = []
        avAngle= []
        spreadAngle = []
        fit1 =[]
        fit2 = [] 
        olda = []
        moms = []
        momx = []
        momy =[]
        momz =[]

    
        for i in dset: #g-2 energy cut
            if i.E > 2000:
                high_E.append(i.decay_time)     
        
            #apply EDM momentum cuts
            ang = i.oA%(2*np.pi)
            
            if i.P > 0 and i.P < 5000: #and abs(i.oA) <= np.pi/2:  
                olda.append(i.oA)
                Px = i.P*np.sin(i.angle)*np.cos(i.vangle)
                Py = i.P*np.sin(i.vangle)
                Pz = i.P*np.cos(i.angle)*np.cos(i.vangle)
                moms.append(i.P)
                momx.append(Px)
                momy.append(Py)
                momz.append(Pz)
                modtime = i.decay_time%T
                times.append(modtime)
                vangles.append(i.vangle)   

            
        #bin the counts/angles into time bins
        counts, edges = np.histogram(high_E,bins)    
        uncert = np.sqrt(counts)
        
        digitized = np.digitize(times, modbins)
        vangles = np.array(vangles)
        
        for i in range(1,len(modbins)):
            bin_contents = vangles[digitized==i]    
            
            if len(bin_contents) == 0:
                stmean = 0
                sterr = 0
                
            else:
            
                stmean = np.mean(bin_contents)
                std = np.std(bin_contents)
            
                #fit gaussian for average vertical angle and error
                b = np.linspace(stmean-1*std,stmean+1*std,500)   
                n, fitbins = np.histogram(bin_contents,b,density=True)      
            
                centers = (0.5*(b[1:]+b[:-1]))
                #pars, cov = curve_fit(lambda x, mu, sig : norm.pdf(x, loc=mu, scale=sig), centers, n, p0=[1,1])
                   
                #m = pars[0]
                #e = np.sqrt(cov[0,0])    
                
                #stmean = np.mean(bin_contents)
                stderr = np.std(bin_contents)/np.sqrt(len(bin_contents))
            
            
            avAngle.append(stmean)
            spreadAngle.append(stderr)
            #fit1.append(m)
            #fit2.append(e)
            '''        
            plt.hist(bin_contents,b,histtype='step',normed=True)
            plt.xlabel('Phi angle post-boost')
            plt.plot(b,norm.pdf(b,m,pars[1]),label = 'fit')
            plt.legend()
            '''
        #avAngle = [i/10 for i in avAngle]
        #spreadAngle = [i/10 for i in spreadAngle]           
        avAngle = np.array(avAngle)    
        spreadAngle = np.array(spreadAngle)
        
        
        #Plot wiggles for checking outputs: not nice plots, use the plotting code!
        plt.figure(1)    
        plt.fill_between(bins[:-1],counts-uncert,counts+uncert,alpha=0.7,color='#ff7f0e')
        plt.scatter(bins[:-1],counts,marker='.',label='',color='k')       
        plt.xlabel("Time [ns]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
        plt.ylabel("Count with E>2000 MeV/300ns",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0, labelpad = 20)
        #plt.legend()    
        ###33cccc
        ##942ed2 #purple
        ##668edd blue
        
    


        
        plt.figure(2)
        plt.fill_between(modbins[:-1],avAngle-spreadAngle,avAngle+spreadAngle,alpha=0.7,color='r',label='40-50 MeV')
        plt.scatter(modbins[:-1],avAngle, marker='.',color='k')     
        plt.xlim(0,max(modbins))
        
        
        plt.xlabel("Time [ns]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
        plt.ylabel("Average vertical angle [rad] /300ns",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0, labelpad = 20)     
        #plt.legend()
        
        name_stampG = str('GM2-')+str(timestr)
        name_stampE = str('delEDM-')+str(timestr)
        
        ns = str('EDM-acc')+str(cnt)+str('.txt')    
    
    
        
        #f=open('%s.txt' %name_stampG,'w')
        f = open('gm2.txt','w')
        for time,count,err in zip(bins[:-1],counts,uncert):
            f.write(str(time)+','+str(count)+','+str(err)+'\n')
        f.close()  
        
        n = 'edmcheck'
        
        f=open('mrf.txt','w')
        for time,angle,err in zip(modbins[:-1],avAngle,spreadAngle):
            f.write(str(time)+','+str(angle)+','+ str(err)+'\n')
        f.close()
        cnt += 1
        
        mombins = np.arange(-100,100,10)
        
        plt.figure(3)
        plt.hist(momy,bins=mombins,histtype='step',density=True, label = 'Unboosted')
        plt.xlabel("Positron total momentum [MeV]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
        plt.ylabel("Count/bin [arbitrary units]",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0, labelpad = 20)
        plt.ylim(bottom=0)
        
        print(np.std(vangles))
        
        
        

       
if option == "test":
    
    #this area is more of a sandbox, good for checking distributions
    #it plots things but doesn't save any data 


    vangles = []
    ypos = []
    xpos = []
    
    oldy = []
    
    ooy = []
    arc = []
    
    mom = []
    path = []
    oldx = []
    
    th = []


    vpre = []
    vforward = []
    vbackward= []
    
    T = (2*np.pi/1.5e-3)-np.pi/2      
    precut = []
    for i in genmuons:
        if i.decay_time > 5000 and i.decay_time > 4000: 
            precut.append(i.P)
            ooy.append(i.oy)
            th.append(i.angle)
            if abs(i.oA) <= np.pi/2:
                vforward.append(i.vangle)
            if abs(i.oA) > np.pi/2:
                vbackward.append(i.vangle)
            vpre.append(i.vangle)
        
        
    x_max = 6.9
    x_min = 6.7
    y_min = -1
    y_max = 0
    

    #data = ring_decay.hit_tracker(genmuons,x_min,x_max,y_min,y_max)[0]
    data = detectorAcceptanceCut(genmuons,50) 
    #data = genmuons
    
    #get information out of all muon decays
    for i in data:    
        if i.beamy > -30 and i.beamy < -20 and i.vangle>-1000:

            vangles.append(i.vangle)
            ypos.append(i.beamy)
            xpos.append(i.beamx)
            oldy.append(i.oy)
            
            
            mom.append(i.P)
            path.append(i.pathlength)
            oldx.append(i.ox)
        

        
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
        
    '''            
    #These plots can be used to check the samples follow the required distributions
    #Not needed unless you've tampered with those
      
    #define checking plots for distributions
    es = range(0,54,1)
    espec = [energyDecaySpectrum(i) for i in es]
    especN = [i*1/float(spint.simps(espec,dx=1)) for i in espec]
    
    angles = np.arange(-np.pi,np.pi,0.1)
    dist = [1 + (1/3)*np.cos(theta) for theta in angles]
    distN = [i*1/float(spint.simps(dist,angles)) for i in dist]
    
    #plot positron energy spectrum
    plt.figure(1)
    plt.hist(olde,bins=50,histtype='step',normed=True, label = 'Unboosted')
    plt.hist(energies,bins=50,histtype='step',normed=True, label = 'Boosted')
    #plt.scatter(es,especN,marker='.', label = 'Expected functional form')
    plt.xlabel("Positron energy [MeV]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
    plt.ylabel("Count/bin [arbitrary units]",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0, labelpad = 20)
    plt.ylim(bottom=0)
        
    #plot positron angle spectrum
    plt.figure(2)
    #plt.hist(pangles,bins=50,histtype='step',label='Boosted')
    plt.hist((olda),bins=50,histtype='step',label='MC data',normed=True)
    #plt.scatter(angles,distN,marker='.', label = 'Expected functional form')
    plt.xlabel('cos(angle), pre boost')
    plt.ylabel('Count [arbitrary units]')
    plt.xlabel(r"Azimuthal angle $\theta$ [rad]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
    plt.ylabel("Count/bin [arbitrary units]",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0, labelpad = 20)
    plt.ylim(bottom=0)
    plt.legend(loc='upper left')
    plt.ylim(0,0.30)
    
    
    #plot momentum spectrum
    plt.figure(3)
    plt.hist(momenta,bins=50,histtype='step',normed=True)
    plt.xlabel('Positron momentum [MeV/c]')
    plt.ylabel('Count [arbitrary units]')
    
    
    #plot 2d angle vs energy 
    plt.figure(4)
    plt.hist2d(olde, energies, (50, 50), cmap=plt.cm.jet)
    plt.colorbar()
    plt.xlabel(r"Angle from momentum direction [rad]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
    plt.ylabel("Positron energy [MeV]",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0, labelpad = 20)
    #plt.xlim(0.925,1.)
    

    
    #plot vertical angle distribution
    plt.figure(1)
    minx = -40
    maxx = 40
    bins = np.linspace(minx,maxx,30)

    plt.hist(precut,bins=bins,histtype='step',label="No cut")
    plt.hist(ypos,bins=bins,histtype='step',label="With cut")
    plt.xlabel("Position [mm]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
    plt.ylabel("Counts/bin",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=1.0, labelpad=20)
    
    binsv = np.linspace(-100e-3,100e-3,50)
    
    pre = np.histogram(vpre,binsv)[0]  
    post = np.histogram(vangles,binsv)[0]
    
    eff = post/pre
    
    ma = max(eff)
    
    eff = np.array(eff)/ma
    
    eff = np.array([ 0.        ,  0.        ,  0.04458911,  0.        ,  0.03369712,
        0.07357204,  0.        ,  0.07524413,  0.06588541,  0.03618297,
        0.0728436 ,  0.05714333,  0.19939965,  0.25928977,  0.26404576,
        0.29974309,  0.41034023,  0.50267478,  0.67688351,  0.79258451,
        1.        ,  1, 0.99106955,  0.99480762,  0.85602642,
        0.86721218,  0.83254371,  0.71364973,  0.5275036,  0.1965231,
        0.0806269,  0.01        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.        ,  0.        ,  0.        ,  0.        ])

    
    plt.scatter(binsv[:-1],eff,label='This MC',color='C0')
    
    #gas gun comparsions

    #effjv = [0.74,0.84,0.918,0.96,0.98,0.98,0.96,0.918,0.84,0.74]
    #effjv = [i*1.02 for i in effjv]
    #binstot = np.linspace(-40e-3,40e-3,10)
    
    
    #fit = interp1d(binsv, effjv,kind='cubic')

    mbins = np.linspace(0,3070,100)
    #js = [fit(i) for i in mbins]
    #plt.plot(mbins,js,color='C0',label='Gas Gun')
    

    
    effjv = [0,
             0,
             0.02,
             0.09,
             0.18,
             0.31,
             0.6,
             0.84,
             1,
             1,
             0.95,
             0.51,
             0.05,
             0,
             0,
             0,
             0,
             0,
             0,
             0           

                            
            
            ]
    
    binz = np.linspace(-100e-3,100e-3,20)
    effj = np.array([0,0.2,0.7,0.9,0.95,0.91,0.8,0.7,0.6,0])
    
    mj = max(effjv)
    effjv = [i/mj for i in effjv]    
    
    plt.scatter(binz,effjv, color = 'C1', label = 'Gas Gun')
    plt.legend()
    

    
    plt.xlabel("Vertical angle [rad]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
    plt.ylabel("Efficiency",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=1.0, labelpad=20)
    
    
    binsv = np.linspace(-100e-3,100e-3,50)
    fig, ax1 = plt.subplots()
    factor = 0.5
    plt.figure(2)
    plt.xlabel("Vertical angle [rad]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
    plt.ylabel("Counts/bin",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=1.0, labelpad=20)
    ax1.tick_params(axis='y', labelcolor='k')
    ax2 = ax1.twinx()
    ax2.tick_params(labelcolor='w')
    ax1.hist(vpre, binsv, histtype='step', label = 'All positrons',density=True, color='C0')

    ax2.hist(vangles, binsv, histtype='step',weights=factor*np.ones_like(vangles), label = 'Passing acceptance',density=True, color = 'C1')
    
    
    
 
    ax1.hist(a, histtype='step', color='b', lw=2, density=True)
    ax1.tick_params(axis='y', labelcolor='b')
    ax2 = ax1.twinx()
    ax2.hist(b, histtype='step', color='r', lw=2, density=True)
    ax2.tick_params(axis='y', labelcolor='r')
    '''   
    '''
  
    plt.xlabel("Momentum [MeV]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
    plt.ylabel("Counts/bin",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=1.0, labelpad=20)
    
    plt.figure(3)
    bin2  = np.linspace(-100,100,50)
    plt.hist(ooy,bin2,histtype='step',label='Width = 20.0')
    plt.hist(xpos,bin2,histtype='step',label='Width = 20.0')
    plt.xlabel("Positron x position at trackers [mm]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
    plt.ylabel("Counts/bin",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=1.0, labelpad=20)
    
    mean = np.mean(xpos)
    std = np.std(xpos)
    
    plt.text(-100,10000,'Mean: %.2f, Std: %.2f' %(mean,std))

    bins = np.linspace(-80,80,50)
    
    plt.figure(4)
    my_cmap = plt.cm.jet
    my_cmap.set_under('w',1)
    plt.xlim(-80,80)
    plt.ylim(-80,80)
    plt.hist2d(xpos,ypos,bins=(bins,bins),cmap=my_cmap,vmin=1.)  
    plt.colorbar()  

        
    #plt.xlim(minx,maxx)
    #plt.ylim(0,21)
    
    plt.xlabel("Radial position [mm]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
    plt.ylabel("Vertical position [mm]",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=1.0, labelpad=20)
    
    #plt.xlabel("Beam x position [mm]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
    #plt.ylabel("Beam y position [mm]",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=1.0, labelpad=20)
    '''
    

    
    '''
    cx = [0,500,1000,1500,2000,2500,2700,2850,3200]
    cy = [0.7,1.1,1.5,2.1,3,4.5,5.7,8.5,30]
    cxd = [0,500,1000,1500,2000,2500,3000,3200]
    cyd = [1.,1.3,1.7,2.1,2.8,3.8,6.1,100]
    
    fit = interp1d(cxd, cyd,kind='linear')
    #p = np.poly1d(fit)
    
    #cofs = [1.15217966e-18,-8.46978380e-15,2.35601047e-11,-3.01736447e-08,1.73910224e-05,-1.42723959e-03,1.39080843e-03]
    
    #p = np.poly1d(cofs)
    

    
    
    
    xp = np.linspace(0,3200,100)
    #plt.scatter(mom,path,marker='.',color='k',label = 'Simulated data')
    
    plt.scatter(cxd,cyd,marker='o',color='r')
    plt.plot(xp,fit(xp),color='r')
    
    


    #bins = 50
    
      
    n, bins = np.histogram(oldv,bins,density=True)
    
    centers = (0.5*(bins[1:]+bins[:-1]))
    pars, cov = curve_fit(lambda x, mu, sig : norm.pdf(x, loc=mu, scale=sig), centers, n, p0=[1,1]) #can use cauchy
    
    print(pars[0],pars[1])
    
    err = np.sqrt(cov[0,0])
    
    x = np.linspace(minx,maxx,200)
    #plt.ylim(0,0.7)
    
    plt.hist(oldv,bins=bins,histtype='step',label="Unboosted",density=True)
    plt.hist(vangles,bins=bins,histtype='step',label="Boosted",density=True)
    #plt.scatter(bins,[0.5*verticalAngleSpectrum(b) for b in bins],marker='.', label = 'Expected functional form')
    
    
    plt.legend()
    

    #plt.plot(x,norm.pdf(x,pars[0],pars[1]))

    #plt.xlim(-0.5,0.5)
    plt.xlabel("Vertical angle [rad]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
    plt.ylabel("Counts/bin [arbitrary units]",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=1.0, labelpad=20)

    
    plt.figure(6)
    plt.hist(times,bins=50,histtype='step',density=True)
    plt.show()
    
    aa = np.linspace(-np.pi,np.pi,100)
    ens = [10,20,30,40,50]
    plt.figure(7)
    plt.xlabel("Angle from momentum direction [rad]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
    plt.ylabel("Expected counts/bin [arbitrary units]",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=1.0, labelpad=20)

    for i in ens:
        plt.plot(aa,angleSpectrum(i,aa),label='%i MeV' %i)
        
      
    Brad = 60
        
    Berr = 60
    
    cnt = np.linspace(0,100,n_events)
        
    sample = [np.random.normal(Brad,Berr) for i in cnt]
    
    plt.figure(7)
    plt.hist(sample,bins=50)
    '''
    
    
    