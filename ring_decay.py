#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 10:38:12 2020

@author: dvasilkova
"""

import numpy as np  

def centre(x,y,r,theta):
    #x,y of decay point, find centre of tangent circle 
    N = np.sqrt(x**2+y**2)
    vect = (-x/N,-y/N)
    newx = x + r*vect[0]
    newy = y + r*vect[1]
    return rotate((x,y),(newx,newy),theta)


def track(x,y,r,theta):
    #calculate the circle of a track from a decay at x,y
    #x,y of tangent point, r radius of curvature
    angles = np.linspace(-np.pi,np.pi,1000)
    sx, sy = centre(x,y,r,theta)
    

    xp = np.cos(-angles)*r+sx
    yp = np.sin(-angles)*r+sy              
       
    return xp,yp,sx,sy

def r_c(P):
    #calculate radius of curvature in m, P = positron momentum in MeV
    B = 1.45 
    e = 0.303
    return (P*1e-3)/(B*e)


def hit_tracker(decays,trackerx_min,trackerx_max,trackery_min,trackery_max):
    #check if the track hits the tracker, defined to exist between xmin-xmax, ymin-ymax
    #decays should be a list of muon decay-style objects with attributes for ring position and momentum
    passed = []
    
    
    for d in decays:
        rad = r_c(d.P)
        circx, circy, cx, cy = track(d.x,d.y,rad,d.angle)
        
        hitcount = 0 
    
        for i in range(len(circx)):
            xpos = circx[i]
            ypos = circy[i]
            #print(xpos,ypos)                   
            
            if xpos > trackerx_min and xpos < trackerx_max and ypos < trackery_max and ypos > trackery_min:
                  #print('hit!')
                
                  d.nhits += 1
                  d.hits.append((xpos,ypos))
                  
                  
                  if hitcount == 0:                                
                  
                      passed.append(d)                           
                            
                      d.hitTracker = True
                    
                      xh = xpos
                      yh = ypos
                    
                      x1 = d.x
                      y1 = d.y
                    
                      a = np.sqrt((x1-xh)**2+(y1-yh)**2)
                      b = np.sqrt((x1-cx)**2+(y1-cy)**2)
                      c = np.sqrt((cx-xh)**2+(cy-yh)**2)
                    
                      angle = np.arccos((b**2+c**2-a**2)/(2*b*c))        
            
                    
                      d.arclength = angle*rad

                    
                      d.xhit = xh
                      d.yhit = yh
                      
                      hitcount = 1
                

            
    return passed, circx,circy, cx, cy


def rotate(origin, point, angle):
    #helper function to rotate tracks by the azimuthal decay angle
    ox, oy = origin
    px, py = point    
    
    qx = np.cos(angle)*(px-ox) - np.sin(angle)*(py-oy) + ox
    qy = np.sin(angle)*(px-ox) + np.cos(angle)*(py-oy) + oy

    return qx,qy



    
    
 
if __name__ == "__main__":   
    
    import matplotlib.pyplot as plt
    from scipy.interpolate import interp1d


    plt.style.use('gm2.mplstyle')
    plt.ion()
    
    
    #these functions are copies of the main MC ones to allow this module to be run solo
    
    def angleSpectrum(E,theta,Pu=1):
        muonM = 105.66
    
        Emax = muonM/2
        
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
        
    def acc_filter(data,function,xs):
        
        scan = [function(i) for i in xs]    
        maxfn = max(scan)   

        
        samples = []
        
        for d in data:
            u = np.random.uniform(0, maxfn)
            if u <= momAcceptanceFunction(d.P):
                samples.append(d)
        return samples
    
    def mom_acc_cut(data):
        output = acc_filter(data,momAcceptanceFunction,np.linspace(0.511,3094,50))
        
        return output
    
    def momAcceptanceFunction(x):
        #applies the volumes hit momentum cut
        xs = np.linspace(0.511,3094,10)
        y = [0,0,14,23,26,25,23,18,10,0]
        
        fit = interp1d(xs, y,kind='linear')
        
        return fit(x)  
  

    
    class decay:
        
        def __init__(self):
            radius = 7.1
            self.ringAngle = np.random.uniform(-np.pi,np.pi)
            self.x = np.cos(self.ringAngle)*radius
            self.y = np.sin(self.ringAngle)*radius
            self.P = np.random.uniform(0.511,3094)
            self.E = np.sqrt(self.P**2+0.51**2)
            self.arclength = radius*self.ringAngle
            self.hitTracker = False
            self.nhits = 0
            self.angle = sample_angle(1,self.P)
            

            
            self.xhit = 0
            self.yhit = 0
            
            self.hits = []
            
            
        def boost(self):        
        #method to boost the muon into the detector frame 
            gamma = 29.3 
            muonM = 105.66        
            mu_momentum = 3.094e3 #3.094 GeV/c
            mu_v = -mu_momentum/(gamma*muonM) #-tive because boosting back into lab frame
            
            Px = self.P*np.sin(self.angle)
            Pz = self.P*np.cos(self.angle)
            
            boostedPz = gamma*(Pz-mu_v*self.E)
            #boostedE = gamma*(self.E-mu_v*Pz) 
            
            #self.P = np.sqrt(boostedPz**2+Px**2+Py**2)
            #self.E = boostedE
            self.angle = np.arctan(Px/boostedPz)



    
    decays = []
    
    for i in range(100000):
        d = decay()
        d.boost()
        decays.append(d)
        
       
    x = []
    y = []
    p = []
    arcs = []
    for i in decays:
        x.append(i.x)
        y.append(i.y)
        
    
    fig, ax = plt.subplots(figsize=(6, 6))
    
    #define tracker position and boundaries
    x_max = 6.9
    x_min = 6.7
    y_min = -1
    y_max = 0
    
    nplanes = 32
    plane_cut = 5
    
    cutdist = plane_cut/nplanes
    
    plane = 1/32

    ring = plt.Circle((0,0),radius=7.1,alpha=0.5,fill=False)
    hits, circx,circy, cx, cy= hit_tracker(decays,x_min,x_max,y_min,y_max)
    
    
    nocuthits = hits #all hits
    
    hits = mom_acc_cut(hits) #hits with acceptance cuts
    
    pos_track = []
    
    for i,j in zip(circx,circy):
        pos_track.append((i,j))
        
    arcs = []
    col = 'r'
    if len(hits) > 0:
        col = 'green'
        
    nohits = []
        
    xhits = []
    yhits = []  
    
    momp = []
    ang = []
    vert = [] 
    
    nocutp = []
    
    for i in nocuthits:
        nocutp.append(i.P)

        
    
    for i in hits:
        
        x.append(i.x)
        y.append(i.y)
        xhits.append(i.xhit)
        yhits.append(i.yhit)
        p.append(i.P)

        
        ang.append(i.angle)
        
        
        #perform the nHits cut - compare yin/yout to see how many planes hit
        y_first = i.hits[0][1]
        y_last = i.hits[-1][1]
        
        xhit = i.hits[0][0]
        
        ydiff = abs(y_last-y_first)

        
        if ydiff >= 11*plane and i.P > 400:
            nohits.append(np.floor(ydiff/plane))
            momp.append(i.P)
            arcs.append(i.arclength)
        elif ydiff >= 5*plane and ydiff < 11*plane and i.P > 400:
            
            ran = np.random.uniform(5,11)

            if np.floor(ydiff/plane) >= ran:
                
                nohits.append(np.floor(ydiff/plane))
                momp.append(i.P)
                arcs.append(i.arclength)
                
                
    '''
    plt.figure(1)
    #this plots the actual vertices and tracks
    #for unknown reasons it can't handle plotting more than a few tracks at once
    
    ax.add_artist(ring)
    ax.scatter(x,y)
    ax.plot([i[0] for i in pos_track],[i[1] for i in pos_track],color=col)
    
    
    plt.scatter(xhits,yhits,marker='x')

    
    
    plt.plot(np.linspace(x_min,x_max,3),[y_max,y_max,y_max],color='k')
    plt.plot([x_max,x_max,x_max],np.linspace(y_min,y_max,3),color='k')
    plt.plot([x_min,x_min,x_min],np.linspace(y_min,y_max,3),color='k')
    plt.plot(np.linspace(x_min,x_max,3),[y_min,y_min,y_min],color='k')
    
    
    plt.xlabel("x position [m]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
    plt.ylabel("y position [m]",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0,labelpad=20)
    
    
    plt.xlim(-8,8)
    plt.ylim(-8,8)
    '''
    plt.figure(2)
    bins = np.linspace(0,3090,30)    
  
    
    plt.hist(p,bins,histtype='step', label='All hits')
    plt.hist(momp,bins,histtype='step', label = 'Trackable hits')
    
    pre = np.histogram(nocutp,bins)[0]
    post = np.histogram(momp,bins)[0]   
    ratio = post/pre

    plt.xlabel("Momentum [MeV]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
    plt.ylabel("Counts/bin",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0,labelpad=20)
    plt.legend()
    

    plt.figure(3)
    plt.hist(arcs,bins=50,histtype='step')
    #plt.hist2d(p,ang,bins=(20,np.linspace(-0.01,0.01,20)))
    plt.xlabel("Arc length [m]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
    plt.ylabel("Counts/bin",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0,labelpad=20)
    
    plt.figure(4)
    binz = np.linspace(0,3100,10)
    binz = np.insert(binz,1,300)
    corr = np.array([ 0.        ,  0, 0.14418263,  0.6452906 ,  0.92512741,  0.99520434,
    0.97328187,  0.880436  ,  0.72599633,  0.62138323,  0.        ])
    plt.scatter(binz,corr/1.5,label='Gas gun correction')
    plt.scatter(bins[:-1],ratio, label='nHits + high P correction')
    plt.xlabel("Momentum [MeV]",horizontalalignment='right', x=1.0, verticalalignment='bottom', y=0.0)
    plt.ylabel("Efficiency",horizontalalignment='right', y=1.0, verticalalignment='bottom', x=0.0,labelpad=20)
    
    plt.legend()

    


    

