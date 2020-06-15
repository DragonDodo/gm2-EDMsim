#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 11:45:35 2020

@author: dvasilkova
"""
#python hadd for .txt files from gm2sim

import os,glob,sys
import numpy as np

folder_path = sys.argv [1]
saveto = sys.argv [2]
print("Reading from folder ", folder_path , "and  hadd-ing  to", saveto)

av_array = np.zeros(99)
er_array = np.zeros(99)
counter = 0

for filename in glob.glob(os.path.join(folder_path, '*.txt')):
    data = np.loadtxt(filename,unpack=True,delimiter=',') 
    time = data[0]
    vertical_angle = data[1]
    error = data[2]

    av_array += vertical_angle
    er_array += error**2
    counter +=1
    
av_array = av_array/counter    
er_array = np.sqrt(er_array)/counter
    
    
f = open(saveto,'w')
for i,j,k in zip(time,av_array,er_array):
    f.write(str(i)+','+str(j)+','+str(k)+'\n')
f.close()