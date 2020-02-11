#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 17:51:04 2020
"""
#who needs ROOT anyway?

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

#data here
t = np.linspace(0., 300, 150)
counts = [2815, 2805, 2794, 2611, 2439, 2373, 2212, 2030, 1861, 1631, 1517, 1354, 1197, 1094, 1048, 962, 889, 970, 1029, 1175, 1291, 1425, 1577, 1850, 1953, 2115, 2312, 2406, 2525, 2629, 2774, 2782, 2804, 2738, 2703, 2507, 2430, 2217, 2097, 1928, 1754, 1602, 1375, 1294, 1164, 1040, 963, 975, 936, 1027, 1068, 1186, 1324, 1550, 1698, 1890, 2001, 2242, 2380, 2458, 2665, 2764, 2782, 2823, 2766, 2715, 2554, 2465, 2269, 2102, 1979, 1783, 1586, 1531, 1394, 1233, 1053, 1033, 995, 957, 977, 1071, 1186, 1314, 1440, 1570, 1799, 2028, 2128, 2240, 2469, 2542, 2667, 2794, 2724, 2802, 2670, 2613, 2551, 2372, 2245, 2004, 1834, 1661, 1535, 1341, 1299, 1106, 1045, 1050, 936, 939, 1039, 1048, 1291, 1383, 1532, 1776, 1935, 2037, 2273, 2395, 2545, 2706, 2674, 2812, 2725, 2745, 2635, 2625, 2540, 2326, 2202, 1968, 1813, 1635, 1433, 1327, 1215, 1080, 956, 982, 1013, 995, 1073, 1186, 1261, 1474, 1676, 1816]



# Fit the wiggle with a shifted sine (no exp decay yet)
fitfunc = lambda p, x: p[0]*np.sin(p[1]*x+p[2]) + p[3] # Target function
errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
p0 = [1000, 0.1, 9.0,2000] # Initial guess for the parameters
p1, success = optimize.leastsq(errfunc, p0[:], args=(t, counts))

time = np.linspace(t.min(), t.max(), 1000)
plt.plot(t, counts, "k.", time, fitfunc(p1, time)) # Plot of the data and the fit
plt.xlabel('Time [ns]')
plt.ylabel('Number of positrons with E>2000 MeV')
plt.text(260,3000,'w from fit: %f' %p1[1])
plt.legend(('Wiggle', 'Fit'))

plt.show()