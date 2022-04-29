#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 16:35:03 2022

@author: luishatashita
"""

import matplotlib.pyplot as plt

x_uni = [0.   , 0.025, 0.05 , 0.075, 0.1  , 0.125, 0.15 , 0.175, 0.2  ,
       0.225, 0.25 , 0.275, 0.3  , 0.325, 0.35 , 0.375, 0.4  , 0.425,
       0.45 , 0.475, 0.5  , 0.525, 0.55 , 0.575, 0.6  , 0.625, 0.65 ,
       0.675, 0.7  , 0.725, 0.75 , 0.775, 0.8  , 0.825, 0.85 , 0.875,
       0.9  , 0.925, 0.95 , 0.975, 1.   ]

# x_uni.pop(0)
# x_uni.pop()

x_adpt = [0.    , 0.2   , 0.2125, 0.225 , 0.25  , 0.275 , 0.3   , 0.325 ,
       0.35  , 0.375 , 0.4   , 0.425 , 0.45  , 0.475 , 0.5   , 0.525 ,
       0.55  , 0.575 , 0.6   , 0.625 , 0.65  , 0.675 , 0.7   , 0.725 ,
       0.75  , 0.775 , 0.7875, 0.8   , 1.    ]

# x_adpt.pop(0)
# x_adpt.pop()

n_bins = 41

plt.figure(1)
plt.subplot(221)
plt.hist(x_uni, bins=n_bins, facecolor='#ff7f0e', edgecolor='#e0e0e0', linewidth=0.5)
plt.axis([0, 1, 0, 2.5])
plt.xlabel('$x$')
plt.ylabel('Number of Elements')
plt.title('(a) Uniform Refinement')
plt.subplot(222)
plt.hist(x_adpt, bins=n_bins, facecolor='#1f77b4', edgecolor='#e0e0e0', linewidth=0.5)
plt.axis([0, 1, 0, 2.5])
plt.xlabel('$x$')
plt.ylabel('Number of Elements')
plt.title('(b) Adaptive Refinement')
plt.subplot(223)
plt.hist(x_uni, bins=n_bins, facecolor='#ff7f0e', edgecolor='#e0e0e0', linewidth=0.5)
plt.axis([0.1, 0.3, 0, 2.5])
plt.xlabel('$x$')
plt.ylabel('Number of Elements')
plt.title('(c) Estimator border zoom')
plt.subplot(224)
plt.hist(x_adpt, bins=n_bins, facecolor='#1f77b4', edgecolor='#e0e0e0', linewidth=0.5)
plt.axis([0.1, 0.3, 0, 2.5])
plt.xlabel('$x$')
plt.ylabel('Number of Elements')
plt.title('(d) Estimator border zoom')

plt.suptitle('Final Domain Histogram')

plt.show()
