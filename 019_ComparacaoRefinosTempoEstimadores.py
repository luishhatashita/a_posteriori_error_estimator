#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 14:18:25 2022

@author: luishatashita
"""

import matplotlib.pyplot as plt

uniform_1 = [0.0220, 0.0163, 0.0883, 0.0126, 0.0169]
adaptive_1 = [0.0129, 0.0189, 0.0332, 0.0136, 0.0780]

uniform_1_2 = [0.0091, 0.0044, 0.0068, 0.0016, 0.0019]
adaptive_1_2 = [0.0017, 0.0018, 0.0074, 0.0104, 0.0114]

uniform_no_outlier_1 = [0.0220, 0.0163, 0.0126, 0.0169]
adaptive_no_outlier_1 = [0.0129, 0.0189, 0.0332, 0.0136]

data_1 = [uniform_1, adaptive_1, uniform_1_2, adaptive_1_2]

data_1_2 = [uniform_1_2, adaptive_1_2]

data_no_outlier_1 = [uniform_no_outlier_1, adaptive_no_outlier_1, uniform_1_2, adaptive_1_2]

plt.figure(1)
plt.subplot(121)
plt.boxplot(data_no_outlier_1, showmeans=True)
plt.xlabel('Refinement/Estimator Type')
plt.xticks([1, 2, 3, 4], ['Uniform 1', 'Adaptive 1', 'Uniform 2', 'Adaptive 2'])
plt.ylabel('Time [$s$]')
plt.title('(a) Estimator Comparison')
plt.subplot(122)
plt.boxplot(data_1_2, showmeans=True)
plt.xlabel('Refinement Type')
plt.xticks([1, 2], ['Uniform 2', 'Adaptive 2'])
plt.ylabel('Time [$s$]')
plt.title('(b) Second Estimator')

plt.suptitle('Estimator Comparison - Elapsed Time')

uniform_2_1 = [0.0914, 0.1687, 0.1141, 0.1160, 0.1165]
adaptive_2_1 = [0.6081, 0.4992, 0.5638, 0.6721, 0.5533]

uniform_2_2 = [0.0036, 0.0197, 0.0034, 0.0055, 0.0034]
adaptive_2_2 = [0.0141, 0.0092, 0.0341, 0.0096, 0.0157]

uniform_no_outlier_2_1 = [0.0220, 0.0163, 0.0126, 0.0169]
adaptive_no_outlier_2_1 = [0.0129, 0.0189, 0.0332, 0.0136]

data_2 = [uniform_2_1, adaptive_2_1, uniform_2_2, adaptive_2_2]

data_2_2 = [uniform_2_2, adaptive_2_2]

data_no_outlier_2_1 = [uniform_no_outlier_2_1, adaptive_no_outlier_2_1, uniform_1_2, adaptive_1_2]

plt.figure(2)
plt.subplot(121)
plt.boxplot(data_2, showmeans=True)
plt.xlabel('Refinement/Estimator Type')
plt.xticks([1, 2, 3, 4], ['Uniform 1', 'Adaptive 1', 'Uniform 2', 'Adaptive 2'])
plt.ylabel('Time [$s$]')
plt.title('(a) Estimator Comparison')
plt.subplot(122)
plt.boxplot(data_1_2, showmeans=True)
plt.xlabel('Refinement Type')
plt.xticks([1, 2], ['Uniform 2', 'Adaptive 2'])
plt.ylabel('Time [$s$]')
plt.title('(b) Second Estimator')

plt.suptitle('Estimator Comparison - Elapsed Time')
 
# uniform_2 = [0.0283, 0.0185, 0.0201, 0.0315, 0.0297]
# adaptive_2 = [0.0457, 0.0361, 0.0530, 0.0769, 0.0485]

# adaptive_no_outlier_2 = [0.0457, 0.0361, 0.0530, 0.0485]

# data_2 = [uniform_2, adaptive_2]

# data_no_outlier_2 = [uniform_2, adaptive_no_outlier_2]

# uniform_3 = [0.2233, 0.2125, 0.1895, 0.1861, 0.1845]
# adaptive_3 = [0.4098, 0.3688, 0.4356, 0.4750, 0.4175]

# adaptive_no_outlier_3 = [0.4098, 0.4356, 0.4175]

# data_3 = [uniform_3, adaptive_3]

# data_no_outlier_3 = [uniform_3, adaptive_no_outlier_3]

# plt.figure(2)
# plt.subplot(221)
# plt.boxplot(data_2, showmeans=True)
# # plt.xlabel('Refinement Type')
# plt.xticks([1, 2], ['Uniform', 'Adaptive'])
# plt.ylabel('Time [$s$]')
# plt.title('(2) (a) 5 samples')
# plt.subplot(222)
# plt.boxplot(data_no_outlier_2, showmeans=True)
# # plt.xlabel('Refinement Type')
# plt.xticks([1, 2], ['Uniform', 'Adaptive'])
# # plt.ylabel('Time [$s$]')
# plt.title('(2) (b) Without Outlier')
# plt.subplot(223)
# plt.boxplot(data_3, showmeans=True)
# plt.xlabel('Refinement Type')
# plt.xticks([1, 2], ['Uniform', 'Adaptive'])
# plt.ylabel('Time [$s$]')
# plt.title('(3) (c) 5 samples')
# plt.subplot(224)
# plt.boxplot(data_no_outlier_3, showmeans=True)
# plt.xlabel('Refinement Type')
# plt.xticks([1, 2], ['Uniform', 'Adaptive'])
# # plt.ylabel('Time [$s$]')
# plt.title('(3) (d) Without Outliers')

# plt.suptitle('Mesh Comparison - Elapsed Time')

plt.show()
