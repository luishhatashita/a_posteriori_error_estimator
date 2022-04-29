#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 15:20:53 2022

@author: luishatashita
"""

import numpy as np
# from numpy import *
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set()

x_ne = np.array([
    5,
    5,
    5,
    5,
    5,
    10,
    10,
    10,
    10,
    10,
    25,
    25,
    25,
    25,
    25,
    50,
    50,
    50,
    50,
    50,
    100,
    100,
    100,
    100,
    100    
])

y_tol_e = np.array([
    0.01,
    0.005,
    0.001,
    0.0005,
    0.0001,
    0.01,
    0.005,
    0.001,
    0.0005,
    0.0001,
    0.01,
    0.005,
    0.001,
    0.0005,
    0.0001,
    0.01,
    0.005,
    0.001,
    0.0005,
    0.0001,
    0.01,
    0.005,
    0.001,
    0.0005,
    0.0001
])

z_etaK_uni = np.array([
    0.1503,
    0.1503,
    0.0781,
    0.0781,
    0.0394,
    0.1503,
    0.1503,
    0.0781,
    0.0781,
    0.0394,
    0.1843,
    0.1358,
    0.0700,
    0.0700,
    0.0353,
    0.1358,
    0.1358,
    0.0700,
    0.0700,
    0.0353,
    0.0980,
    0.0980,
    0.0700,
    0.0700,
    0.0353    
])

z_etaK_adpt = np.array([
    0.1072,
    0.1002,
    0.0486,
    0.0478,
    0.0239,
    0.1323,
    0.1311,
    0.0646,
    0.0600,
    0.0317,
    0.1843,
    0.1309,
    0.0665,
    0.0657,
    0.0327,
    0.1358,
    0.1358,
    0.0684,
    0.0684,
    0.0341,
    0.0980,
    0.0980,
    0.0694,
    0.0694,
    0.0348
])

z_time_uni = np.array([
   0.0164,
   0.0087,
   0.0896,
   0.0281,
   0.1171,
   0.0103,
   0.0155,
   0.0390,
   0.0281,
   0.1213,
   0.0022,
   0.0153,
   0.0364,
   0.0430,
   0.2247,
   0.0109,
   0.0039,
   0.0319,
   0.0294,
   0.3044,
   0.0138,
   0.0066,
   0.0590,
   0.0335,
   0.2362 
])

z_time_adpt = np.array([
   0.0313,
   0.0347,
   0.3150,
   0.4661,
   0.8750,
   0.0115,
   0.0069,
   0.0928,
   0.1989,
   0.5793,
   0.0034,
   0.0079,
   0.0520,
   0.0880,
   0.9468,
   0.0117,
   0.0085,
   0.0482,
   0.0576,
   0.5705,
   0.0073,
   0.0052,
   0.0248,
   0.0318,
   0.4214
])

df_etaK_uni = pd.DataFrame.from_dict(np.array([x_ne,y_tol_e,z_etaK_uni]).T)
df_etaK_uni.columns = ['Initial Number of Elements','Individual Error Tolerance','Total Error']

pivotted_etaK_uni = df_etaK_uni.pivot('Individual Error Tolerance','Initial Number of Elements','Total Error')

df_etaK_adpt = pd.DataFrame.from_dict(np.array([x_ne,y_tol_e,z_etaK_adpt]).T)
df_etaK_adpt.columns = ['Initial Number of Elements','Individual Error Tolerance','Total Error']

pivotted_etaK_adpt = df_etaK_adpt.pivot('Individual Error Tolerance','Initial Number of Elements','Total Error')

df_time_uni = pd.DataFrame.from_dict(np.array([x_ne,y_tol_e,z_time_uni]).T)
df_time_uni.columns = ['Initial Number of Elements','Individual Error Tolerance','Elapsed Time']

pivotted_time_uni = df_time_uni.pivot('Individual Error Tolerance','Initial Number of Elements','Elapsed Time')

df_time_adpt = pd.DataFrame.from_dict(np.array([x_ne,y_tol_e,z_time_adpt]).T)
df_time_adpt.columns = ['Initial Number of Elements','Individual Error Tolerance','Elapsed Time']

pivotted_time_adpt = df_time_adpt.pivot('Individual Error Tolerance','Initial Number of Elements','Elapsed Time')

# plt.figure(1)
# plt.subplot(121)
# sns.heatmap(pivotted_etaK_uni, cmap='YlGnBu', linewidths=.5)
# plt.title("(a) Uniform Mesh")
# plt.subplot(122)
# sns.heatmap(pivotted_etaK_adpt, cmap='YlGnBu', linewidths=.5)
# plt.title("(b) Adaptive Mesh")

# plt.suptitle("Total Error - $\Sigma \eta_K$")

fig1, ax1 = plt.subplots(1,2)
sns.heatmap(pivotted_etaK_uni, cmap='YlGnBu', linewidths=.5, ax=ax1[0], vmin=min(z_etaK_adpt), vmax=max(z_etaK_adpt))
ax1[0].set_title("(a) Uniform Mesh")
sns.heatmap(pivotted_etaK_adpt, cmap='YlGnBu', linewidths=.5, ax=ax1[1])
ax1[1].set_title("(b) Adaptive Mesh")

fig1.suptitle("Total Error - $\Sigma \eta_K$")

fig2, ax2 = plt.subplots(1,2)
sns.heatmap(pivotted_time_uni, cmap='rocket_r', linewidths=.5, ax=ax2[0], vmin=0, vmax=max(z_time_adpt))
ax2[0].set_title("(a) Uniform Mesh")
sns.heatmap(pivotted_time_adpt, cmap='rocket_r', linewidths=.5, ax=ax2[1], vmin=0, vmax=max(z_time_adpt))
ax2[1].set_title("(b) Adaptive Mesh")

fig2.suptitle("Elapsed Time [s]")

plt.show()

# --- heatmap example ---
# df1 = pd.DataFrame.from_dict(np.array([x_ne,y_tol_e,z_etaK_adpt]).T)
# df1.columns = ['X_value','Y_value','Z_value']

# pivotted_1 = df1.pivot('Y_value','X_value','Z_value')
# sns.heatmap(pivotted_1,cmap='YlGnBu', linewidths=.5)