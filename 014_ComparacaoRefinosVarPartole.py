#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 16 16:18:41 2022

@author: luishatashita
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

tol_e = [0.01, 0.005, 0.001, 0.0005, 0.0001]

# Ne_1 = 5
sum_etaK_Ne_uni_1 = [
    0.1503,
    0.1503,
    0.0781,
    0.0781,
    0.0394
]

sum_etaK_Ne_adpt_1 = [
    0.1072,
    0.1002,
    0.0486,
    0.0478,
    0.0239
]

final_ne_Ne_uni_1 = [40, 40, 160, 160, 640]

final_ne_Ne_adpt_1 = [28, 32, 110, 114, 406]

mean_time_Ne_uni_1 = [
    0.0164,
    0.0087,
    0.0896,
    0.0281,
    0.1171
]

mean_time_Ne_adpt_1 = [
    0.0313,
    0.0347,
    0.3150,
    0.4661,
    0.8750
]

final_it_Ne_uni_1 = [3, 3, 5, 5, 7]

final_it_Ne_adpt_1 = [4, 6, 11, 13, 17]

# Ne_2 = 25
sum_etaK_Ne_uni_2 = [
    0.1843,
    0.1358,
    0.0700,
    0.0700,
    0.0353
]

sum_etaK_Ne_adpt_2 = [
    0.1843,
    0.1309,
    0.0665,
    0.0657,
    0.0327
]

final_ne_Ne_uni_2 = [25, 50, 200, 200, 800]

final_ne_Ne_adpt_2 = [25, 48, 188, 192, 750]

mean_time_Ne_uni_2 = [
    0.0022,
    0.0153,
    0.0364,
    0.0430,
    0.2247
]

mean_time_Ne_adpt_2 = [
    0.0034,
    0.0079,
    0.0520,
    0.0880,
    0.9468
]

final_it_Ne_uni_2 = [0, 1, 3, 3, 5]

final_it_Ne_adpt_2 = [0, 1, 4, 6, 11]

# Ne_3 = 100
sum_etaK_Ne_uni_3 = [
    0.0980,
    0.0980,
    0.0700,
    0.0700,
    0.0353
]

sum_etaK_Ne_adpt_3 = [
    0.0980,
    0.0980,
    0.0694,
    0.0694,
    0.0348
]

final_ne_Ne_uni_3 = [100, 100, 200, 200, 800]

final_ne_Ne_adpt_3 = [100, 100, 198, 198, 790]

mean_time_Ne_uni_3 = [
    0.0138,
    0.0066,
    0.0590,
    0.0335,
    0.2362
]

mean_time_Ne_adpt_3 = [
    0.0073,
    0.0052,
    0.0248,
    0.0318,
    0.4214
]

final_it_Ne_uni_3 = [0, 0, 1, 1, 3]

final_it_Ne_adpt_3 = [0, 0, 1, 1, 5]

plt.figure(1)
plt.subplot(221)
plt.plot(tol_e, sum_etaK_Ne_uni_1, color="#ff7f0e", marker="s", ls='--', label='(uni) Ne_0 = 5')
plt.plot(tol_e, sum_etaK_Ne_uni_2, color="#ff7f0e", marker="^", ls='--', label='(uni) Ne_0 = 25')
plt.plot(tol_e, sum_etaK_Ne_uni_3, color="#ff7f0e", marker="D", ls='--', label='(uni) Ne_0 = 100')
plt.plot(tol_e, sum_etaK_Ne_adpt_1, color="#1f77b4", marker="s", ls='--', label='(adpt) Ne_0 = 5')
plt.plot(tol_e, sum_etaK_Ne_adpt_2, color="#1f77b4", marker="^", ls='--', label='(adpt) Ne_0 = 25')
plt.plot(tol_e, sum_etaK_Ne_adpt_3, color="#1f77b4", marker="D", ls='--', label='(adpt) Ne_0 = 100')
# plt.xlabel("Individual Error Tolerance")
plt.xscale('log')
# plt.xlim([0.01025, -0.00025])
plt.gca().invert_xaxis()
plt.ylabel("$\Sigma\eta_K$")
plt.title("(a) Total error")
plt.legend(loc="best")
plt.grid(True)

# plt.figure(2)
plt.subplot(222)
plt.plot(tol_e, final_ne_Ne_uni_1, color="#ff7f0e", marker="s", ls='--', label='(uni) Ne_0 = 5')
plt.plot(tol_e, final_ne_Ne_uni_2, color="#ff7f0e", marker="^", ls='--', label='(uni) Ne_0 = 25')
plt.plot(tol_e, final_ne_Ne_uni_3, color="#ff7f0e", marker="D", ls='--', label='(uni) Ne_0 = 100')
plt.plot(tol_e, final_ne_Ne_adpt_1, color="#1f77b4", marker="s", ls='--', label='(adpt) Ne_0 = 5')
plt.plot(tol_e, final_ne_Ne_adpt_2, color="#1f77b4", marker="^", ls='--', label='(adpt) Ne_0 = 25')
plt.plot(tol_e, final_ne_Ne_adpt_3, color="#1f77b4", marker="D", ls='--', label='(adpt) Ne_0 = 100')
# plt.xlabel("Individual Error Tolerance")
plt.xscale('log')
plt.gca().invert_xaxis()
plt.ylabel("Final Number of Elements")
plt.title("(b) Final number of elements")
# plt.legend(loc="best")
plt.grid(True)

# plt.figure(3)
plt.subplot(223)
plt.plot(tol_e, mean_time_Ne_uni_1, color="#ff7f0e", marker="s", ls='--', label='(uni) Ne_0 = 5')
plt.plot(tol_e, mean_time_Ne_uni_2, color="#ff7f0e", marker="^", ls='--', label='(uni) Ne_0 = 25')
plt.plot(tol_e, mean_time_Ne_uni_3, color="#ff7f0e", marker="D", ls='--', label='(uni) Ne_0 = 100')
plt.plot(tol_e, mean_time_Ne_adpt_1, color="#1f77b4", marker="s", ls='--', label='(adpt) Ne_0 = 5')
plt.plot(tol_e, mean_time_Ne_adpt_2, color="#1f77b4", marker="^", ls='--', label='(adpt) Ne_0 = 25')
plt.plot(tol_e, mean_time_Ne_adpt_3, color="#1f77b4", marker="D", ls='--', label='(adpt) Ne_0 = 100')
plt.xlabel("Individual Error Tolerance")
plt.xscale('log')
plt.gca().invert_xaxis()
plt.ylabel("Time [$s$]")
plt.title("(c) Elapsed time")
# plt.legend(loc="best")
plt.grid(True)

# plt.figure(4)
plt.subplot(224)
plt.plot(tol_e, final_it_Ne_uni_1, color="#ff7f0e", marker="s", ls='--', label='(uni) Ne_0 = 5')
plt.plot(tol_e, final_it_Ne_uni_2, color="#ff7f0e", marker="^", ls='--', label='(uni) Ne_0 = 25')
plt.plot(tol_e, final_it_Ne_uni_3, color="#ff7f0e", marker="D", ls='--', label='(uni) Ne_0 = 100')
plt.plot(tol_e, final_it_Ne_adpt_1, color="#1f77b4", marker="s", ls='--', label='(adpt) Ne_0 = 5')
plt.plot(tol_e, final_it_Ne_adpt_2, color="#1f77b4", marker="^", ls='--', label='(adpt) Ne_0 = 25')
plt.plot(tol_e, final_it_Ne_adpt_3, color="#1f77b4", marker="D", ls='--', label='(adpt) Ne_0 = 100')
plt.xlabel("Individual Error Tolerance")
plt.xscale('log')
plt.gca().invert_xaxis()
plt.ylabel("Iterations")
plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
plt.title("(d) Number of iterations")
# plt.legend(loc="best")
plt.grid(True)

plt.suptitle("Influence of individual error tolerance")

# plt.figure(2)
# plt.subplot(121)
# plt.plot(tol_e, sum_etaK_tol_e_uni_1, color="#ff7f0e", marker="s", ls='--', label='tol_e = 0.01')
# plt.plot(tol_e, sum_etaK_tol_e_uni_2, color="#ff7f0e", marker="^", ls='--', label='tol_e = 0.001')
# plt.plot(tol_e, sum_etaK_tol_e_uni_3, color="#ff7f0e", marker="D", ls='--', label='tol_e = 0.0001')
# plt.xlabel("Number of Elements")
# plt.ylabel("$\Sigma\eta_K^2$")
# plt.title("Uniform Mesh")
# plt.legend(loc="best")
# plt.grid(True)
# plt.subplot(122)
# plt.plot(tol_e, sum_etaK_tol_e_adpt_1, color="#1f77b4", marker="s", ls='--', label='tol_e = 0.01')
# plt.plot(tol_e, sum_etaK_tol_e_adpt_2, color="#1f77b4", marker="^", ls='--', label='tol_e = 0.001')
# plt.plot(tol_e, sum_etaK_tol_e_adpt_3, color="#1f77b4", marker="D", ls='--', label='tol_e = 0.0001')
# plt.xlabel("Number of Elements")
# plt.ylabel("$\Sigma\eta_K^2$")
# plt.title("Adaptive Mesh")
# plt.legend(loc="best")
# plt.grid(True)

# plt.suptitle("Influence of initial number of elements")

plt.show()
