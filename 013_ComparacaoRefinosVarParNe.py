#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 16 16:18:41 2022

@author: luishatashita
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

initial_ne = [5, 10, 25, 50, 100]

# tol_e_1 = 0.01
sum_etaK_tol_e_uni_1 = [
    0.1503,
    0.1503,
    0.1843,
    0.1358,
    0.0980
]

sum_etaK_tol_e_adpt_1 = [
    0.1072,
    0.1323,
    0.1843,
    0.1358,
    0.0980
]

final_ne_tol_e_uni_1 = [40, 40, 25, 50, 100]

final_ne_tol_e_adpt_1 = [28, 34, 25, 50, 100]

mean_time_tol_e_uni_1 = [
    0.0164,
    0.0103,
    0.0022,
    0.0109,
    0.0138
]

mean_time_tol_e_adpt_1 = [
    0.0313,
    0.0115,
    0.0034,
    0.0117,
    0.0073
]

final_it_tol_e_uni_1 = [3, 2, 0, 0, 0]

final_it_tol_e_adpt_1 = [4, 2, 0, 0, 0]

# tol_e_2 = 0.001
sum_etaK_tol_e_uni_2 = [
    0.0781,
    0.0781,
    0.0700,
    0.0700,
    0.0700
]

sum_etaK_tol_e_adpt_2 = [
    0.0486,
    0.0646,
    0.0665,
    0.0684,
    0.0694
]

final_ne_tol_e_uni_2 = [160, 160, 200, 200, 200]

final_ne_tol_e_adpt_2 = [110, 138, 188, 194, 198]

mean_time_tol_e_uni_2 = [
    0.0896,
    0.0390,
    0.0364,
    0.0319,
    0.0590
]

mean_time_tol_e_adpt_2 = [
    0.3150,
    0.0928,
    0.0520,
    0.0482,
    0.0248
]

final_it_tol_e_uni_2 = [5, 4, 3, 2, 1]

final_it_tol_e_adpt_2 = [11, 8, 4, 2, 1]

# tol_e_3 = 0.0001
sum_etaK_tol_e_uni_3 = [
    0.0394,
    0.0394,
    0.0353,
    0.0353,
    0.0353
]

sum_etaK_tol_e_adpt_3 = [
    0.0239,
    0.0317,
    0.0327,
    0.0341,
    0.0348
]

final_ne_tol_e_uni_3 = [640, 640, 800, 800, 800]

final_ne_tol_e_adpt_3 = [406, 532, 750, 778, 790]

mean_time_tol_e_uni_3 = [
    0.1171,
    0.1213,
    0.2247,
    0.3044,
    0.2362
]

mean_time_tol_e_adpt_3 = [
    0.8750,
    0.5793,
    0.9468,
    0.5705,
    0.4214
]

final_it_tol_e_uni_3 = [7, 6, 5, 4, 3]

final_it_tol_e_adpt_3 = [17, 15, 11, 8, 5]

plt.figure(1)
plt.subplot(221)
plt.plot(initial_ne, sum_etaK_tol_e_uni_1, color="#ff7f0e", marker="s", ls='--', label='(uni) tol_e = 0.01')
plt.plot(initial_ne, sum_etaK_tol_e_uni_2, color="#ff7f0e", marker="^", ls='--', label='(uni) tol_e = 0.001')
plt.plot(initial_ne, sum_etaK_tol_e_uni_3, color="#ff7f0e", marker="D", ls='--', label='(uni) tol_e = 0.0001')
plt.plot(initial_ne, sum_etaK_tol_e_adpt_1, color="#1f77b4", marker="s", ls='--', label='(adpt) tol_e = 0.01')
plt.plot(initial_ne, sum_etaK_tol_e_adpt_2, color="#1f77b4", marker="^", ls='--', label='(adpt) tol_e = 0.001')
plt.plot(initial_ne, sum_etaK_tol_e_adpt_3, color="#1f77b4", marker="D", ls='--', label='(adpt) tol_e = 0.0001')
# plt.xlabel("Initial Number of Elements")
plt.ylabel("$\Sigma\eta_K$")
plt.title("(a) Total error")
# plt.legend(loc="best")
plt.grid(True)

plt.subplot(222)
plt.plot(initial_ne, final_ne_tol_e_uni_1, color="#ff7f0e", marker="s", ls='--', label='(uni) tol_e = 0.01')
plt.plot(initial_ne, final_ne_tol_e_uni_2, color="#ff7f0e", marker="^", ls='--', label='(uni) tol_e = 0.001')
plt.plot(initial_ne, final_ne_tol_e_uni_3, color="#ff7f0e", marker="D", ls='--', label='(uni) tol_e = 0.0001')
plt.plot(initial_ne, final_ne_tol_e_adpt_1, color="#1f77b4", marker="s", ls='--', label='(adpt) tol_e = 0.01')
plt.plot(initial_ne, final_ne_tol_e_adpt_2, color="#1f77b4", marker="^", ls='--', label='(adpt) tol_e = 0.001')
plt.plot(initial_ne, final_ne_tol_e_adpt_3, color="#1f77b4", marker="D", ls='--', label='(adpt) tol_e = 0.0001')
# plt.xlabel("Initial Number of Elements")
plt.ylabel("Final Number of Elements")
plt.title("(b) Final number of elements")
# plt.legend(loc="best")
plt.grid(True)

plt.subplot(223)
plt.plot(initial_ne, mean_time_tol_e_uni_1, color="#ff7f0e", marker="s", ls='--', label='(uni) tol_e = 0.01')
plt.plot(initial_ne, mean_time_tol_e_uni_2, color="#ff7f0e", marker="^", ls='--', label='(uni) tol_e = 0.001')
plt.plot(initial_ne, mean_time_tol_e_uni_3, color="#ff7f0e", marker="D", ls='--', label='(uni) tol_e = 0.0001')
plt.plot(initial_ne, mean_time_tol_e_adpt_1, color="#1f77b4", marker="s", ls='--', label='(adpt) tol_e = 0.01')
plt.plot(initial_ne, mean_time_tol_e_adpt_2, color="#1f77b4", marker="^", ls='--', label='(adpt) tol_e = 0.001')
plt.plot(initial_ne, mean_time_tol_e_adpt_3, color="#1f77b4", marker="D", ls='--', label='(adpt) tol_e = 0.0001')
plt.xlabel("Initial Number of Elements")
plt.ylabel("Time [$s$]")
plt.title("(c) Elapsed time")
# plt.legend(loc="best")
plt.grid(True)

plt.subplot(224)
plt.plot(initial_ne, final_it_tol_e_uni_1, color="#ff7f0e", marker="s", ls='--', label='(uni) tol_e = 0.01')
plt.plot(initial_ne, final_it_tol_e_uni_2, color="#ff7f0e", marker="^", ls='--', label='(uni) tol_e = 0.001')
plt.plot(initial_ne, final_it_tol_e_uni_3, color="#ff7f0e", marker="D", ls='--', label='(uni) tol_e = 0.0001')
plt.plot(initial_ne, final_it_tol_e_adpt_1, color="#1f77b4", marker="s", ls='--', label='(adpt) tol_e = 0.01')
plt.plot(initial_ne, final_it_tol_e_adpt_2, color="#1f77b4", marker="^", ls='--', label='(adpt) tol_e = 0.001')
plt.plot(initial_ne, final_it_tol_e_adpt_3, color="#1f77b4", marker="D", ls='--', label='(adpt) tol_e = 0.0001')
plt.xlabel("Initial Number of Elements")
plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
plt.ylabel("Iterations")
plt.title("(d) Number of iterations")
plt.legend(loc="best")
plt.grid(True)

plt.suptitle("Influence of initial number of elements")

plt.figure(2)
plt.subplot(121)
plt.plot(initial_ne, sum_etaK_tol_e_uni_1, color="#ff7f0e", marker="s", ls='--', label='tol_e = 0.01')
plt.plot(initial_ne, sum_etaK_tol_e_uni_2, color="#ff7f0e", marker="^", ls='--', label='tol_e = 0.001')
plt.plot(initial_ne, sum_etaK_tol_e_uni_3, color="#ff7f0e", marker="D", ls='--', label='tol_e = 0.0001')
plt.xlabel("Number of Elements")
plt.ylabel("$\Sigma\eta_K^2$")
plt.title("Uniform Mesh")
plt.legend(loc="best")
plt.grid(True)
plt.subplot(122)
plt.plot(initial_ne, sum_etaK_tol_e_adpt_1, color="#1f77b4", marker="s", ls='--', label='tol_e = 0.01')
plt.plot(initial_ne, sum_etaK_tol_e_adpt_2, color="#1f77b4", marker="^", ls='--', label='tol_e = 0.001')
plt.plot(initial_ne, sum_etaK_tol_e_adpt_3, color="#1f77b4", marker="D", ls='--', label='tol_e = 0.0001')
plt.xlabel("Number of Elements")
plt.ylabel("$\Sigma\eta_K^2$")
plt.title("Adaptive Mesh")
plt.legend(loc="best")
plt.grid(True)

plt.suptitle("Influence of initial number of elements")

plt.show()
