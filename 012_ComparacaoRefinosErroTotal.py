#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 21:21:27 2022

@author: luishatashita
"""

import matplotlib.pyplot as plt

Ne_adpt_2 = [25, 48, 94, 186, 188]

sum_etaK_adpt_2 = [0.18429416485607997,
                0.13085162675731335,
                0.09355657618769571,
                0.06684245196054997,
                0.06646245692181807
                ]

Ne_uni_2 = [25, 50, 100, 200]

sum_etaK_uni_2 = [0.18429416485607997,
                0.13581879693179438,
                0.09800979951004901,
                0.07000532140487609,
                ]

Ne_adpt_3 = [100, 198, 394, 786, 788, 790]

sum_etaK_adpt_3 = [0.09800979951004901,
                0.06938244194936206,
                0.049199180561130776,
                0.034879140679898765,
                0.03483153451489044,
                0.03478032543079019
                ]

Ne_uni_3 = [100, 200, 400, 800]

sum_etaK_uni_3 = [0.09800979951004901,
                0.07000532140487609,
                0.04975031093652833,
                0.0352670058162465,
                ]

plt.figure(1)
plt.subplot(121)  
plt.plot(Ne_adpt_2, sum_etaK_adpt_2, "s", label="Adaptive mesh")
plt.plot(Ne_uni_2, sum_etaK_uni_2, "^", label="Uniform mesh")
# plt.xscale("log")
# plt.yscale("log")
plt.xlabel("Number of Elements")
plt.ylabel("$\Sigma\eta_K$")
plt.title("(2) $Ne = 25$ and $tol_e = 0.001$")
# plt.legend(loc="best")
plt.grid(True)
plt.subplot(122)  
plt.plot(Ne_adpt_3, sum_etaK_adpt_3, "s", label="Adaptive mesh")
plt.plot(Ne_uni_3, sum_etaK_uni_3, "^", label="Uniform mesh")
# plt.xscale("log")
# plt.yscale("log")
plt.xlabel("Number of Elements")
plt.ylabel("$\Sigma\eta_K$")
plt.title("(3) $Ne = 100$ and $tol_e = 0.0001$")
plt.legend(loc="best")
plt.grid(True)

plt.suptitle('Paremeters comparison - Total $\eta_K$ (sum)')

plt.figure(2)
plt.plot(Ne_adpt_2, sum_etaK_adpt_2, color="#1f77b4", marker="s",ls="--", label="(2) - Adaptive mesh")
plt.plot(Ne_uni_2, sum_etaK_uni_2, color="#ff7f0e", marker="^", ls="--", label="(2) - Uniform mesh")
plt.plot(Ne_adpt_3, sum_etaK_adpt_3, color="#b41c74", marker="*", ls="--", label="(3) - Adaptive mesh")
plt.plot(Ne_uni_3, sum_etaK_uni_3, color="#74b41c", marker="D", ls="--", label="(3) - Uniform mesh")
# plt.xscale("log")
# plt.yscale("log")
plt.xlabel("Number of Elements")
plt.ylabel("$\Sigma\eta_K$")
plt.title("Parameters - Total $\eta_K$ (Sum)")
plt.legend(loc="best")
plt.grid(True)

plt.show()