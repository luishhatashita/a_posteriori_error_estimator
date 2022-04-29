#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 31 10:07:14 2021

@author: luishatashita
"""

# support libraries

import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt

# initial hypothesis

Ne = 10

N = Ne + 1
Lx = 1
dx = Lx/(N-1)

# domain

x = np.linspace(0, Lx, N) # split the domain into Ne elements

# matrices for the linear system Au = b

A = np.zeros((N, N))
b = np.zeros((N, 1))

# populate the rigidity matrix and the charge vector
# python utilises zero indexing

for i in range(0, N-1):
    Be_i = x[i+1] - x[i] # jacobian
    A[i, i] = A[i,i] + np.power(Be_i,-1)
    A[i, i+1] = A[i, i+1] - np.power(Be_i,-1)
    A[i+1, i] = A[i+1, i] - np.power(Be_i,-1)
    A[i+1, i+1] = A[i+1, i+1] + np.power(Be_i,-1)
    b[i, 0] = b[i, 0] - Be_i/2
    b[i+1, 0] = b[i+1, 0] - Be_i/2

# boundary values
# assure the boundary values readequating the rigidity matrix first

A[0, 0] = 1
A[0, 1] = 0
A[N-1, N-1] = 1
A[N-1, N-2] = 0
b[0, 0] = 1 # u(0) = 1
b[N-1, 0] = 1 # u(1) = 1

# solving the system, through u = A^-1.b

u = linalg.inv(A).dot(b)

# plot

plt.plot(x, u)
plt.xlabel("$x$")
plt.ylabel("$u(x)$")
plt.title("$-u^{''}(x) = -1$ with $u(0)=u(1)=1$")
plt.grid(True)
