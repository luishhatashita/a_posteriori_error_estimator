#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 31 15:14:01 2021

@author: luishatashita
"""

# support libraries

import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt

# initial hypothesis

Ne = 5

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
    B_e_i = x[i+1] - x[i] # jacobian
    A[i, i] = A[i,i] + np.power(B_e_i,-1)
    A[i, i+1] = A[i, i+1] - np.power(B_e_i,-1)
    A[i+1, i] = A[i+1, i] - np.power(B_e_i,-1)
    A[i+1, i+1] = A[i+1, i+1] + np.power(B_e_i,-1)
    b[i, 0] = b[i, 0] - B_e_i/2
    b[i+1, 0] = b[i+1, 0] - B_e_i/2

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

plt.figure(1)
plt.plot(x, u, label="Iteration 0")


#-----------------------------------------------------------------------------
# error estimator
#-----------------------------------------------------------------------------

# initial hypothesis

tol_e = 0.01

# creating empty arrays for error vectors

xe = np.zeros((N-3))
r2 = np.zeros((N-3))
R2 = np.zeros((N-3))
etaK = np.zeros((N-3))

# loop through elements where the error estimator expression is valid
# from elements 1 to Ne-1

for i in range(1,N-2):
    B_e_i = x[i+1] - x[i] # jacobian, in this case the h_K
    xe[i-1] = (x[i+1] - x[i])/2 + x[i]
    r2[i-1] = B_e_i
    R2[i-1] = np.power((-u[i-1,0] + 2*u[i,0] - u[i+1,0]),2) + np.power((-u[i,0] + 2*u[i+1,0] - u[i+2,0]),2)
    etaK[i-1] = np.sqrt(B_e_i**2*r2[i-1] + B_e_i*R2[i-1]) 
    # print(i, xe[i-1], etaK[i-1])
    
plt.plot(xe, r2, label="$||r||^2$")
plt.plot(xe, R2, label="$||R||^2$")
plt.bar(xe, etaK, width=0.7*dx, label="$\eta_K$")
plt.axis([x[0], x[N-1], 0, r2.max()+0.01])
plt.xlabel("$x$")
plt.title("Estimador de Erro")
plt.legend(loc="best")
plt.grid(True)

plt.show()

#-----------------------------------------------------------------------------
# adaptive meshing
#-----------------------------------------------------------------------------

iteration = 1
x_tmp_base = x
etaK_tmp_base = etaK
N_tmp = N

while iteration <=5:
    
    # mesh refinement
    x_tmp = x_tmp_base
    for i in range(1,N_tmp-2):
        if etaK_tmp_base[i-1] > tol_e:
            index_tmp = list(x_tmp).index(x_tmp_base[i])
            # print(index_tmp)
            x_new_tmp = (x_tmp[index_tmp+1] - x_tmp[index_tmp])/2 + x_tmp[index_tmp]
            x_tmp = np.insert(x_tmp, index_tmp+1, x_new_tmp)
            # print(x_tmp)
    
    # break condition, all elements' errors are under the tolerance
    if np.size(x_tmp) == np.size(x_tmp_base):
        break
    
    # solution of the new system equation
    N_tmp = np.size(x_tmp)
    A_tmp = np.zeros((N_tmp, N_tmp))
    b_tmp = np.zeros((N_tmp, 1))
    
    for i in range(0, N_tmp-1):
        B_e_i = x_tmp[i+1] - x_tmp[i] # jacobian
        A_tmp[i, i] = A_tmp[i,i] + np.power(B_e_i,-1)
        A_tmp[i, i+1] = A_tmp[i, i+1] - np.power(B_e_i,-1)
        A_tmp[i+1, i] = A_tmp[i+1, i] - np.power(B_e_i,-1)
        A_tmp[i+1, i+1] = A_tmp[i+1, i+1] + np.power(B_e_i,-1)
        b_tmp[i, 0] = b_tmp[i, 0] - B_e_i/2
        b_tmp[i+1, 0] = b_tmp[i+1, 0] - B_e_i/2
    
    A_tmp[0, 0] = 1
    A_tmp[0, 1] = 0
    A_tmp[N_tmp-1, N_tmp-1] = 1
    A_tmp[N_tmp-1, N_tmp-2] = 0
    b_tmp[0, 0] = 1 # u(0) = 1
    b_tmp[N_tmp-1, 0] = 1 # u(1) = 1
    
    u_tmp = linalg.inv(A_tmp).dot(b_tmp)
    
    print("u #", iteration, u_tmp)
    
    str_plot = "Iteration " + str(iteration)
    
    plt.figure(1)
    plt.plot(x_tmp, u_tmp, label=str_plot)
    plt.legend()
    
    # calculate the error estimator again through the mesh
    xe_tmp = np.zeros((N_tmp-3))
    r2_tmp = np.zeros((N_tmp-3))
    R2_tmp = np.zeros((N_tmp-3))
    etaK_tmp = np.zeros((N_tmp-3))
    
    for i in range(1,N_tmp-2):
        B_e_i = x_tmp[i+1] - x_tmp[i] # jacobian, in this case the h_K
        # xe[i-1] = (x[i+1] - x[i])/2 + x[i]
        r2_tmp[i-1] = B_e_i
        R2_tmp[i-1] = np.power((-u_tmp[i-1,0] + 2*u_tmp[i,0] - u_tmp[i+1,0]),2) + np.power((-u_tmp[i,0] + 2*u_tmp[i+1,0] - u_tmp[i+2,0]),2)
        etaK_tmp[i-1] = np.sqrt(B_e_i**2*r2_tmp[i-1] + B_e_i*R2_tmp[i-1]) 
        # print(i, xe[i-1], etaK[i-1])
    
    print("etaK #", iteration, etaK_tmp)
    print("Media #", iteration, np.average(etaK_tmp))
    
    x_tmp_base = x_tmp
    etaK_tmp_base = etaK_tmp
    
    iteration += 1
    
plt.figure(1)
plt.xlabel("$x$")
plt.ylabel("$u(x)$")
plt.title("$-u^{''}(x) = -1$ with $u(0)=u(1)=1$")
plt.grid(True)
plt.show()