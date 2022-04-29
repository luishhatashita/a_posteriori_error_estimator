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
import pandas as pd

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

plt.figure(2)
plt.plot(x, u, label="Iteration 0")

indexdp_values = [0]

#-----------------------------------------------------------------------------
# error estimator
#-----------------------------------------------------------------------------

# initial hypothesis

tol_e = 0.01
max_iteration = 5

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

plt.figure(3)  
plt.plot(xe, r2, label="$||r||^2$")
plt.plot(xe, R2, label="$||R||^2$")
plt.bar(xe, etaK, width=0.7*dx, label="$\eta_K$")
plt.axis([x[0], x[N-1], 0, r2.max()+0.01])
plt.xlabel("$x$")
plt.title("Estimador de Erro")
plt.legend(loc="best")
plt.grid(True)

array_results = [[Ne, np.average(etaK), np.sum(etaK), 1, 0, Ne, np.average(etaK), np.sum(etaK), 1, 0]]
dict_results = {
        '0': {
            'x': x,
            'u': u,
            'xe': xe,
            'r2': r2,
            'R2': R2,
            'etaK': etaK,
        },
    }

#-----------------------------------------------------------------------------
# adaptive meshing
#-----------------------------------------------------------------------------

iteration = 1
x_tmp_base_adpt = x
etaK_tmp_base_adpt = etaK
N_tmp_adpt = N

while iteration <= max_iteration:
    
    # mesh refinement - adaptive
    ref_tmp_adpt = np.zeros((N_tmp_adpt-3))
    x_tmp_adpt = x_tmp_base_adpt
    for i in range(1,N_tmp_adpt-2):
        if etaK_tmp_base_adpt[i-1] > tol_e:
            index_tmp = list(x_tmp_adpt).index(x_tmp_base_adpt[i])
            x_new_tmp = (x_tmp_adpt[index_tmp+1] - x_tmp_adpt[index_tmp])/2 + x_tmp_adpt[index_tmp]
            x_tmp_adpt = np.insert(x_tmp_adpt, index_tmp+1, x_new_tmp)
            ref_tmp_adpt[i-1] = 1
    
    print(ref_tmp_adpt)
    
    dict_results[str(iteration-1)]['ref_adpt'] = ref_tmp_adpt
    
    # str_plot = "Iteration " + str(iteration-1)
    # x_ref_adpt = list(dict_results[str(iteration-1)]['x'])
    # x_ref_adpt.pop(0)
    # x_ref_adpt.pop()
    # x_ref_adpt.pop()
    
    # plt.figure(8)
    # plt.plot(x_ref_adpt, ref_tmp_adpt, "s", label=str_plot)
    # plt.legend()
    
    # mesh refinement - uniform
    N_tmp_uni = 2**iteration*Ne+1
    x_tmp_uni = np.linspace(0, Lx, N_tmp_uni)
   
    # break condition, all elements' errors are under the tolerance
    if np.size(x_tmp_adpt) == np.size(x_tmp_base_adpt):
        # dict_results[str(iteration)]['ref_adpt'] = ref_tmp_adpt
        break
    
    # solution of the new system equation - adaptive
    N_tmp_adpt = np.size(x_tmp_adpt)
    A_tmp_adpt = np.zeros((N_tmp_adpt, N_tmp_adpt))
    b_tmp_adpt = np.zeros((N_tmp_adpt, 1))
    
    for i in range(0, N_tmp_adpt-1):
        B_e_i = x_tmp_adpt[i+1] - x_tmp_adpt[i] # jacobian
        A_tmp_adpt[i, i] = A_tmp_adpt[i,i] + np.power(B_e_i,-1)
        A_tmp_adpt[i, i+1] = A_tmp_adpt[i, i+1] - np.power(B_e_i,-1)
        A_tmp_adpt[i+1, i] = A_tmp_adpt[i+1, i] - np.power(B_e_i,-1)
        A_tmp_adpt[i+1, i+1] = A_tmp_adpt[i+1, i+1] + np.power(B_e_i,-1)
        b_tmp_adpt[i, 0] = b_tmp_adpt[i, 0] - B_e_i/2
        b_tmp_adpt[i+1, 0] = b_tmp_adpt[i+1, 0] - B_e_i/2
    
    A_tmp_adpt[0, 0] = 1
    A_tmp_adpt[0, 1] = 0
    A_tmp_adpt[N_tmp_adpt-1, N_tmp_adpt-1] = 1
    A_tmp_adpt[N_tmp_adpt-1, N_tmp_adpt-2] = 0
    b_tmp_adpt[0, 0] = 1 # u(0) = 1
    b_tmp_adpt[N_tmp_adpt-1, 0] = 1 # u(1) = 1
    
    u_tmp_adpt = linalg.inv(A_tmp_adpt).dot(b_tmp_adpt)
    
    str_plot = "Iteration " + str(iteration)
    
    plt.figure(1)
    plt.plot(x_tmp_adpt, u_tmp_adpt, label=str_plot)
    plt.legend()
    
    # solution of the new system equation - uniform
    A_tmp_uni = np.zeros((N_tmp_uni, N_tmp_uni))
    b_tmp_uni = np.zeros((N_tmp_uni, 1))
     
    for i in range(0, N_tmp_uni-1):
        B_e_i = x_tmp_uni[i+1] - x_tmp_uni[i] # jacobian
        A_tmp_uni[i, i] = A_tmp_uni[i,i] + np.power(B_e_i,-1)
        A_tmp_uni[i, i+1] = A_tmp_uni[i, i+1] - np.power(B_e_i,-1)
        A_tmp_uni[i+1, i] = A_tmp_uni[i+1, i] - np.power(B_e_i,-1)
        A_tmp_uni[i+1, i+1] = A_tmp_uni[i+1, i+1] + np.power(B_e_i,-1)
        b_tmp_uni[i, 0] = b_tmp_uni[i, 0] - B_e_i/2
        b_tmp_uni[i+1, 0] = b_tmp_uni[i+1, 0] - B_e_i/2
     
    A_tmp_uni[0, 0] = 1
    A_tmp_uni[0, 1] = 0
    A_tmp_uni[N_tmp_uni-1, N_tmp_uni-1] = 1
    A_tmp_uni[N_tmp_uni-1, N_tmp_uni-2] = 0
    b_tmp_uni[0, 0] = 1 # u(0) = 1
    b_tmp_uni[N_tmp_uni-1, 0] = 1 # u(1) = 1
     
    u_tmp_uni = linalg.inv(A_tmp_uni).dot(b_tmp_uni)
    
    plt.figure(2)
    plt.plot(x_tmp_uni, u_tmp_uni, label=str_plot)
    plt.legend()
    
    # calculate the error estimator again through the mesh - adaptive
    xe_tmp_adpt = np.zeros((N_tmp_adpt-3))
    r2_tmp_adpt = np.zeros((N_tmp_adpt-3))
    R2_tmp_adpt = np.zeros((N_tmp_adpt-3))
    etaK_tmp_adpt = np.zeros((N_tmp_adpt-3))
    
    for i in range(1,N_tmp_adpt-2):
        B_e_i = x_tmp_adpt[i+1] - x_tmp_adpt[i] # jacobian, in this case the h_K
        xe_tmp_adpt[i-1] = (x_tmp_adpt[i+1] - x_tmp_adpt[i])/2 + x_tmp_adpt[i]
        r2_tmp_adpt[i-1] = B_e_i
        R2_tmp_adpt[i-1] = np.power((-u_tmp_adpt[i-1,0] + 2*u_tmp_adpt[i,0] - u_tmp_adpt[i+1,0]),2) + np.power((-u_tmp_adpt[i,0] + 2*u_tmp_adpt[i+1,0] - u_tmp_adpt[i+2,0]),2)
        etaK_tmp_adpt[i-1] = np.sqrt(B_e_i**2*r2_tmp_adpt[i-1] + B_e_i*R2_tmp_adpt[i-1]) 
    
    # calculate the error estimator again through the mesh - uniform
    # xe_tmp_uni = np.zeros((N_tmp_uni-3))
    r2_tmp_uni = np.zeros((N_tmp_uni-3))
    R2_tmp_uni = np.zeros((N_tmp_uni-3))
    etaK_tmp_uni = np.zeros((N_tmp_uni-3))

    for i in range(1,N_tmp_uni-2):
        B_e_i = x_tmp_uni[i+1] - x_tmp_uni[i] # jacobian, in this case the h_K
        # xe_tmp_uni[i-1] = (x_tmp_uni[i+1] - x_tmp_uni[i])/2 + x_tmp_uni[i]
        r2_tmp_uni[i-1] = B_e_i
        R2_tmp_uni[i-1] = np.power((-u_tmp_uni[i-1,0] + 2*u_tmp_uni[i,0] - u_tmp_uni[i+1,0]),2) + np.power((-u_tmp_uni[i,0] + 2*u_tmp_uni[i+1,0] - u_tmp_uni[i+2,0]),2)
        etaK_tmp_uni[i-1] = np.sqrt(B_e_i**2*r2_tmp_uni[i-1] + B_e_i*R2_tmp_uni[i-1]) 
    
    # save compilation of results to array_results
    # the results included are as follows for both meshes:
    # 1 - number of elements
    # 2 - average of individual error estimates
    # 3 - total error (sum of individual error estimates)
    # 4 - absolute error drop of total error (3, always comparing to the iteration 0)
    # 5 - relative error of total error (3, comparing to previous iteration)
    Ne_tmp_uni = N_tmp_uni-1
    avrg_tmp_uni = np.average(etaK_tmp_uni)
    sum_tmp_uni = np.sum(etaK_tmp_uni)
    error_drop_rel_tmp_uni = abs(array_results[iteration-1][2] - sum_tmp_uni)/array_results[iteration-1][2] # relative error to previous iteration
    error_drop_abs_tmp_uni = 1 - abs(array_results[0][2] - sum_tmp_uni)/array_results[0][2] # absolute error comparing to iteration 0
    Ne_tmp_adpt = N_tmp_adpt-1
    avrg_tmp_adpt = np.average(etaK_tmp_adpt)
    sum_tmp_adpt = np.sum(etaK_tmp_adpt)
    error_drop_rel_tmp_adpt = abs(array_results[iteration-1][7] - sum_tmp_adpt)/array_results[iteration-1][7] # relative error to previous iteration
    error_drop_abs_tmp_adpt = 1 - abs(array_results[0][7] - sum_tmp_adpt)/array_results[0][7] # absolute error comparing to iteration 0
    
    array_results.append([Ne_tmp_uni, avrg_tmp_uni, sum_tmp_uni, error_drop_abs_tmp_uni, error_drop_rel_tmp_uni, Ne_tmp_adpt, avrg_tmp_adpt, sum_tmp_adpt, error_drop_abs_tmp_adpt, error_drop_rel_tmp_adpt])
    
    # save intermidiate results to dict_results
    # the choice of a dictionary was mainly due to the lack of uniformity of the
    #   results (size(x) is different than size(xe), for example)
    dict_results[str(iteration)] = {
            'x': x_tmp_adpt,
            'u': u_tmp_adpt,
            'xe': xe_tmp_adpt,
            'r2': r2_tmp_adpt,
            'R2': R2_tmp_adpt,
            'etaK': etaK_tmp_adpt,
        }
    
    # reset of comparing arrays, necessary to be able to include new mesh nodes
    #   without messing with what was already evaluated
    x_tmp_base_adpt = x_tmp_adpt
    etaK_tmp_base_adpt = etaK_tmp_adpt
    
    indexdp_values.append(iteration)
    
    iteration += 1

# create results dataframe (df_results)
columndp_values = ["Uni N elmt", "Uni Avrg etaK", "Uni Sum etaK", "Uni Error Drop Abs", "Uni Error Drop Rel", "Adpt N elmt", "Adpt Avrg etaK", "Adpt Sum etaK", "Adpt Error Drop Abs", "Adpt Error Drop Rel"]
df_results = pd.DataFrame(data = array_results, index = indexdp_values, columns = columndp_values)

print(df_results)

plt.figure(1)
plt.xlabel("$x$")
plt.ylabel("$u(x)$")
plt.title("$-u^{'''}(x) = -1$ with $u(0)=u(1)=1$ - Adaptive Mesh")
plt.grid(True)

plt.figure(2)
plt.xlabel("$x$")
plt.ylabel("$u(x)$")
plt.title("$-u^{'''}(x) = -1$ with $u(0)=u(1)=1$ - Uniform Mesh")
plt.grid(True)

np_array_results = np.array(array_results)

plt.figure(4)  
plt.plot(df_results['Adpt N elmt'], df_results['Adpt Avrg etaK'], "s", label="Adaptive mesh")
plt.plot(df_results['Uni N elmt'], df_results['Uni Avrg etaK'], "^", label="Uniform mesh")
plt.xlabel("Number of Elements")
plt.ylabel("Average $\eta_K$")
plt.title("Mesh Comparison - Average $\eta_K$")
plt.legend(loc="best")
plt.grid(True)

plt.figure(5)  
plt.plot(df_results['Adpt N elmt'], df_results['Adpt Sum etaK'], "s", label="Adaptive mesh")
plt.plot(df_results['Uni N elmt'], df_results['Uni Sum etaK'], "^", label="Uniform mesh")
# plt.xscale("log")
plt.xlabel("Number of Elements")
plt.ylabel("Sum $\eta_K$")
plt.title("Mesh Comparison - Total $\eta_K$ (Sum)")
plt.legend(loc="best")
plt.grid(True)

width = 0.25

plt.figure(6)
plt.subplot(1, 2, 1)  
# plt.plot(df_results.index.values, np_array_results[:,7], "s", label="Adaptive mesh")
# plt.plot(df_results.index.values, np_array_results[:,3], "^", label="Uniform mesh")
plt.bar(df_results.index.values - width/2, df_results['Adpt Error Drop Abs'], width, label="Adaptive mesh")
plt.bar(df_results.index.values + width/2, df_results['Uni Error Drop Abs'], width, label="Uniform mesh")
plt.xlabel("Iteration")
plt.ylabel("Absolute Error for Sum $\eta_K$")
# plt.title("Absolute Error Drop Total $\eta_K$ (Sum)")
# plt.legend(loc="best")
# plt.grid(True)

plt.subplot(1, 2, 2)   
plt.plot(df_results['Adpt N elmt'], df_results['Adpt Error Drop Abs'], "s", label="Adaptive mesh")
plt.plot(df_results['Uni N elmt'], df_results['Uni Error Drop Abs'], "^", label="Uniform mesh")
plt.xscale("log")
plt.xlabel("Number or Elements")
# plt.ylabel("Absolute Error for Sum $\eta_K$")
plt.suptitle("Absolute Error Drop Total $\eta_K$ (Sum)")
plt.legend(loc="best")
plt.grid(True)

plt.figure(7)
plt.subplot(1, 2, 1)  
# plt.plot(df_results.index.values, np_array_results[:,7], "s", label="Adaptive mesh")
# plt.plot(df_results.index.values, np_array_results[:,3], "^", label="Uniform mesh")
plt.bar(df_results.iloc[1:, :].index.values - width/2, df_results.iloc[1:, :]['Adpt Error Drop Rel'], width, label="Adaptive mesh")
plt.bar(df_results.iloc[1:, :].index.values + width/2, df_results.iloc[1:, :]['Uni Error Drop Rel'], width, label="Uniform mesh")
plt.xlabel("Iteration")
plt.ylabel("Relative Error for Sum $\eta_K$")
# plt.title("Absolute Error Drop Total $\eta_K$ (Sum)")
plt.legend(loc="best")
# plt.grid(True)

plt.subplot(1, 2, 2)   
plt.plot(df_results.iloc[1:, :]['Adpt N elmt'], df_results.iloc[1:, :]['Adpt Error Drop Rel'], "s", label="Adaptive mesh")
plt.plot(df_results.iloc[1:, :]['Uni N elmt'], df_results.iloc[1:, :]['Uni Error Drop Rel'], "^", label="Uniform mesh")
# plt.xscale("log")
# plt.yscale("log")
plt.xlabel("Number or Elements")
# plt.ylabel("Absolute Error for Sum $\eta_K$")
plt.suptitle("Relative Error Total $\eta_K$ (Sum)")
# plt.legend(loc="best")
plt.grid(True)

# preparation for error estimate analysis behavior
# histogram plot

x_end_adpt = list(dict_results[max(dict_results.keys())]['x'])
x_end_adpt.pop(0)
x_end_adpt.pop()
x_end_adpt.pop()
ref_end_adpt = [0 for i in range(np.size(x_end_adpt))]

for iteration, information in dict_results.items():
    # print(iteration)
    x_end_tmp = list(information['x'])
    x_end_tmp.pop(0)
    x_end_tmp.pop()
    x_end_tmp.pop()
    for x in x_end_tmp:
        index_x = x_end_tmp.index(x)
        index_x_end = x_end_adpt.index(x)
        ref_end_adpt[index_x_end] += list(information['ref_adpt'])[index_x]

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

number_bins_hist = 4
chunks_x_end_adpt = list(chunks(x_end_adpt, int(np.size(x_end_adpt)/number_bins_hist)))
agregated_ref = list(np.zeros(np.size(chunks_x_end_adpt)))

for idx, chunk in enumerate(chunks_x_end_adpt):
    for x in chunk:
        index_x_end = x_end_adpt.index(x)
        # index_x_chunk = chunks_x_end_adpt.index(x)
        # print(index_x_chunk)
        agregated_ref[idx] += ref_end_adpt[index_x_end]

print(agregated_ref)

# plt.figure(8)
# plt.xlabel("$x$")
# plt.ylabel("Refinement Count")
# plt.title("Refinement Count - Adaptive Mesh")
# plt.grid(True)

# plt.figure(9)
# plt.plot(x_end_adpt, ref_end_adpt, "s")
# plt.xlabel("$x$")
# plt.ylabel("Refinement Count")
# plt.title("Refinement Count - Adaptive Mesh")

plt.figure(10)
plt.plot(agregated_ref)
plt.xlabel("$x$")
plt.ylabel("Refinement Count")
plt.title("Refinement Count - Adaptive Mesh")

print(ref_end_adpt)
plt.show()