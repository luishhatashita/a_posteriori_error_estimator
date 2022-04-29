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
# from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import time

# ----------------------------------------------------------------------------
#                               Observations
# ----------------------------------------------------------------------------
# 1 - This script diverges from the first one by the usage of another error
# estimator
# 2 - However, the variables will still be the same to reduce the effort 
# necessary for implementation
#
# ----------------------------------------------------------------------------
#                       final results matrix scheme
# ----------------------------------------------------------------------------
# results_matrix = {
#        Ne: {
#            tol_e: {
#                final_iteration: ... ,
#                elapsed_time: ... ,
#                etaK: ... ,
#                total_error: ... ,
#                },
#            },
#     }
# ----------------------------------------------------------------------------

ne_array = [5, 10, 25, 50, 100]
tol_e_array = [0.01, 0.005, 0.001, 0.0005, 0.0001]

# initialize eval_matrix dictionary
eval_matrix = {}
eval_array = []

for Ne in ne_array:
    for tol_e in tol_e_array:
        
        success_flag = False
        
        eval_matrix.setdefault(str(Ne), {}).setdefault(str(tol_e), {})
        
        print(f"\n--- Ne: {Ne} --- tol_e: {tol_e} ---")
        
        # initial hypothesis
        
        acc_time_all = 0
        acc_time_uni = 0
        
        tic_all = time.perf_counter()
        
        # Ne = 5
        
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
        
        toc_all = time.perf_counter()
        
        acc_time_all += toc_all - tic_all
        
        indexdp_values = [0]
        
        #-----------------------------------------------------------------------------
        # error estimator
        #-----------------------------------------------------------------------------
        
        tic_all = time.perf_counter()
        
        # initial hypothesis
        
        # tol_e = 0.001
        max_iteration = 100
        
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
            etaK[i-1] = np.sqrt(B_e_i**4*r2[i-1] + B_e_i**3*R2[i-1]) # NEW etaK, h^4*r^2 and h^3*R^2
        
        toc_all = time.perf_counter()
        
        acc_time_all += toc_all - tic_all
        
        array_results = [[Ne, np.average(etaK), np.sum(etaK), 1, 0]]
        
        #-----------------------------------------------------------------------------
        # adaptive meshing
        #-----------------------------------------------------------------------------
        
        iteration = 1
        
        etaK_tmp_uni = etaK
        sum_tmp_uni = np.sum(etaK_tmp_uni)
        Ne_tmp_uni = Ne
        
        while iteration <= max_iteration:
            
            tic_uni = time.perf_counter()
            
            if max(etaK_tmp_uni) < tol_e:
                final_it = iteration - 1
                print(f"Iteracao final: {final_it}")
                toc_uni = time.perf_counter()
                acc_time_uni += toc_uni - tic_uni
                total_time = acc_time_all + acc_time_uni
                print(f"Total error: {sum_tmp_uni}")
                print(f"Elapsed time: {total_time}s")
                
                success_flag = True
                
                eval_matrix[str(Ne)][str(tol_e)] = {
                                'final_iteration': final_it,
                                'elapsed_time': total_time,
                                'etaK': etaK_tmp_uni,
                                'total_error': sum_tmp_uni,
                                }
                
                eval_array.append([Ne, tol_e, final_it, Ne_tmp_uni, sum_tmp_uni, total_time, total_time/sum_tmp_uni])
                
                break
            
            # mesh refinement - uniform
            N_tmp_uni = 2**iteration*Ne+1
            x_tmp_uni = np.linspace(0, Lx, N_tmp_uni)
            
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
            
            toc_uni = time.perf_counter()
            acc_time_uni += toc_uni - tic_uni
            
            tic_uni = time.perf_counter()
            
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
                etaK_tmp_uni[i-1] = np.sqrt(B_e_i**4*r2_tmp_uni[i-1] + B_e_i**3*R2_tmp_uni[i-1]) # NEW etaK, h^4*r^2 and h^3*R^2
            
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
            # error_drop_rel_tmp_uni = abs(array_results[iteration-1][2] - sum_tmp_uni)/array_results[iteration-1][2] # relative error to previous iteration
            # error_drop_abs_tmp_uni = 1 - abs(array_results[0][2] - sum_tmp_uni)/array_results[0][2] # absolute error comparing to iteration 0
            
            # array_results.append([Ne_tmp_uni, avrg_tmp_uni, sum_tmp_uni, error_drop_abs_tmp_uni, error_drop_rel_tmp_uni])
                
            # indexdp_values.append(iteration)
            
            # if max(etaK_tmp_uni) < tol_e:
            #     print(f"Iteracao final: {iteration}")
            #     toc_uni = time.perf_counter()
            #     acc_time_uni += toc_uni - tic_uni
            #     total_time = acc_time_all + acc_time_uni
            #     print(f"Total error: {sum_tmp_uni}")
            #     print(f"Elapsed time: {total_time}s")
                
            #     success_flag = True
                
            #     eval_matrix[str(Ne)][str(tol_e)] = {
            #                     'final_iteration': iteration,
            #                     'elapsed_time': total_time,
            #                     'etaK': etaK_tmp_uni,
            #                     'total_error': sum_tmp_uni,
            #                     }
                
            #     eval_array.append([Ne, tol_e, iteration, Ne_tmp_uni, sum_tmp_uni, total_time, total_time/sum_tmp_uni])
                
            #     break
            
            iteration += 1
            
            toc_uni = time.perf_counter()
            acc_time_uni += toc_uni - tic_uni
        
        # create results dataframe (df_results)
        # columndp_values = ["Uni N elmt", "Uni Avrg etaK", "Uni Sum etaK", "Uni Error Drop Abs", "Uni Error Drop Rel"]
        # df_results = pd.DataFrame(data = array_results, index = indexdp_values, columns = columndp_values)
        
        # print(df_results)
        
        if success_flag == False:
            print(f"\tMax iteration: {max_iteration} INSUFICIENT for >>> Ne: {Ne} --- tol_e: {tol_e} <<<")

columns_eval_array = ['Initial Ne', 'Individual tol_e', 'Final Iteration', 'Final Ne', 'Final Total Error', 'Elapsed Time', 'Cost_Benefit']

df = pd.DataFrame(data=eval_array, columns=columns_eval_array)

fig1 = plt.figure(1, figsize=plt.figaspect(0.5))
# plt.title('Comparison of Initial Parameters')

x = np.array(eval_array)[:,0] # Initial Ne
y = np.array(eval_array)[:,1] # Individual tol_e
z1 = np.array(eval_array)[:,4] # Total error > sum_tmp_uni

# Subplot 1 - Scatter data
ax1 = fig1.add_subplot(1, 2, 1, projection='3d')
ax1.scatter(x, y, z1, color='b', marker='s')
ax1.set_xlabel('Initial Number of Elements')
ax1.set_ylabel('Indivual Error Tolerance')
ax1.set_zlabel('Total Error - $\Sigma{\eta_K}$')

# Subplot 2 - Approximated surface
ax1 = fig1.add_subplot(1, 2, 2, projection='3d')
ax1.plot_trisurf(x, y, z1)
ax1.set_xlabel('Initial Number of Elements')
ax1.set_ylabel('Indivual Error Tolerance')
ax1.set_zlabel('Total Error - $\Sigma{\eta_K}$')

fig2 = plt.figure(2, figsize=plt.figaspect(0.5))
# plt.title('Comparison of Initial Parameters')

z2 = np.array(eval_array)[:,5] # Elapsed time > total_time

# Subplot 1 - Scatter data
ax2 = fig2.add_subplot(1, 2, 1, projection='3d')
ax2.scatter(x, y, z2, color='b', marker='s')
ax2.view_init(25, 60)
ax2.set_xlabel('Initial Number of Elements')
ax2.set_ylabel('Indivual Error Tolerance')
ax2.set_zlabel('Elapsed Time')

# Subplot 2 - Approximated surface
ax2 = fig2.add_subplot(1, 2, 2, projection='3d')
ax2.plot_trisurf(x, y, z2)
ax2.view_init(25, 60)
ax2.set_xlabel('Initial Number of Elements')
ax2.set_ylabel('Indivual Error Tolerance')
ax2.set_zlabel('Elapsed Time')

fig3 = plt.figure(3, figsize=plt.figaspect(0.5))
# plt.title('Comparison of Initial Parameters')

z3 = np.array(eval_array)[:,6] # Cost/Benefit ratio > total_time/sum_tmp_uni

# Subplot 1 - Scatter data
ax3 = fig3.add_subplot(1, 2, 1, projection='3d')
ax3.scatter(x, y, z3, color='b', marker='s')
ax3.view_init(25, 60)
ax3.set_xlabel('Initial Number of Elements')
ax3.set_ylabel('Indivual Error Tolerance')
ax3.set_zlabel('Cost/Benefit')

# Subplot 2 - Approximated surface
ax3 = fig3.add_subplot(1, 2, 2, projection='3d')
ax3.plot_trisurf(x, y, z3)
ax3.view_init(25, 60)
ax3.set_xlabel('Initial Number of Elements')
ax3.set_ylabel('Indivual Error Tolerance')
ax3.set_zlabel('Cost/Benefit')

plt.show()