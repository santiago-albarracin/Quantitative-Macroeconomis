#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 20:05:43 2020

@author: santiago
"""
#Importing packages
import pandas as pd
import numpy as np
import mpmath as mp
import sympy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import LinearConstraint, minimize, fsolve
import seaborn as sb

#Parameters
A=1
rho=1.1
k=0.2
w=20
gamma=0.9
i_0=0.2
N=1

x_points, y_points= np.meshgrid(np.arange(0,101,1), np.arange(0,101,1), indexing="xy")
x_points = x_points * 0.01
y_points = y_points * 0.01
grid = np.array([x_points, y_points])
grid_2 = grid[0]

results = grid.copy()
output = grid_2.copy()
wel = grid_2.copy()
infec = grid_2.copy()
det = grid_2.copy()

for z in np.arange(0,101,1):
    for y in np.arange(0,101,1):
        c_tw, B = grid[:, z, y]
        Y = lambda H_f, H_nf:(A * (H_f**((rho-1)/rho))+c_tw*A*(H_nf**((rho-1)/rho)))**(rho/(rho-1))
        I = lambda H_f: B * ((i_0*H_f)/N) * H_f
        D = lambda H_f: (1-gamma) * I(H_f)
        
        def solution_func(x):
            H_f = x[0]
            H_nf = x[1]
            return -1*(Y(H_f, H_nf) - k * H_f - k * H_nf - w * D(H_f))
       
        #1st guess
        x0 = np.array([0.5, 0.5])
        
        #Such that (s.t.)
        g_1 = np.array([-1, -1])
        g_2 = np.array([-N])
        s_t = [(0, 1), (0, 1)]
        
        #Solving
        const={"type" : "ineq", "fun" : lambda x: g_1 @ x - g_2}
    
        solution = minimize(solution_func, x0, method="SLSQP", bounds=s_t, constraints=const)
        results[:, z, y] = solution.x 
        
        #In order to plot the required figures
        output[z, y] = (A * (results[0, z, y] * ((rho - 1) / rho)) + c_tw * A * (results[1, z, y] ** ((rho - 1) / rho))) ** (rho / (rho - 1))
        infec[z, y] = B * ((i_0*results[0, z, y]) / N) * results[0, z, y]
        det[z, y] = (1 - gamma) * infec[z, y]
        wel[z, y] = output[z, y] - k * results[0, z, y] - k * results[1, z, y] - w * det[z, y]
    
        H_f = results[0]
        H_nf = results[1]
        
        H = np.array(H_f + H_nf)
        H_f_H = np.array(H_f/H)
        
        #Ploting the results
        plt.figure(figsize=(10,6))
        sb.set(font_scale=1.2)
        sb.heatmap(H_f,cmap="Greens" , cbar_kws={'label': ''})
        plt.xlabel('c(TW)')
        plt.ylabel('beta(HC)')
        plt.title('Optimal allocation of H')
        plt.show()

        plt.figure(figsize=(10,6))
        sb.set(font_scale=1.2)
        sb.heatmap(H_f,cmap="Green"  , cbar_kws={'label': ''}, vmin=0, vmax=1)
        plt.xlabel('c(TW)')
        plt.ylabel('beta(HC)')
        plt.title('Optimal allocation of H_f')
        plt.show()
        
        plt.figure(figsize=(10,6))
        sb.set(font_scale=1.2)
        sb.heatmap(H_nf,cmap="Green"  , cbar_kws={'label': ''}, vmin=0, vmax=1)
        plt.xlabel('c(TW)')
        plt.ylabel('beta(HC)')
        plt.title('Optimal allocation of H_nf')
        plt.show()
        
        plt.figure(figsize=(10,6))
        sb.set(font_scale=1.2)
        sb.heatmap(H_f_H,cmap="Green"  , cbar_kws={'label': ''}, vmin=0, vmax=1)
        plt.xlabel('c(TW)')
        plt.ylabel('beta(HC)')
        plt.title('Optimal allocation of H_f_H')
        plt.show()
        
        plt.figure(figsize=(10,6))
        sb.set(font_scale=1.2)
        sb.heatmap(output,cmap="Greens"  , cbar_kws={'label': ''}, vmin=0, vmax=1)
        plt.xlabel('c(TW)')
        plt.ylabel('beta(HC)')
        plt.title('Optimal allocation of Output')
        plt.show()
        
        plt.figure(figsize=(10,6))
        sb.set(font_scale=1.2)
        sb.heatmap(wel,cmap="Greens"  , cbar_kws={'label': ''}, vmin=0, vmax=1)
        plt.xlabel('c(TW)')
        plt.ylabel('beta(HC)')
        plt.title('Optimal allocation of Welfare')
        plt.show()
        
        plt.figure(figsize=(10,6))
        sb.set(font_scale=1.2)
        sb.heatmap(infec,cmap="Greens"  , cbar_kws={'label': ''}, vmin=0, vmax=1)
        plt.xlabel('c(TW)')
        plt.ylabel('beta(HC)')
        plt.title('Optimal allocation of Amount of infections')
        plt.show()
        
        plt.figure(figsize=(10,6))
        sb.set(font_scale=1.2)
        sb.heatmap(det,cmap="Greens"  , cbar_kws={'label': ''}, vmin=0, vmax=1)
        plt.xlabel('c(TW)')
        plt.ylabel('beta(HC)')
        plt.title('Optimal allocation of Deaths')        
        plt.show()

        
        
        
        
        