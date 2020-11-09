#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 23:22:21 2020

@author: Santiago Albarracin
"""
import math
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve
import pandas as pd


#=============### QUESTION 1 ###==============================================#

# PART (a) #-------------------------------------------------------------------

# Parameters
theta = 0.67
h = 0.31
y_1=1 #normslization assumption

# Given the ratios information in HW2:
i_1=0.25*y_1
k_1=4*y_1

# Computations to obtain the values of the rest of the variables in the SS
c_1=y_1-i_1
z_1=(y_1/(k_1**(1-theta)))**(1/theta)/h
delta=i_1/k_1
beta=1/((1-theta)*(h*z_1)**theta*k_1**(-theta)-delta+1)

# Output of the first SS:
SS_1={'Variable name': ['y_1', 'k_1', 'c_1', 'i_1', 'z_1', 'h_1', 'beta', 'theta', 'delta'], 'Values at  SS_1': [y_1, k_1, c_1, i_1, z_1, h, beta, theta, delta]}
data_SS_1 = pd.DataFrame(SS_1)
print(data_SS_1)

# PART (b) #-------------------------------------------------------------------

# New z:
z_2 = z_1*2
k_2 = (((1/beta)+(delta-1))/((1-theta)*((h*z_2)**theta)))**(-(1/theta))
i_2 = delta*k_2
y_2=k_2**(1-theta)*(z_2*h)**theta
c_2=y_2-i_2

# Output of the first and second SS:
SS_1_2={'Variable name': ['y', 'k', 'c', 'i', 'z', 'h', 'beta', 'theta', 'delta'], 'Values at  SS_1': [y_1, k_1, c_1, i_1, z_1, h, beta, theta, delta], 'Values at  SS_2': [y_2, k_2, c_2, i_2, z_2, h, beta, theta, delta]}
data_SS_1_2 = pd.DataFrame(SS_1_2)
print(data_SS_1_2)

# PART (c) #-------------------------------------------------------------------

# Values of the variables at the initial time (when shock hits):
c = 1.024
t = 0
c_t = []
c_t.append(c)
i_t = []
i_t.append(i_1)
y_t = []
y_t.append(y_1)
k_t = []
k_t.append(k_1)

# Computing the transitioni loop until k fids its SS (k=k_2):
while k_2 - k_t[t] > 0.1:
    k_t.append(k_t[t]**(1-theta)*(z_2*h)**theta + (1-delta)*k_t[t]- c_t[t])
    c_t.append(((1-theta)*(h*z_2)**theta*k_t[t+1]**(-theta) + (1 - delta))*c_t[t]*beta)
    y_t.append(k_t[t+1]**(1-theta) * (z_2*h)**theta)
    i_t.append(y_t[t+1] - c_t[t+1])
    t += 1

# Plot time paths:
plt.figure()
plt.subplot(221)
plt.plot(c_t)
plt.title('Consumption path till SS')
plt.ylabel('C')
plt.xlabel('Time')
plt.subplot(222)
plt.plot(k_t)
plt.title('Capital path till SS')
plt.ylabel('K')
plt.xlabel('Time')
plt.subplot(223)
plt.plot(y_t)
plt.title('Output path till SS')
plt.ylabel('Y')
plt.xlabel('Time')
plt.subplot(224)
plt.plot(i_t)
plt.title('Investment path till SS')
plt.ylabel('I')
plt.xlabel('Time')
plt.subplots_adjust(top=2, bottom=0.08, left=0, right=2, hspace=0.3, wspace=0.2)
plt.show()

print('The economy reaches the SS after '+str(t)+' periods')

# PART (d) #-------------------------------------------------------------------

# Values of the variables at the initial time (when first shock hits):
c = 1.024
t = 0
c_t = []
c_t.append(c)
i_t = []
i_t.append(i_1)
y_t = []
y_t.append(y_1)
k_t = []
k_t.append(k_1)

# Two loops: 
    # The fisrt one explaining how the economy behaves when not expecting
        #another shock;
    # The second one explainig how the economy behaves when in period t=10 the
        #new unexpected shock hits.

while (k_2 - k_t[t] > 0.1) and (t < 10):
    k_t.append(k_t[t]**(1-theta)*(z_2*h)**theta + (1-delta)*k_t[t]- c_t[t])
    c_t.append(((1-theta)*(h*z_2)**theta*k_t[t+1]**(-theta) + (1 - delta))*c_t[t]*beta)
    y_t.append(k_t[t+1]**(1-theta) * (z_2*h)**theta)
    i_t.append(y_t[t+1] - c_t[t+1])
    t += 1
while (k_t[t] - k_1 > 0.1) and (t > 9):
    if t == 10:
        c_t.append(0.93) # ???????????????????????????????????????????????????
        k_t.append(k_t[t]**(1-theta)*(z_1*h)**theta + (1-delta)*k_t[t]- c_t[t])
        y_t.append(k_t[t+1]**(1-theta) * (z_1*h)**theta)
        i_t.append(y_t[t+1] - c_t[t+1])
        t += 1
    else:
        k_t.append(k_t[t]**(1-theta)*(z_1*h)**theta + (1-delta)*k_t[t]- c_t[t])
        c_t.append(((1-theta)*(h*z_1)**theta*k_t[t+1]**(-theta) + (1 - delta))*c_t[t]*beta)
        y_t.append(k_t[t+1]**(1-theta) * (z_1*h)**theta)
        i_t.append(y_t[t+1] - c_t[t+1])
        t += 1

# Plot time paths:
plt.figure()
plt.subplot(221)
plt.plot(c_t)
plt.title('Consumption path till SS')
plt.ylabel('C')
plt.xlabel('Time')
plt.subplot(222)
plt.plot(k_t)
plt.title('Capital path till SS')
plt.ylabel('K')
plt.xlabel('Time')
plt.subplot(223)
plt.plot(y_t)
plt.title('Output path till SS')
plt.ylabel('Y')
plt.xlabel('Time')
plt.subplot(224)
plt.plot(i_t)
plt.title('Investment path till SS')
plt.ylabel('I')
plt.xlabel('Time')
plt.subplots_adjust(top=2, bottom=0.08, left=0, right=2, hspace=0.3, wspace=0.2)
plt.show()

print('The economy reaches the SS after '+str(t)+' periods')


