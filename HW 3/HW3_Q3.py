# QUANTITATIVE MACROECONOMICS: 
# HWK 3 - VALUE FUNCTION ITERATION
    # Santiago Albarracin
    # Victor Caballero
    # Lorenzo de Carolis

import numpy as np
import matplotlib.pyplot as plt
import quantecon as qe
import math

# =============================================================================
# EXERCISE 3: Chebyshev
# =============================================================================

# parameters 
theta = 0.679
beta = 0.988
delta = 0.013
h = 1
kappa = 5.24
nu = 2.0

kSS = ((1/beta-1+delta)/(1-theta))**(-1/theta)

print("The steady state is "+str(kSS))

qe.tic()
nk = 200
k = np.linspace(1,2*kSS,nk)

V0 = np.zeros(nk)

def f(k,h):
    return k**(1-theta)*h**theta

def u(c,h):
    return math.log(c) - kappa*h**(1+1/nu)/(1+1/nu)

def re_M(k1,k2):
    c = f(k1,h) + (1-delta)*k1 - k2
    if c>0:
        return u(c,h)
    else:
        return -1000000

M = np.empty([nk,nk])
i=0
while i<=nk-1:
    j=0
    while j<=nk-1:
        M[i,j] = re_M(k[i],k[j]) 
        j = j+1
    i = i+1

# Chebyshev 

def cheb_nodes(x,a,b):
    k = []
    z = []
    for j in range(1,x+1):   
        z_k=-np.cos(np.pi*(2*j-1)/(2*x))   
        k_ch=(z_k+1)*((b-a)/2)+a  
        z.append(z_k)
        k.append(k_ch)
    return np.array(z), np.array(k)

#
#
#
# We tried different thing but nothing from this point on worked.
# We couldn't complete Q3
#
#
#


T = qe.toc()

plt.plot(k, V, label='V(k)', color="red")    
plt.xlabel('Capital (k)')
plt.ylabel('Value (Vk)')
plt.show()

print("Time elapse in seconds: " + str(T))
print("Number of iterations: " + str(j))
