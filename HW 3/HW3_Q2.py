# QUANTITATIVE MACROECONOMICS: 
# HWK 3 - VALUE FUNCTION ITERATION
    # Santiago Albarracin
    # Victor Caballero
    # Lorenzo de Carolis
    
import numpy as np
import matplotlib.pyplot as plt
import math
import quantecon as qe

# =============================================================================
# EXERCISE 2: VFI with continous labour choice
# =============================================================================

theta = 0.679
beta = 0.988
delta = 0.013
h = 1
k = 1
kappa = 5.24
nu = 2.0

def m(ki, kj):
    grid = 200
    kappa = 5.24
    theta = 0.679
    delta = 0.013
    nu = 2.0
    h = np.linspace(0.01, 1, grid)
    u = np.zeros(grid)
    c = 0
    for i in range(grid):
        c = ki**(1-theta)*h[i]**theta+(1-delta)-kj
        if c >= 0:
            u[i] = math.log(c) - kappa*(h[i]**(1+1/nu)/(1+1/nu))
        else:
            u[i] = -1000000
        
    return max(u)


#%%(2.a)----------------------------------------------------------------#


qe.tic()
k_min = 0.01
k_max = 2
nk = 100

k = np.linspace(k_min, k_max, nk)

V = []
V.append(np.zeros(nk))
V.append(np.zeros(nk))
k_dec = []
k_dec.append(np.zeros(nk))
kj_dec = []
kj_dec.append(np.zeros(nk))


M = [0]*nk
Chi = [0]*nk
for i in range(nk):
    M[i] = [0] * nk
    Chi[i] = [0] * nk
    
for i in range(nk):
    for j in range(nk):
        M[i][j] = m(k[i], k[j])

for i in range(nk):
    for j in range(nk):
        Chi[i][j] = M[i][j] + beta* V[-1][j]

V_new = [0]*nk
k_dec_new = [0]*nk
for i in range(nk):
    V_new[i]=max(Chi[i])

V.append(V_new)

while (np.linalg.norm(np.asarray(V[-1]) - np.asarray(V[-2])) > 0.01):
    for i in range(nk):
        for j in range(nk):
            Chi[i][j] = M[i][j] + beta* V[-1][j]
    
    V_new = [0]*nk
    k_dec_new = [0]*nk
    for i in range(nk):
        V_new[i]=max(Chi[i])

    V.append(V_new)
    
T = qe.toc()

# Plotting:
plt.plot(k, V[-1], label='V(k)', color="red")    
plt.xlabel('Capital (k)')
plt.ylabel('Value (Vk)')
plt.show()

print("Time elapse in seconds: " + str(T))
print("Number of iterations: " + str(j))


#%%(2.b)----------------------------------------------------------------#

qe.tic()
k_min = 0.01
k_max = 2
nk = 100

k = np.linspace(k_min, k_max, nk) 

V = []
V.append(np.zeros(nk))
V.append(np.zeros(nk))
k_dec = []
k_dec.append(np.zeros(nk))

kj_dec = []
kj_dec.append(np.zeros(nk))

M = [0]*nk
Chi = [0]*nk
for i in range(nk):
    M[i] = [0] * nk
    Chi[i] = [0] * nk
    
for i in range(nk):
    for j in range(nk):
        M[i][j] = m(k[i], k[j])

for i in range(nk):
    for j in range(nk):
        Chi[i][j] = M[i][j] + beta* V[-1][j]

V_new = [0]*nk
k_dec_new = [0]*nk
for i in range(nk):
    V_new[i]=max(Chi[i])

V.append(V_new)

while (np.linalg.norm(np.asarray(V[-1]) - np.asarray(V[-2])) > 0.01):
    for i in range(nk):
        for j in range(nk):
            if k[j]>=k_dec[-1][j]:
                Chi[i][j] =M[i][j] + beta* V[-1][j]

    V_new = [0]*nk
    k_dec_new = [0]*nk
    for i in range(nk):
        V_new[i]=max(Chi[i])

    V.append(V_new)

T = qe.toc()

plt.plot(k, V[-1], label='V(k)', color="blue")    
plt.xlabel('Capital (k)')
plt.ylabel('Value (Vk)')
plt.show()

print("Time elapse in seconds: " + str(T))
print("Number of iterations: " + str(j))


#%%(2.c)----------------------------------------------------------------#

qe.tic()
Cond = True
k_min = 0.01
k_max = 2
nk = 100

k = np.linspace(k_min, k_max, nk) 

V = []
V.append(np.zeros(nk))
V.append(np.zeros(nk))
k_dec = []
k_dec.append(np.zeros(nk))

kj_dec = []
kj_dec.append(np.zeros(nk))


M = [0]*nk
Chi = [0]*nk
for i in range(nk):
    M[i] = [0] * nk
    Chi[i] = [0] * nk
    
for i in range(nk):
    for j in range(nk):
        M[i][j] = m(k[i], k[j])

for i in range(nk):
    for j in range(nk):
        Chi[i][j] = M[i][j] + beta* V[-1][j]

V_new = [0]*nk
k_dec_new = [0]*nk
for i in range(nk):
    V_new[i]=max(Chi[i])

V.append(V_new)

while (np.linalg.norm(np.asarray(V[-1]) - np.asarray(V[-2])) > 0.01): 
    for i in range(nk):
        j=0
        while j<nk:
            Chi[i][j] =M[i][j] + beta* V[-1][j]
            if j>0:
                if Chi[i][j]<Chi[i][j-1]:
                    j=nk
            j = j+1
            
    V_new = [0]*nk
    k_dec_new = [0]*nk
    for i in range(nk):
        V_new[i]=max(Chi[i])

    V.append(V_new)

T = qe.toc()

plt.plot(k, V[-1], label='V(k)', color="green")    
plt.xlabel('Capital (k)')
plt.ylabel('Value (Vk)')
plt.show()

print("Time elapse in seconds: " + str(T))
print("Number of iterations: " + str(j))


#%%(2.d)----------------------------------------------------------------#

qe.tic()

Cond = True
k_min = 0.01
k_max = 2
nk = 100

k = np.linspace(k_min, k_max, nk) 

V = []
V.append(np.zeros(nk)) 
V.append(np.zeros(nk))
k_dec = []
k_dec.append(np.zeros(nk))

kj_dec = []
kj_dec.append(np.zeros(nk))

M = [0]*nk
Chi = [0]*nk
for i in range(nk):
    M[i] = [0] * nk
    Chi[i] = [0] * nk
    
for i in range(nk):
    for j in range(nk):
        M[i][j] = m(k[i], k[j])

for i in range(nk):
    for j in range(nk):
        Chi[i][j] = M[i][j] + beta* V[-1][j]

V_new = [0]*nk
k_dec_new = [0]*nk
for i in range(nk):
    V_new[i]=max(Chi[i])
    k_dec_new[i] = Chi[i].index(max(Chi[i]))

V.append(V_new)
k_dec.append(k_dec_new)

while (np.linalg.norm(np.asarray(V[-1]) - np.asarray(V[-2])) > 0.01): 
    for i in range(nk):
        for j in range(k_dec_new[i]-1,k_dec_new[i]+1):
            Chi[i][j] =M[i][j] + beta* V[-1][j]

    V_new = [0]*nk
    k_dec_new = [0]*nk
    for i in range(nk):
        V_new[i]=max(Chi[i])
        k_dec_new[i] = Chi[i].index(max(Chi[i]))

    V.append(V_new)
    k_dec.append(k_dec_new)

T = qe.toc()

plt.plot(k, V[-1], label='V(k)', color="cyan")    
plt.xlabel('Capital (k)')
plt.ylabel('Value (Vk)')
plt.show()

print("Time elapse in seconds: " + str(T))
print("Number of iterations: " + str(j))


#%%(2.e)----------------------------------------------------------------#

qe.tic()
k_min = 0.01
k_max = 2
nk = 100

k = np.linspace(k_min, k_max, nk) 

V = []
V.append(np.zeros(nk)) 
V.append(np.zeros(nk))
k_dec = []
k_dec.append(np.zeros(nk))

kj_dec = []
kj_dec.append(np.zeros(nk))

M = [0]*nk
Chi = [0]*nk
for i in range(nk):
    M[i] = [0] * nk
    Chi[i] = [0] * nk
    
for i in range(nk):
    for j in range(nk):
        M[i][j] = m(k[i], k[j])

for i in range(nk):
    for j in range(nk):
        Chi[i][j] = M[i][j] + beta* V[-1][j]

V_new = [0]*nk
k_dec_new = [0]*nk
for i in range(nk):
    V_new[i]=max(Chi[i])

V.append(V_new)

while (np.linalg.norm(np.asarray(V[-1]) - np.asarray(V[-2])) > 0.01): 
    for i in range(nk):
        j=0
        while j<nk:
            if k[j]>=k_dec[-1][j]:
                Chi[i][j] =M[i][j] + beta* V[-1][j]
                if j>0:
                    if Chi[i][j]<Chi[i][j-1]:
                        j=nk
                j = j+1
    
    V_new = [0]*nk
    k_dec_new = [0]*nk
    for i in range(nk):
        V_new[i]=max(Chi[i])

    V.append(V_new)

T = qe.toc()

plt.plot(k, V[-1], label='V(k)', color="magenta")    
plt.xlabel('Capital (k)')
plt.ylabel('Value (Vk)')
plt.show()

print("Time elapse in seconds: " + str(T))
print("Number of iterations: " + str(j))


