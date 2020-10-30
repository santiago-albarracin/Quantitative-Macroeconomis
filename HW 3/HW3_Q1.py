# QUANTITATIVE MACROECONOMICS: 
# HWK 3 - VALUE FUNCTION ITERATION
    # Santiago Albarracin
    # Victor Caballero
    # Lorenzo de Carolis
    
import numpy as np
import matplotlib.pyplot as plt
import quantecon as qe


# =============================================================================
# EXERCISE 1: VFI with inelastic labour supply
# =============================================================================

# parameters 
theta = 0.679
beta = 0.988
delta = 0.013
h = 1

kSS = ((1/beta-1+delta)/(1-theta))**(-1/theta)

print("The steady state is "+str(kSS))

#%%(1.a)----------------------------------------------------------------#
# Brute force Value Function Iteration 


# step 1
qe.tic()
nk = 200
k = np.linspace(1, 2*kSS, nk)   # nk evenly spaced points

# step 2
M = np.empty([nk,nk])
for i,ki in enumerate(k):
        for j, kj in enumerate(k):
            if kj <= ki**(1- theta) + (1-delta)*ki:   # using only feasible values (c>0)
                M[i,j]  = np.log(ki**(1- theta) + (1-delta)*ki - kj)
            else:
                M[i,j] = -1000000

# step 3
epsilon = 0.01 
convergence = False
s = 0
S = 400
X = np.empty([nk, nk])
g = np.ones([nk, 1])
Vi = np.zeros([nk,1])
Vj = np.empty([nk,1])

# step 4: Value function iteration:
while (s<S) and (convergence == False):
    for i in range(nk):                   #step 4.1
        for j in range(nk):
            X[i,j] = M[i,j] + beta*Vi[j]
    for i in range(nk):                   #step 4.2
        Vj[i] = np.max(X[i,:])
        g[i] = np.argmax(X[i,:])          #step 4.4
    if np.max(Vj - Vi) >= epsilon:        #step 4.3
        Vi = np.copy(Vj)
        Vj = np.empty([nk,1])
        s += 1
    else:
        convergence = True

# step 5:
gk = np.empty(nk)
gc = np.empty(nk)
for i in range(nk):
    gk[i] = k[int(g[i])]
    gc[i] = k[i]**(1- theta) + (1-delta)*k[i] - gk[i]


T = qe.toc()

# Plotting:
plt.plot(k, Vj, label='V(k)', color="red")    
plt.xlabel('Capital (k)')
plt.ylabel('Value (Vk)')
plt.show()

plt.plot(k,gk, label = "Policy function", color="red")
plt.plot(k,k, label = "45 degree", color="black")
plt.axvline(kSS,label = "kSS", color='black', linestyle='dashed')
plt.xlabel('kt')
plt.ylabel('kt+1')
plt.legend()
plt.show()

print("Time elapse in seconds: " + str(T))
print("Number of iterations: " + str(s))


#%%(1.b)----------------------------------------------------------------#
# Monotonicity Value Function Iteration 


qe.tic()
nk = 200
k = np.linspace(1, 2*kSS, nk)   # nk evenly spaced points

# step 2
M = np.empty([nk,nk])
for i,ki in enumerate(k):
        for j, kj in enumerate(k):
            if kj <= ki**(1- theta) + (1-delta)*ki:   # using only feasible values (c>0)
                M[i,j]  = np.log(ki**(1- theta) + (1-delta)*ki - kj)
            else:
                M[i,j] = -1000000

# step 3
epsilon = 0.01 
convergence = False
s = 0
S = 500
X = np.empty([nk, nk])
g = np.ones([nk, 1])
Vi = np.zeros([nk,1])
Vj = np.empty([nk,1])

# step 4: Value function iteration:
while (s<S) and (convergence == False):
    for i in range(nk):                   #step 4.1
        for j in range(nk):
            gLB = int(g[i])    
            if j+gLB < nk:    
                X[i,j+gLB] = M[i,j+gLB] + beta*Vi[j+gLB]
            else:
                continue
    for i in range(nk):                   #step 4.2
        Vj[i] = np.max(X[i,:])
        g[i] = np.argmax(X[i,:])          #step 4.4
    if np.max(Vj - Vi) >= epsilon:        #step 4.3
        Vi = np.copy(Vj)
        Vj = np.empty([nk,1])
        s += 1
    else:
        convergence = True

# step 5:
gk = np.empty(nk)
gc = np.empty(nk)
for i in range(nk):
    gk[i] = k[int(g[i])]
    gc[i] = k[i]**(1- theta) + (1-delta)*k[i] - gk[i]


T = qe.toc()

# Plotting:
plt.plot(k, Vj, label='V(k)', color="blue")       
plt.xlabel('Capital (k)')
plt.ylabel('Value (Vk)')
plt.show()

plt.plot(k,gk, label = "Policy function", color="blue")
plt.plot(k,k, label = "45 degree", color="black")
plt.axvline(kSS,label = "kSS", color='black', linestyle='dashed')
plt.xlabel('kt')
plt.ylabel('kt+1')
plt.legend()
plt.show()


print("Time elapse in seconds: " + str(T))
print("Number of iterations: " + str(s))


#%%(1.c)----------------------------------------------------------------#
# concavity Value Function Iteration 


qe.tic()
nk = 200
k = np.linspace(1, 2*kSS, nk)   # nk evenly spaced points

# step 2
M = np.empty([nk,nk])
for i,ki in enumerate(k):
        for j, kj in enumerate(k):
            if kj <= ki**(1- theta) + (1-delta)*ki:   # using only feasible values (c>0)
                M[i,j]  = np.log(ki**(1- theta) + (1-delta)*ki - kj)
            else:
                M[i,j] = -1000000

# step 3
epsilon = 0.01 
convergence = False
s = 0
S = 500
X = np.empty([nk, nk])
g = np.ones([nk, 1])
Vi = np.zeros([nk,1])
Vj = np.empty([nk,1])

# step 4: Value function iteration:
while (s<S) and (convergence == False):
    for i in range(nk):                   #step 4.1
        for j in range(nk):
            X[i,j] = M[i,j] + beta*Vi[j]
            if X[i,j] < X[i,j-1]:  
                break
    for i in range(nk):                   #step 4.2
        Vj[i] = np.max(X[i,:])
        g[i] = np.argmax(X[i,:])          #step 4.4
    if np.max(Vj - Vi) >= epsilon:        #step 4.3
        Vi = np.copy(Vj)
        Vj = np.empty([nk,1])
        s += 1
    else:
        convergence = True

# step 5:
gk = np.empty(nk)
gc = np.empty(nk)
for i in range(nk):
    gk[i] = k[int(g[i])]
    gc[i] = k[i]**(1- theta) + (1-delta)*k[i] - gk[i]


T = qe.toc()

# Plotting:
plt.plot(k, Vj, label='V(k)', color="green")       
plt.xlabel('Capital (k)')
plt.ylabel('Value (Vk)')
plt.show()

plt.plot(k,gk, label = "Policy function", color="green")
plt.plot(k,k, label = "45 degree", color="black")
plt.axvline(kSS,label = "kSS", color='black', linestyle='dashed')
plt.xlabel('kt')
plt.ylabel('kt+1')
plt.legend()
plt.show()


print("Time elapse in seconds: " + str(T))
print("Number of iterations: " + str(s))


#%%(1.d)----------------------------------------------------------------#
# Local search Value Function Iteration 


qe.tic()
nk = 200
k = np.linspace(1, 2*kSS, nk)   # nk evenly spaced points

# step 2
M = np.empty([nk,nk])
for i,ki in enumerate(k):
        for j, kj in enumerate(k):
            if kj <= ki**(1- theta) + (1-delta)*ki:   # using only feasible values (c>0)
                M[i,j]  = np.log(ki**(1- theta) + (1-delta)*ki - kj)
            else:
                M[i,j] = -1000000

# step 3
epsilon = 0.01 
convergence = False
s = 0
S = 500
X = np.empty([nk, nk])
g = np.ones([nk, 1])
Vi = np.zeros([nk,1])
Vj = np.empty([nk,1])

# step 4: Value function iteration:
while (s<S) and (convergence == False):
    for i in range(nk):                   #step 4.1
        for j in range(nk):
            if (j >= g[i]) and (j <= g[i] + 5) :
                X[i,j] = M[i,j] + beta*Vi[j]
    for i in range(nk):                   #step 4.2
        Vj[i] = np.max(X[i,:])
        g[i] = np.argmax(X[i,:])          #step 4.4
    if np.max(Vj - Vi) >= epsilon:        #step 4.3
        Vi = np.copy(Vj)
        Vj = np.empty([nk,1])
        s += 1
    else:
        convergence = True

# step 5:
gk = np.empty(nk)
gc = np.empty(nk)
for i in range(nk):
    gk[i] = k[int(g[i])]
    gc[i] = k[i]**(1- theta) + (1-delta)*k[i] - gk[i]


T = qe.toc()

# Plotting:
plt.plot(k, Vj, label='V(k)', color="cyan")       
plt.xlabel('Capital (k)')
plt.ylabel('Value (Vk)')
plt.show()

plt.plot(k,gk, label = "Policy function", color="cyan")
plt.plot(k,k, label = "45 degree", color="black")
plt.axvline(kSS,label = "kSS", color='black', linestyle='dashed')
plt.xlabel('kt')
plt.ylabel('kt+1')
plt.legend()
plt.show()


print("Time elapse in seconds: " + str(T))
print("Number of iterations: " + str(s))


#%%(1.e)----------------------------------------------------------------#
# Concavity and monotonicity Value Function Iteration 


qe.tic()
nk = 200
k = np.linspace(1, 2*kSS, nk)   # nk evenly spaced points

# step 2
M = np.empty([nk,nk])
for i,ki in enumerate(k):
        for j, kj in enumerate(k):
            if kj <= ki**(1- theta) + (1-delta)*ki:   # using only feasible values (c>0)
                M[i,j]  = np.log(ki**(1- theta) + (1-delta)*ki - kj)
            else:
                M[i,j] = -1000000

# step 3
epsilon = 0.01 
convergence = False
s = 0
S = 500
X = np.empty([nk, nk])
g = np.ones([nk, 1])
Vi = np.zeros([nk,1])
Vj = np.empty([nk,1])

# step 4: Value function iteration:
while (s<S) and (convergence == False):
    for i in range(nk):                   #step 4.1
        for j in range(nk):
            gLB = int(g[i])   
            if j+gLB < nk:     
                X[i,j+gLB] = M[i,j+gLB] + beta*Vi[j+gLB]
                if X[i,j] < X[i,j-1]:
                    break
            else:
                continue
    for i in range(nk):                   #step 4.2
        Vj[i] = np.max(X[i,:])
        g[i] = np.argmax(X[i,:])          #step 4.4
    if np.max(Vj - Vi) >= epsilon:        #step 4.3
        Vi = np.copy(Vj)
        Vj = np.empty([nk,1])
        s += 1
    else:
        convergence = True

# step 5:
gk = np.empty(nk)
gc = np.empty(nk)
for i in range(nk):
    gk[i] = k[int(g[i])]
    gc[i] = k[i]**(1- theta) + (1-delta)*k[i] - gk[i]


T = qe.toc()

# Plotting:
plt.plot(k, Vj, label='V(k)', color="magenta")       
plt.xlabel('Capital (k)')
plt.ylabel('Value (Vk)')
plt.show()

plt.plot(k,gk, label = "Policy function", color="magenta")
plt.plot(k,k, label = "45 degree", color="black")
plt.axvline(kSS,label = "kSS", color='black', linestyle='dashed')
plt.xlabel('kt')
plt.ylabel('kt+1')
plt.legend()
plt.show()


print("Time elapse in seconds: " + str(T))
print("Number of iterations: " + str(s))
