# Santiago Albarracin
#   Quantitative Macroeconomics
#   Project

# Import packages
import numpy as np
import matplotlib.pyplot as plt
import random

###############################################################################
#First part
PI = [ [7/8, 1/8], 
       [1/8, 7/8] ]

# unemployment bad around 20% 
# unemployment good around 10% 
PI_v = [0.31166667,0.055,0.01125,0.02267593,0.56333333,0.82,0.11375,0.10232407,0.11375,0.00682292,0.575,0.07,0.01125,0.11817708,0.3,0.805]
PI_z = np.reshape(PI_v,[4,4], order = 'F')

# Parameters
beta=0.95
delta=0.0025
z=[1.01, 0.99]
alfa=0.36
L=[0.96, 0.9]

# Initial values
def v1good(k,K):
    return np.log(alfa*z[0]*(K/L[0])**(alfa-1)*k + (1-alfa)*z[0]*(K/L[0])**(alfa)-delta*k)/(1-beta)  
def v1bad(k,K):
    return np.log(alfa*z[1]*(K/L[1])**(alfa-1)*k + (1-alfa)*z[1]*(K/L[1])**(alfa)-delta*k)/(1-beta)  
def v0good(k,K):
    return np.log(alfa*z[0]*(K/L[0])**(alfa-1)*k -delta*k)/(1-beta)  
def v0bad(k,K):
    return np.log(alfa*z[1]*(K/L[1])**(alfa-1)*k -delta*k)/(1-beta) 

# Grid

k_grid = np.linspace(0.1,60,120)

x=1
i = 0
K_grid=[]
while x==1:
    K_grid.append(round(16+i*0.15, 2))
    i = i+1
    if 16+i*0.15 >= 45:
        x = 0
        
#Value function
V1good = [0]*len(k_grid)
V1bad = [0]*len(k_grid)
V0good = [0]*len(k_grid)
V0bad = [0]*len(k_grid)
V1goodt = [0]*len(k_grid)
V1badt = [0]*len(k_grid)
V0goodt = [0]*len(k_grid)
V0badt = [0]*len(k_grid)

for i in range(len(k_grid)):
    V1good[i] = [0]*len(K_grid)
    V1bad[i] = [0]*len(K_grid)
    V0good[i] = [0]*len(K_grid)
    V0bad[i] = [0]*len(K_grid)
    V1goodt[i] = [0]*len(K_grid)
    V1badt[i] = [0]*len(K_grid)
    V0goodt[i] = [0]*len(K_grid)
    V0badt[i] = [0]*len(K_grid)     

for i in range(len(k_grid)):
    for j in range(len(K_grid)):
        V1good[i][j]= v1good(k_grid[i],K_grid[j])
        V1bad[i][j]= v1bad(k_grid[i],K_grid[j])
        V0good[i][j]= v0good(k_grid[i],K_grid[j])
        V0bad[i][j]= v0bad(k_grid[i],K_grid[j])
        
V1good2 = [0]*len(k_grid)
V1bad2 = [0]*len(k_grid)
V0good2 = [0]*len(k_grid)
V0bad2 = [0]*len(k_grid)
V1goodt2 = [0]*len(k_grid)
V1badt2 = [0]*len(k_grid)
V0goodt2 = [0]*len(k_grid)
V0badt2 = [0]*len(k_grid)

for i in range(len(k_grid)):
    V1good2[i] = [0]*len(K_grid)
    V1bad2[i] = [0]*len(K_grid)
    V0good2[i] = [0]*len(K_grid)
    V0bad2[i] = [0]*len(K_grid)
    V1goodt2[i] = [0]*len(K_grid)
    V1badt2[i] = [0]*len(K_grid)
    V0goodt2[i] = [0]*len(K_grid)
    V0badt2[i] = [0]*len(K_grid)     

for i in range(len(k_grid)):
    for j in range(len(K_grid)):
        V1good2[i][j]= v1good(k_grid[i],K_grid[j])
        V1bad2[i][j]= v1bad(k_grid[i],K_grid[j])
        V0good2[i][j]= v0good(k_grid[i],K_grid[j])
        V0bad2[i][j]= v0bad(k_grid[i],K_grid[j])

## Law of motion  
b0good=0
b1good=1
b0bad=0
b1bad=1

iter_b = 1
stop = 0
while iter_b<15 and stop ==0:
    def H(K,zi):
        return np.exp( (b0good+b1good*np.log(K))*zi+ (b0bad+b1bad*np.log(K))*(1-zi))
    
#Consumption problem
    def c(i,I,e,g):
        vec = np.dot(alfa*z[g]*(K_grid[I]/L[g])**(alfa-1),k_grid)+(1-alfa)*z[g]*(K_grid[I]/L[g])**(alfa)*e+(1-delta)*k_grid[i]- k_grid
        for i in range(len(vec)):
            if vec[i]<0:
                vec[i]=0
        return(vec)    
    
    def column(matrix, i):
        return [row[i] for row in matrix]

    iter = 1
    STOP = 0
    while iter < 15 and STOP == 0:
        i=0
        for i in range(len(k_grid)):
            I = 0
            for I in range(len(K_grid)):
                values = abs(K_grid-H(K_grid[I],1))
                Ip= np.argmin(values)
                V0goodt[i][I] = max(np.log(c(i,I,0,0)) + beta*(np.matmul(PI_z[0],[column(V0good, Ip), column(V1good, Ip), column(V0bad,Ip), column(V1bad,Ip)])))
                V1goodt[i][I] = max(np.log(c(i,I,1,0)) + beta*(np.matmul(PI_z[1],[column(V0good, Ip), column(V1good, Ip), column(V0bad,Ip), column(V1bad,Ip)])))
                 
                values = abs(K_grid-H(K_grid[I],0))
                Ip= np.argmin(values)
                V0badt[i][I] = max(np.log(c(i,I,0,1)) + beta*(np.matmul(PI_z[2],[column(V0good, Ip), column(V1good, Ip), column(V0bad,Ip), column(V1bad,Ip)])))
                V1badt[i][I] = np.nanmax(np.log(c(i,I,1,1)) + beta*(np.matmul(PI_z[3],[column(V0good, Ip), column(V1good, Ip), column(V0bad,Ip), column(V1bad,Ip)])))
        
        temp=np.array([abs(np.array(V0goodt)-np.array(V0good)),abs(np.array(V1goodt)-np.array(V1good)),abs(np.array(V0badt)-np.array(V0bad)),abs(np.array(V1badt)-np.array(V1bad))])
        dev = temp.max()        
              
        for i in range(len(k_grid)):
            for j in range(len(K_grid)):
                V0good[i][j]=V0goodt[i][j]
                V1good[i][j]=V1goodt[i][j]
                V0bad[i][j]=V0badt[i][j]
                V1bad[i][j]=V1badt[i][j]      
                
        for i in range(len(k_grid)):
            for j in range(len(K_grid)):
                V0good2[i][j]=V0goodt2[i][j]
                V1good2[i][j]=V1goodt2[i][j]
                V0bad2[i][j]=V0badt2[i][j]
                V1bad2[i][j]=V1badt2[i][j]  
        
        if dev<0.1:
            STOP = 1
        iter = iter + 1

    #Policy function
    pf = [0]*len(k_grid)
    for i in range(len(k_grid)):
        pf[i]=[0]*len(K_grid)
        for I in range(len(K_grid)):
            pf[i][I]=[0]*2
            for x in range(2):
                pf[i][I][x] = [0]*2                
            
    for i in range(len(k_grid)):
        for I in range(len(K_grid)):
            values = abs(K_grid-H(K_grid[I],1))
            Ip= np.argmin(values)           
            pf[i][I][1][0] = np.argmax(np.log(c(i,I,0,0)) + beta*(np.matmul(PI_z[0],[column(V0good, Ip), column(V1good, Ip), column(V0bad,Ip), column(V1bad,Ip)])))
            pf[i][I][0][0] = np.argmax(np.log(c(i,I,1,0)) + beta*(np.matmul(PI_z[1],[column(V0good, Ip), column(V1good, Ip), column(V0bad,Ip), column(V1bad,Ip)])))
            
            values = abs(K_grid-H(K_grid[I],0))
            Ip= np.argmin(values)
            pf[i][I][1][1] = np.argmax(np.log(c(i,I,0,1)) + beta*(np.matmul(PI_z[2],[column(V0good, Ip), column(V1good, Ip), column(V0bad,Ip), column(V1bad,Ip)])))
            pf[i][I][0][1] = np.argmax(np.log(c(i,I,1,1)) + beta*(np.matmul(PI_z[3],[column(V0good, Ip), column(V1good, Ip), column(V0bad,Ip), column(V1bad,Ip)])))
    
    ar = np.array(pf)
      
    if iter_b==1:
        zt = [0]*50          
        zt[0]=0            
        for t in range(1,50):
            draw=random.randint(0,100)/100
            zt[t]= int(draw>=PI[zt[t-1]][0])
                
        i_zgood = []
        i_zbad = []      
        for t in range(len(zt)):
            if zt[t]==0:
                i_zgood.append(t)  
            else:
                i_zbad.append(t)           
        
        # Assets and employment
        state = [40]*1000
        for i in range(len(state)):
            state[i] = [0]*2
            for j in range(2):
                state[i][j] = [0]*50
        
        for i in range(len(state)):
            if i <880: 
                for j in range(2):
                    state[i][j][0] = (1-j)*40+j*0
            else:
                for j in range(2):
                    state[i][j][0] = (1-j)*40+j*1
        
        state = np.array(state)
        
        # Aggregate capital    
        K_ind = [0]*50
        K_ind[0]=3
        for t in range(1,50):
            temp_sum = 0
            for n in range(1000):
                state[n][0][t]=pf[state[n][0][t-1]][K_ind[t-1]][state[n][1][t-1]][zt[t-1]]
                if zt[t-1]==0 and zt[t]==0 and state[n][1][t-1]==0: 
                    state[n][1][t]= 1-((random.randint(0,100)/100)>=PI_z[1][0]/PI[0][0])    
                if zt[t-1]==0 and zt[t]==1 and state[n][1][t-1]==0: 
                    state[n][1][t]= 1-((random.randint(0,100)/100)>=PI_z[1][2]/PI[0][1])
                if zt[t-1]==1 and zt[t]==0 and state[n][1][t-1]==0: 
                    state[n][1][t]= 1-((random.randint(0,100)/100)>=PI_z[3][0]/PI[1][0])
                if zt[t-1]==1 and zt[t]==1 and state[n][1][t-1]==0: 
                    state[n][1][t]= 1-((random.randint(0,100)/100)>=PI_z[3][2]/PI[1][1])
                if zt[t-1]==0 and zt[t]==0 and state[n][1][t-1]==1: 
                    state[n][1][t]= 1-((random.randint(0,100)/100)>=PI_z[0][0]/PI[0][0])
                if zt[t-1]==0 and zt[t]==1 and state[n][1][t-1]==1: 
                    state[n][1][t]= 1-((random.randint(0,100)/100)>=PI_z[0][2]/PI[0][1])
                if zt[t-1]==1 and zt[t]==0 and state[n][1][t-1]==1: 
                    state[n][1][t]= 1-((random.randint(0,100)/100)>=PI_z[2][0]/PI[1][0])
                if zt[t-1]==1 and zt[t]==1 and state[n][1][t-1]==1: 
                    state[n][1][t]= 1-((random.randint(0,100)/100)>=PI_z[2][2]/PI[1][1])
                temp_sum = temp_sum + state[n][0][t]
            K_ind[t]=np.argmin(abs(k_grid[int(round(temp_sum/1000))]-K_grid))
    
    else:
        for t in range(1,50):
            temp_sum = 0
            for n in range(1000):
                state[n][0][t]=pf[state[n][0][t-1]][K_ind[t-1]][state[n][1][t-1]][zt[t-1]]
                temp_sum = temp_sum + state[n][0][t]
            K_ind[t]=np.argmin(abs(k_grid[int(round(temp_sum/1000))]-K_grid))   
    
    #state    
    state_temp = [0]*1000
    for i in range(1000):
        state_temp[i]=[0]*2       
    for n in range(1000):
        state_temp[n][0]=state[n][0][49]
        state_temp[n][1]=state[n][1][49]       
    state_old = state
        
    # Regressions 
    Yg = [0]*len(i_zgood)
    Xg = [0]*len(i_zgood)
    for i in range(len(i_zgood)):
        Xg[i] = [0]*2
        
    for i in range(20,len(i_zgood)):
        Yg[i] = np.log(K_grid[K_ind[i_zgood[i]-1]])
        for j in range(2):
            Xg[i][0] = 1
            Xg[i][1] = np.log(K_grid[K_ind[i_zgood[i]-2]])
    
    Bg = np.linalg.lstsq(Xg,Yg)
    b0goodp=Bg[0][0]
    b1goodp=Bg[0][1]
    
    Yb = [0]*len(i_zbad)
    Xb = [0]*len(i_zbad)
    for i in range(len(i_zbad)):
        Xb[i] = [0]*2
        
    for i in range(len(i_zbad)):
        Yb[i] = np.log(K_grid[K_ind[i_zbad[i]-1]])
        for j in range(2):
            Xb[i][0] = 1
            Xb[i][1] = np.log(K_grid[K_ind[i_zbad[i]-2]])
   
    Bb = np.linalg.lstsq(Xb,Yb)
    b0badp=Bb[0][0]
    b1badp=Bb[0][1]
    
    dev_b=max([abs(b0good-b0goodp), abs(b1good-b1goodp), abs(b0bad-b0badp), abs(b1bad-b1badp)])
    
    if dev_b<=0.01:
        stop = 1
    
    b0good=0.1*b0goodp+0.9*b0good
    b1good=0.1*b1goodp+0.9*b1good
    b0bad=0.1*b0badp+0.9*b0bad
    b1bad=0.1*b1badp+0.9*b1bad    
    iter_b = iter_b + 1

#Plot
plt.plot(k_grid, k_grid[ar[:,22,1,0]], label='Unemployment, good')
plt.plot(k_grid, k_grid[ar[:,24,1,1]], label='Unemployment, bad')
plt.plot(k_grid, k_grid, label='45`')
plt.xlim([0,5]) 
plt.ylim([0,5])
plt.legend() 
plt.show()

##############################################################################
#Second part

PI = [ [7/8, 1/8], 
       [1/8, 7/8] ] 

# unemployment bad around 10% 
# unemployment good around 5% 

PI_v = [0.27166667,0.025,0.01125,0.00267593,0.60333333,0.85,0.11375,0.12232407,0.11375,0.00682292,0.425,0.04,0.01125,0.11817708,0.45,0.835]
PI_z = np.reshape(PI_v,[4,4], order = 'F')

# Parameters
beta=0.95
delta=0.0025
z=[1.01, 0.99]
alfa=0.36
L=[0.96, 0.9]

# Initial values
def v1good2(k,K):
    return np.log(alfa*z[0]*(K/L[0])**(alfa-1)*k + (1-alfa)*z[0]*(K/L[0])**(alfa)-delta*k)/(1-beta) 
def v1bad2(k,K):
    return np.log(alfa*z[1]*(K/L[1])**(alfa-1)*k + (1-alfa)*z[1]*(K/L[1])**(alfa)-delta*k)/(1-beta) 
def v0good2(k,K):
    return np.log(alfa*z[0]*(K/L[0])**(alfa-1)*k -delta*k)/(1-beta) 
def v0bad2(k,K):
    return np.log(alfa*z[1]*(K/L[1])**(alfa-1)*k -delta*k)/(1-beta)

# Grid

k_grid = np.linspace(0.1,60,120)


x=1
i = 0
K_grid=[]
while x==1:
    K_grid.append(round(16+i*0.15, 2))
    i = i+1
    if 16+i*0.15 >= 45:
        x = 0
        
#Value function
V1good2 = [0]*len(k_grid)
V1bad2 = [0]*len(k_grid)
V0good2 = [0]*len(k_grid)
V0bad2 = [0]*len(k_grid)
V1goodt2 = [0]*len(k_grid)
V1badt2 = [0]*len(k_grid)
V0goodt2 = [0]*len(k_grid)
V0badt2 = [0]*len(k_grid)

for i in range(len(k_grid)):
    V1good2[i] = [0]*len(K_grid)
    V1bad2[i] = [0]*len(K_grid)
    V0good2[i] = [0]*len(K_grid)
    V0bad2[i] = [0]*len(K_grid)
    V1goodt2[i] = [0]*len(K_grid)
    V1badt2[i] = [0]*len(K_grid)
    V0goodt2[i] = [0]*len(K_grid)
    V0badt2[i] = [0]*len(K_grid)     

for i in range(len(k_grid)):
    for j in range(len(K_grid)):
        V1good2[i][j]= v1good2(k_grid[i],K_grid[j])
        V1bad2[i][j]= v1bad2(k_grid[i],K_grid[j])
        V0good2[i][j]= v0good2(k_grid[i],K_grid[j])
        V0bad2[i][j]= v0bad2(k_grid[i],K_grid[j])

## Law of motion  
b0good=0
b1good=1
b0bad=0
b1bad=1

iter_b = 1
stop = 0
while iter_b<15 and stop ==0:
    print(iter_b)
    def H(K,zi):
        return np.exp( (b0good+b1good*np.log(K))*zi+ (b0bad+b1bad*np.log(K))*(1-zi))
    
#Consumption problem
    def c(i,I,e,g):
        vec = np.dot(alfa*z[g]*(K_grid[I]/L[g])**(alfa-1),k_grid)+(1-alfa)*z[g]*(K_grid[I]/L[g])**(alfa)*e+(1-delta)*k_grid[i]- k_grid
        for i in range(len(vec)):
            if vec[i]<0:
                vec[i]=0
        return(vec)
    
    
    def column(matrix, i):
        return [row[i] for row in matrix]
    
    iter = 1
    STOP = 0
    while iter < 15 and STOP == 0:
        i=0
        for i in range(len(k_grid)):
            I = 0
            for I in range(len(K_grid)):
                values = abs(K_grid-H(K_grid[I],1))
                Ip= np.argmin(values)
                V0goodt2[i][I] = max(np.log(c(i,I,0,0)) + beta*(np.matmul(PI_z[0],[column(V0good2, Ip), column(V1good2, Ip), column(V0bad2,Ip), column(V1bad2,Ip)])))
                V1goodt2[i][I] = max(np.log(c(i,I,1,0)) + beta*(np.matmul(PI_z[1],[column(V0good2, Ip), column(V1good2, Ip), column(V0bad2,Ip), column(V1bad2,Ip)])))
                 
                values = abs(K_grid-H(K_grid[I],0))
                Ip= np.argmin(values)
                V0badt2[i][I] = max(np.log(c(i,I,0,1)) + beta*(np.matmul(PI_z[2],[column(V0good2, Ip), column(V1good2, Ip), column(V0bad2,Ip), column(V1bad2,Ip)])))
                V1badt2[i][I] = np.nanmax(np.log(c(i,I,1,1)) + beta*(np.matmul(PI_z[3],[column(V0good2, Ip), column(V1good2, Ip), column(V0bad2,Ip), column(V1bad2,Ip)])))
        
        temp=np.array([abs(np.array(V0goodt2)-np.array(V0good2)),abs(np.array(V1goodt2)-np.array(V1good2)),abs(np.array(V0badt2)-np.array(V0bad2)),abs(np.array(V1badt2)-np.array(V1bad2))])
        dev = temp.max()        
        
        for i in range(len(k_grid)):
            for j in range(len(K_grid)):
                V0good2[i][j]=V0goodt2[i][j]
                V1good2[i][j]=V1goodt2[i][j]
                V0bad2[i][j]=V0badt2[i][j]
                V1bad2[i][j]=V1badt2[i][j]  
        
        if dev<0.1:
            STOP = 1

        iter = iter + 1
      
    #Policy function
    pf2 = [0]*len(k_grid)
    for i in range(len(k_grid)):
        pf2[i]=[0]*len(K_grid)
        for I in range(len(K_grid)):
            pf2[i][I]=[0]*2
            for x in range(2):
                pf2[i][I][x] = [0]*2
                
                
    for i in range(len(k_grid)):
        for I in range(len(K_grid)):
            values = abs(K_grid-H(K_grid[I],1))
            Ip= np.argmin(values)           
            pf2[i][I][1][0] = np.argmax(np.log(c(i,I,0,0)) + beta*(np.matmul(PI_z[0],[column(V0good2, Ip), column(V1good2, Ip), column(V0bad2,Ip), column(V1bad2,Ip)])))
            pf2[i][I][0][0] = np.argmax(np.log(c(i,I,1,0)) + beta*(np.matmul(PI_z[1],[column(V0good2, Ip), column(V1good2, Ip), column(V0bad2,Ip), column(V1bad2,Ip)])))
            
            values = abs(K_grid-H(K_grid[I],0))
            Ip= np.argmin(values)
            pf2[i][I][1][1] = np.argmax(np.log(c(i,I,0,1)) + beta*(np.matmul(PI_z[2],[column(V0good2, Ip), column(V1good2, Ip), column(V0bad2,Ip), column(V1bad2,Ip)])))
            pf2[i][I][0][1] = np.argmax(np.log(c(i,I,1,1)) + beta*(np.matmul(PI_z[3],[column(V0good2, Ip), column(V1good2, Ip), column(V0bad2,Ip), column(V1bad2,Ip)])))
    
    ar2 = np.array(pf2)  
    
    if iter_b==1:
        zt = [0]*50  
        
        zt[0]=0            
        for t in range(1,50):
            draw=random.randint(0,100)/100
            zt[t]= int(draw>=PI[zt[t-1]][0])
               
        i_zgood = []
        i_zbad = []
        
        for t in range(len(zt)):
            if zt[t]==0:
                i_zgood.append(t)
            else:
                i_zbad.append(t)           
                
        # Assets and employment    
        state = [40]*1000
        for i in range(len(state)):
            state[i] = [0]*2
            for j in range(2):
                state[i][j] = [0]*50
        
        for i in range(len(state)):
            state[i][0][0] = state_temp[i][0]
            state[i][1][0] = state_temp[i][1]
        
        state = np.array(state)
        
        # Aggregate capital      
        K_ind = [0]*50
        K_ind[0]=3
        for t in range(1,50):
            temp_sum = 0
            for n in range(1000):
                state[n][0][t]=pf[state[n][0][t-1]][K_ind[t-1]][state[n][1][t-1]][zt[t-1]]
                if zt[t-1]==0 and zt[t]==0 and state[n][1][t-1]==0: 
                    state[n][1][t]= 1-((random.randint(0,100)/100)>=PI_z[1][0]/PI[0][0])    
                if zt[t-1]==0 and zt[t]==1 and state[n][1][t-1]==0: 
                    state[n][1][t]= 1-((random.randint(0,100)/100)>=PI_z[1][2]/PI[0][1])
                if zt[t-1]==1 and zt[t]==0 and state[n][1][t-1]==0: 
                    state[n][1][t]= 1-((random.randint(0,100)/100)>=PI_z[3][0]/PI[1][0])
                if zt[t-1]==1 and zt[t]==1 and state[n][1][t-1]==0: 
                    state[n][1][t]= 1-((random.randint(0,100)/100)>=PI_z[3][2]/PI[1][1])
                if zt[t-1]==0 and zt[t]==0 and state[n][1][t-1]==1: 
                    state[n][1][t]= 1-((random.randint(0,100)/100)>=PI_z[0][0]/PI[0][0])
                if zt[t-1]==0 and zt[t]==1 and state[n][1][t-1]==1: 
                    state[n][1][t]= 1-((random.randint(0,100)/100)>=PI_z[0][2]/PI[0][1])
                if zt[t-1]==1 and zt[t]==0 and state[n][1][t-1]==1: 
                    state[n][1][t]= 1-((random.randint(0,100)/100)>=PI_z[2][0]/PI[1][0])
                if zt[t-1]==1 and zt[t]==1 and state[n][1][t-1]==1: 
                    state[n][1][t]= 1-((random.randint(0,100)/100)>=PI_z[2][2]/PI[1][1])
                temp_sum = temp_sum + state[n][0][t]
            K_ind[t]=np.argmin(abs(k_grid[int(round(temp_sum/1000))]-K_grid))
    
    else:
        for t in range(1,50):
            temp_sum = 0
            for n in range(1000):
                state[n][0][t]=pf[state[n][0][t-1]][K_ind[t-1]][state[n][1][t-1]][zt[t-1]]
                temp_sum = temp_sum + state[n][0][t]
            K_ind[t]=np.argmin(abs(k_grid[int(round(temp_sum/1000))]-K_grid))

    #Regressions
    Yg = [0]*len(i_zgood)
    Xg = [0]*len(i_zgood)
    for i in range(len(i_zgood)):
        Xg[i] = [0]*2
        
    for i in range(20,len(i_zgood)):
        Yg[i] = np.log(K_grid[K_ind[i_zgood[i]-1]])
        for j in range(2):
            Xg[i][0] = 1
            Xg[i][1] = np.log(K_grid[K_ind[i_zgood[i]-2]])
    
    Bg = np.linalg.lstsq(Xg,Yg)
    b0goodp=Bg[0][0]
    b1goodp=Bg[0][1]
    
    Yb = [0]*len(i_zbad)
    Xb = [0]*len(i_zbad)
    for i in range(len(i_zbad)):
        Xb[i] = [0]*2
        
    for i in range(20,len(i_zbad)):
        Yb[i] = np.log(K_grid[K_ind[i_zbad[i]-1]])
        for j in range(2):
            Xb[i][0] = 1
            Xb[i][1] = np.log(K_grid[K_ind[i_zbad[i]-2]])
    
    Bb = np.linalg.lstsq(Xb,Yb)
    b0badp=Bb[0][0]
    b1badp=Bb[0][1]

    
    dev_b=max([abs(b0good-b0goodp), abs(b1good-b1goodp), abs(b0bad-b0badp), abs(b1bad-b1badp)])
    
    if dev_b<=0.01:
        stop = 1
    
    b0good=0.1*b0goodp+0.9*b0good
    b1good=0.1*b1goodp+0.9*b1good
    b0bad=0.1*b0badp+0.9*b0bad
    b1bad=0.1*b1badp+0.9*b1bad
    
    iter_b = iter_b + 1

#Plot
plt.plot(k_grid, k_grid[ar2[:,22,1,0]], label='Unemployment, good')
plt.plot(k_grid, k_grid[ar2[:,24,1,1]], label='Unemployment, bad')
plt.plot(k_grid, k_grid, label='45`')
plt.xlim([0,5]) 
plt.ylim([0,5])
plt.legend() 
plt.show()

##############################################################################
# Comparison
state_join = [40]*1000
for i in range(len(state_join)):
    state_join[i] = [0]*2
    for j in range(2):
        state_join[i][j] = [0]*100

for t in range(100):
    for n in range(1000):
        if t<50:
            state_join[n][0][t] = state_old[n][0][t]
            state_join[n][1][t] = state_old[n][1][t]
        else:
            state_join[n][0][t] = state[n][0][t-50]
            state_join[n][1][t] = state[n][1][t-50]
 
#Plot            
unemp2 = [0]*100
for t in range(100):
    for n in range(1000):   
        unemp2[t] = unemp2[t] + state_join[n][1][t]
plt.plot(range(100), unemp2)
plt.show()

state_join = np.array(state_join)

#Plot
plt.hist(k_grid[state_join[:,0,51]], 40, label='First part (50)')
plt.legend() 
plt.show()

#Plot
plt.hist(k_grid[state_join[:,0,90]], 40, label='Second part (100)')
plt.legend() 
plt.show()