#   Quantitative Macroeconomics
#   HW 4

# Import packages
import numpy as np
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------------------------------------------------------
# I) A simple wealth model

class HW4:
    def __init__(self, rho, r, sigma, c_bar, gamma, sigma_y, y):
        self.rho = rho
        self.r = r
        self.sigma = sigma
        self.c_bar = c_bar
        self.gamma = gamma
        self.sigma_y = sigma_y
        self.beta = 1/(1+rho)
        self.Y = [y-sigma_y, y+sigma_y]
        self.pi = [[(1+gamma)/2, (1-gamma)/2], [(1-gamma)/2, (1+gamma)/2]]
        self.A_min = round(-min(self.Y)*((1+r)/r)/1.5)
        self.p = abs(self.A_min)+20+1
        self.A = np.linspace(self.A_min, 20, self.p)
    
    def u(self, c):
        if c<0:
            return(-1e16)
        else:
            if self.sigma > 0:
                return ((c**(1-self.sigma)-1)/(1-self.sigma)) # CRRA
            else:
                return (-0.5*(c-self.c_bar)**2) # Quadratic

    def der_u(self, c):
        if c<0:
            return(-1e16)
        else:   
            if self.sigma > 0:
                return (c**(-self.sigma)) # CRRA
            else:
                return (-(c-self.c_bar)) # Quadratic
       
    def Consumption(self, assets, y):
        Consumption = []
        for i in range(len(assets)-1):
            Consumption.append(y+(1+self.r)*assets[i]-assets[i+1])
        return(Consumption)


    def m(self, a0,a1, y):
        if a1<self.A_min:
            return -1e16
        else:
            c = y-a1+a0*(1+self.r)
            if c >= 0:
                return self.u(c)
            else:
                return (-1e16)
            
    def mat_M(self):
        M = [0]*(self.p*len(self.Y))
        for i in range(self.p*len(self.Y)):
            M[i] = [0] * self.p
        for z in range(len(self.Y)):
            for i in range(self.p):
                for j in range(self.p):
                    index = i+self.p*z
                    M[index][j] = self.m(self.A[i], self.A[j], self.Y[z])
        return(M)            

    def mat_W(self, V):
        W = [0]*len(self.Y)
        for i in range(len(self.Y)):
            W[i] = [0] * self.p
        for j in range(self.p):
            for i in range(len(self.Y)):
                W[i][j] = 0
                for x in range(len(self.Y)):
                    W[i][j] = W[i][j] + self.pi[i][x]*V[-1][j]
        return(W)
    
    def XI(self, M, W):
        Xi = [0]*(self.p*len(self.Y))  
        for i in range(self.p*len(self.Y)):
            Xi[i] = [0] * self.p
                
        for j in range(self.p):
            for i in range(self.p*len(self.Y)):
                Xi[i][j] = 0
                for x in range(len(self.Y)):
                    Xi[i][j] = M[i][j] + self.beta*W[x][j]   
        return(Xi)    
    

    def new_V(self, V, Xi, a_dec):
        V_new = [0]*(self.p*len(self.Y))
        a_new = [0]*(self.p*len(self.Y))
        for i in range(self.p*len(self.Y)):
            V_new[i]=max(Xi[i])
            a_new[i] = np.argmax(Xi[i])
        V.append(V_new)
        a_dec.append(a_new)
        return(V, a_dec)
        
    def VFI(self):
        V = []
        a_dec = []
        V.append(list(np.zeros(self.p*len(self.Y))))
        V.append(list(np.zeros(self.p*len(self.Y))))
        a_dec.append(list(np.zeros(self.p*len(self.Y))))
        a_dec.append(list(np.zeros(self.p*len(self.Y))))
        M = self.mat_M()
        W = self.mat_W(V) 
        Xi = self.XI(M, W) 
        V, a_dec = self.new_V(V, Xi, a_dec) 
        while (np.linalg.norm(np.asarray(V[-1]) - np.asarray(V[-2])) > 1):
            W = self.mat_W(V) 
            Xi = self.XI(M, W) 
            V, a_dec = self.new_V(V, Xi, a_dec) 
        return(V, a_dec)
    
    def new_C(self):
        V, a_dec = self.VFI()
        c_dec = [0]*len(self.Y)
        a_dec_t = [0]*len(self.Y)
        for i in range(len(self.Y)):
            c_dec[i]=[0]*len(self.A)
            a_dec_t[i]=[0]*len(self.A)
        for y in range(len(self.Y)):
            for a in range(len(self.A)):
                idx=int(a+y*len(self.A))
                c_dec[y][a] = self.Y[y] + (1+self.r)*self.A[a] - self.A[a_dec[-1][idx]]
                a_dec_t[y][a] = self.A[a_dec[-1][idx]]
        return(c_dec, a_dec_t)

#------------------------------------------------------------------------------------------------------------------------------
# II.4 Partial Equilibrium

#   II.4.1 With certainty

CRRA_cer = HW4(0.06, 0.04, 0.1, 100, 0, 0, 1)
CRRA_cer2 = HW4(0.06, 0.04, 0.2, 100, 0, 0, 1)
V_T, a_dec_T = CRRA_cer.VFI()
c_cer, a_cer = CRRA_cer.new_C()
c_cer2, a_cer2 = CRRA_cer2.new_C()

plt.plot(CRRA_cer.A, c_cer[1], label="CRRA, sigma=0.1")
plt.plot(CRRA_cer.A, c_cer2[1], label="CRRA, sigma=0.2")
plt.legend()
plt.show()

quad_cer = HW4(0.06, 0.04, 0, 100, 0, 0, 1)
quad_cer2 = HW4(0.06, 0.04, 0, 200, 0, 0, 1)
V_TEST_quad, a_dec_TEST_quad = quad_cer.VFI()
c_cer_quad, a_cer_quad = quad_cer.new_C()
c_cer_quad2, a_cer_quad2 = quad_cer2.new_C()

plt.plot(quad_cer.A, c_cer_quad[1], label="Quadratic, c_bar=100")
plt.plot(quad_cer.A, c_cer_quad2[1], label="Quadratic, c_bar=200")
plt.legend()
plt.show()

#   II.4.2 With uncertainty
quad_uncer = HW4(0.06, 0.04, 0, 100, 0, 0.1, 1)
V_quad_2, a_uncer_quad = quad_uncer.VFI()
c_cer_quad2, a_cer_quad2 = quad_uncer.new_C()

c_cer_quad2_n = [0]*2
c_cer_quad2_n[0] = np.concatenate(([0]*(len(quad_cer.A)-len(quad_uncer.A)), c_cer_quad2[0]), axis=0)
c_cer_quad2_n[1] = np.concatenate(([0]*(len(quad_cer.A)-len(quad_uncer.A)), c_cer_quad2[1]), axis=0)

plt.plot(quad_uncer.A, c_cer_quad2_n[0][0:(len(quad_uncer.A))], label="Quadratic, sigma_y=0.1, (-) shock")
plt.plot(quad_uncer.A, c_cer_quad2_n[1][0:(len(quad_uncer.A))], label="Quadratic, sigma_y=0.1, (+) shock")
plt.plot(quad_uncer.A, c_cer_quad[1][0:(len(quad_uncer.A))], label="Quadratic, sigma_y=0, certainty")
plt.legend()
plt.show()

a_cer_quad2_n = [0]*2
a_cer_quad2_n[0] = np.concatenate(([0]*(len(quad_cer.A)-len(quad_uncer.A)), a_cer_quad2[0]), axis=0)
a_cer_quad2_n[1] = np.concatenate(([0]*(len(quad_cer.A)-len(quad_uncer.A)), a_cer_quad2[1]), axis=0)

plt.plot(quad_uncer.A, a_cer_quad2[0][0:(len(quad_uncer.A))], label="(-) shock")
plt.plot(quad_uncer.A, a_cer_quad2[1][0:(len(quad_uncer.A))], label="(+) shock")
plt.plot(quad_uncer.A, a_cer_quad[0][0:(len(quad_uncer.A))], label="certainty")
plt.legend()
plt.show()
