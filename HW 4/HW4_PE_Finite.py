#   Quantitative Macroeconomics
#   HW 4

# Import packages
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#------------------------------------------------------------------------------------------------------------------------------
# I) A simple wealth model

class HW4:
    def __init__(self, rho, r, sigma, c_bar, gamma, sigma_y, T):
        self.rho = rho
        self.r = r
        self.sigma = sigma
        self.c_bar = c_bar
        self.gamma = gamma
        self.sigma_y = sigma_y
        self.beta = 1/(1+rho)
        self.Y = [1-sigma_y, 1+sigma_y]	#Change to 1
        self.pi = [[(1+gamma)/2, (1-gamma)/2], [(1-gamma)/2, (1+gamma)/2]]
        if T==0:
            self.A_min = -min(self.Y)*((1+r)/r)
        else:
            self.A_min = -min(self.Y)*((1+r-(r+1)**(-T))/r)
        self.A = np.linspace(0.1, 20, 100)
        self.T = T
            
    def u(self, c):
        if c<0:
            return(-1e16)
        else:
            if self.sigma != 0:
                return ((c**(1-self.sigma)-1)/(1-self.sigma)) # CRRA
            else:
                return (-0.5*(c-self.c_bar)**2) # Quadratic
        
    def der_u(self, c):
        if c<0:
            return(-1e16)
        else:   
            if self.sigma != 0:
                return (c**(-self.sigma)) # CRRA
            else:
                return (-c) # Quadratic
    
    def Eu(self, a0, a1, a_past, y):
        for i in range(len(self.pi)):
            if y == self.Y[i]:
                row = i
        sum = 0
        for i in range(len(self.pi)):
            sum = self.pi[row][i]*self.der_u(self.Y[i]+(1+self.r)*a1-a_past) 
        x = self.der_u(y+(1+self.r)*a0-a1) - self.beta * (1 + self.r) * sum
        return (x)
    
    def solution_Eu(self, a0, a1, y):
        a1_opt = list(fsolve(lambda a_past: self.Eu(a0, a1, a_past, y), a1))
        return(a1_opt[0])
    
    def trans(self, a1, a0, y, T):
        aT = 0 
        A = []
        for i in range(T+1):
            if i == 0:
                a = self.der_u(y + (1+self.r)*a0-a1[i]) - self.beta*(1+self.r)*y*self.der_u(y+(1+self.r)*a1[i]-a1[i+1])        
            elif i == T:
                a = self.der_u(y + (1+self.r)*a1[i-1]-a1[i]) - self.beta*(1+self.r)*y*self.der_u(y+(1+self.r)*a1[i]-aT)           
            else:
                a = self.der_u(y + (1+self.r)*a1[i-1]-a1[i]) - self.beta*(1+self.r)*y*self.der_u(y+(1+self.r)*a1[i]-a1[i+1])        
            A.append(a)
        return(A)
              
    def assets(self, a0, y, aT):
        T = self.T
        if self.T > 0:       
            a10 = np.linspace(a0, aT, T+1) #first guess - linear
            assets = fsolve(lambda a1: self.trans(a1, a0, y, T), a10)
            return(assets)    
        else:
            a1 = []
            a1.append(a0) # init guesses of a"0
            a1.append(self.solution_Eu(a0, a1[-1], y))            
            while abs(a1[-1]-a1[-2])>0.05:
                a1.append(self.solution_Eu(a1[-2], a1[-1], y))
            aT_SS = a1[-1]
            a10 = np.linspace(a0, aT_SS, 30+1) #first guess - linear
            assets = fsolve(lambda a1: self.trans(a1, a0, y, 30), a10)
            return(assets)              
        
    def Consumption(self, assets, y):
        Consumption = []
        for i in range(len(assets)-1):
            Consumption.append(y+(1+self.r)*assets[i]-assets[i+1])
        return(Consumption)


#------------------------------------------------------------------------------------------------------------------------------
# II.4 Partial Equilibrium

#   II.4.1 With certainty

CRRA_cer = HW4(0.06, 0.04, 0.5, 100, 0, 0, 45)
assets_CRRA = CRRA_cer.assets(2, 1, 0)
Consumption_CRRA_cer = CRRA_cer.Consumption(assets_CRRA,1)

plt.plot(range(len(assets_CRRA)), assets_CRRA, label="Assets")
plt.plot(range(len(Consumption_CRRA_cer)),Consumption_CRRA_cer, label="Consumption")
plt.legend()
plt.show()

quad_cer = HW4(0.06, 0.04, 0, 100, 0, 0, 45)
assets_quad = quad_cer.assets(2, 1, 0)
Consumption_quad_cer = CRRA_cer.Consumption(assets_quad,1)

plt.plot(range(len(assets_quad)), assets_quad, label="Assets")
plt.plot(range(len(Consumption_quad_cer)),Consumption_quad_cer, label="Consumption")
plt.legend()
plt.show()

c_5 = [0]*len(CRRA_cer.A)
c_40 = [0]*len(CRRA_cer.A)
for a in range(len(CRRA_cer.A)):
    assets_path = CRRA_cer.assets(CRRA_cer.A[a], 1, 0)
    Consumption_path = CRRA_cer.Consumption(assets_path,1)
    plt.plot(range(45), Consumption_path)
    c_5[a] = Consumption_path[6]
    c_40[a] = Consumption_path[41]
plt.show()
    
plt.plot(CRRA_cer.A, c_5, label="age 5")
plt.plot(CRRA_cer.A, c_40, label="age 40")
plt.legend()
plt.show()

c_5_quad = [0]*len(quad_cer.A)
c_40_quad = [0]*len(quad_cer.A)
for a in range(len(quad_cer.A)):
    assets_path = quad_cer.assets(quad_cer.A[a], 1, 0)   
    Consumption_path = quad_cer.Consumption(assets_path,1)
    plt.plot(range(45), Consumption_path)
    c_5_quad[a] = Consumption_path[6]
    c_40_quad[a] = Consumption_path[41]
plt.show()

plt.plot(quad_cer.A, c_5_quad, label="age 5")
plt.plot(quad_cer.A, c_40_quad, label="age 40")
plt.legend()
plt.show()