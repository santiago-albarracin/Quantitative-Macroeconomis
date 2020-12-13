#   Quantitative Macroeconomics
#   HW 4

# Import packages
import numpy as np
from numba import jit
import matplotlib.pyplot as plt
from quantecon.markov import DiscreteDP

#------------------------------------------------------------------------------------------------------------------------------
# II.5 Aiyagari model

"""
All the used code is based on:
Quantitative Economics with Python. Sargent and Stachurski. 2020
"""


class Household:
    def __init__(self,
                 r=0.04,
                 w=1.0,                 
                 β=1/(1+0.04),               
                 a_min=0,
                 Π=[[0.05, 0.15, 0.2, 0.2, 0.2, 0.15, 0.05], [0.05, 0.15, 0.2, 0.2, 0.2, 0.15, 0.05], [0.05, 0.15, 0.2, 0.2, 0.2, 0.15, 0.05],[0.05, 0.15, 0.2, 0.2, 0.2, 0.15, 0.05],[0.05, 0.15, 0.2, 0.2, 0.2, 0.15, 0.05],[0.05, 0.15, 0.2, 0.2, 0.2, 0.15, 0.05],[0.05, 0.15, 0.2, 0.2, 0.2, 0.15, 0.05]],
                 z_vals=[0.01, 0.2, 0.4, 0.6, 0.8, 1, 1.2],
                 a_max=18,
                 a_size=100):  
        
        self.r, self.w, self.β = r, w, β
        self.a_min, self.a_max, self.a_size = a_min, a_max, a_size

        self.Π = np.asarray(Π)
        self.z_vals = np.asarray(z_vals)
        self.z_size = len(z_vals)

        self.a_vals = np.linspace(a_min, a_max, a_size)
        self.n = a_size * self.z_size

        self.Q = np.zeros((self.n, a_size, self.n))
        self.build_Q()

        self.R = np.empty((self.n, a_size))
        self.build_R()
        
    def set_prices(self, r, w):
        """
        Use this method to reset prices.  Calling the method will trigger a
        re-build of R.
        """
        self.r, self.w = r, w
        self.build_R()

    def build_Q(self):
        populate_Q(self.Q, self.a_size, self.z_size, self.Π)

    def build_R(self):
        self.R.fill(-np.inf)
        populate_R(self.R, self.a_size, self.z_size, self.a_vals, self.z_vals, self.r, self.w)


@jit(nopython=True)
def populate_R(R, a_size, z_size, a_vals, z_vals, r, w):
    n = a_size * z_size
    for s_i in range(n):
        a_i = s_i // z_size
        z_i = s_i % z_size
        a = a_vals[a_i]
        z = z_vals[z_i]
        for new_a_i in range(a_size):
            a_new = a_vals[new_a_i]
            c = w * z + (1 + r) * a - a_new
            if c > 0:
                R[s_i, new_a_i] = np.log(c)

@jit(nopython=True)
def populate_Q(Q, a_size, z_size, Π):
    n = a_size * z_size
    for s_i in range(n):
        z_i = s_i % z_size
        for a_i in range(a_size):
            for next_z_i in range(z_size):
                Q[s_i, a_i, a_i * z_size + next_z_i] = Π[z_i, next_z_i]


@jit(nopython=True)
def asset_marginal(s_probs, a_size, z_size):
    a_probs = np.zeros(a_size)
    for a_i in range(a_size):
        for z_i in range(z_size):
            a_probs[a_i] += s_probs[a_i * z_size + z_i]
    return a_probs

@jit(nopython=True)
def consumption_marginal(s_probs, a_size, z_size):
    c_probs = np.zeros(a_size)
    for a_i in range(a_size):
        for z_i in range(z_size):
            c_probs[a_i] += s_probs[a_i * z_size + z_i]
    return c_probs

#Parameters
A = 1.0 
N = 1.0
α = 0.33 
β = 1/(1+0.04)
δ = 0.05  

def r_to_w(r):
    return A * (1 - α) * (A * α / (r + δ))**(α / (1 - α))

def rd(K):
    return A * α * (N / K)**(1 - α) - δ

def prices_to_capital_stock(am, r):
    w = r_to_w(r)
    am.set_prices(r, w)
    aiyagari_ddp = DiscreteDP(am.R, am.Q, β)
    # Compute the optimal policy
    results = aiyagari_ddp.solve(method='policy_iteration')
    # Compute the stationary distribution
    stationary_probs = results.mc.stationary_distributions[0]
    # Extract the marginal distribution for assets
    asset_probs = asset_marginal(stationary_probs, am.a_size, am.z_size)
    # Return K
    return np.sum(asset_probs * am.a_vals)

am = Household(a_max=30)

am_ddp = DiscreteDP(am.R, am.Q, am.β)

num_points = 20
r_vals = np.linspace(0.005, 0.04, num_points)

k_vals = np.empty(num_points)
for i, r in enumerate(r_vals):
    k_vals[i] = prices_to_capital_stock(am, r)

fig, ax = plt.subplots()
ax.plot(k_vals, r_vals, lw=2, alpha=0.6, label='Supply')
ax.plot(k_vals, rd(k_vals), lw=2, alpha=0.6, label='Demand')
ax.set_xlabel('K')
ax.set_ylabel('r')
ax.legend()
