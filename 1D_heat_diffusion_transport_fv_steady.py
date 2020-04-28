import numpy as np
import matplotlib.pylab as plt

# finite volume method

# pro-processing 

# define domain
L = 2.0 # unit:m

# discretization or "meshing"
N_x = 31 # cell number
x_P_i = np.linspace(0,L,N_x)
delt_x = L/(N_x-1)

phi_i = np.zeros(N_x)

phi_i[0]=373.0 # unit: K
phi_i[-1]=273.0 # unit: K

a_E_i = np.zeros(N_x)
a_W_i = np.zeros(N_x)
a_P_i = np.zeros(N_x)

k_e_i = np.zeros(N_x)
k_w_i = np.zeros(N_x)

# diffusivity
k_P_i = np.ones(N_x)*0.05 # unit:m^2/s

# source S=SC+SP*T
SC = 2.0
SP = 0.01

for i in range(1,len(k_P_i)-1):    
    k_e_i[i] = 2.0*k_P_i[i]*k_P_i[i+1]/(k_P_i[i]+k_P_i[i+1])
    k_w_i[i] = 2.0*k_P_i[i]*k_P_i[i-1]/(k_P_i[i]+k_P_i[i-1])
    a_E_i[i] = k_e_i[i]/delt_x
    a_W_i[i] = k_w_i[i]/delt_x
    a_P_i[i] = a_W_i[i]+a_E_i[i]-SP*delt_x

a_ij = np.zeros([N_x,N_x])

b = np.ones(N_x)*SC*delt_x

b[1] = a_W_i[1]*phi_i[0]
b[-2] = a_E_i[-2]*phi_i[-1]

for i in range(1,N_x-1):
    a_ij[i,i] = a_P_i[i]
    if i!=1:
        a_ij[i,i-1] = -a_W_i[i-1]
    if i!=N_x-2:
        a_ij[i,i+1] = -a_E_i[i+1]

print a_ij

# note: a[-1] means the last element, 
# but a[:-1] means the elements from the first to the last SECOND, NOT the last 
linalg_x = np.linalg.solve(a_ij[1:-1,1:-1], b[1:-1])  
    
phi_i[1:-1] = linalg_x

print phi_i
print a_P_i
    
# Post-processing or visualize the results        
fig,ax = plt.subplots(1,figsize=[6,4])
ax.scatter(x_P_i,phi_i)
ax.set_xlim(0,L)
ax.set_xlabel('Distance[m]')
ax.set_ylabel('Mass fraction')

plt.show()
