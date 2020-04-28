import numpy as np
import matplotlib.pylab as plt

# pro-processing 

# define domain
L = 2.0 # unit:m

# material properties
D = 0.05 # unit:m^2/s
v_x = 0.5 # velocity

# discretization or "meshing"
N_x = 101 # node number
x_i = np.linspace(0,L,N_x)
delt_x = L*1.0/(N_x-1)

# time step and total computational time
delt_t = 0.002 # unit:s
N_t = 2001
t_total=delt_t*N_t
t_i = np.linspace(0,t_total,N_t)

#print 'x_i={0}'.format(x_i)

# field variable mass fraction 
mf_ij = np.zeros([N_x,N_t]) # i for x and j for t
# create a temporary mf_ij for applying the periodic boundary condition
mf_tmp_ij = np.zeros([N_x,N_t]) # i for x and j for t

# set initial conditions
for i in range(N_x):
    mf_ij[i,0]=0.5*np.sin(2*np.pi*delt_x*(i-1))

#print 'mf_i={0}'.format(mf_i)

# solver
for j in range(N_t-1):  
    mf_tmp_ij=mf_ij # store the solution
    
    for i in range(1,N_x-1):
        mf_ij[i,j+1] = D*1.0*(mf_tmp_ij[i+1,j]-2.0*mf_tmp_ij[i,j]+\
            mf_tmp_ij[i-1,j])*delt_t/delt_x**2.0+mf_tmp_ij[i,j]-\
            v_x*1.0*(mf_tmp_ij[i+1,j]-mf_tmp_ij[i-1,j])*delt_t/(2.0*delt_x)
    
    # periodic boundary condition
    mf_ij[-1,j+1] = D*1.0*(mf_tmp_ij[1,j]-2.0*mf_tmp_ij[-1,j]+\
        mf_tmp_ij[-2,j])*delt_t/delt_x**2.0+mf_tmp_ij[i,j]-\
        v_x*1.0*(mf_tmp_ij[1,j]-mf_tmp_ij[-1,j])*delt_t/(2.0*delt_x)
    
    mf_ij[0,j+1] = mf_ij[-1,j+1] 
    
# Post-processing or visualize the results        
fig,ax = plt.subplots(1,figsize=[6,4])
markersize=40
ax.scatter(x_i,mf_ij[:,0],marker='o',color='k',facecolor='none',s=markersize,label='{0:1.2f}s'.format(0))
ax.scatter(x_i,mf_ij[:,(N_t-1)/40],marker='s',color='r',facecolor='none',s=markersize,label='{0:1.2f}s'.format(delt_t*(N_t-1)/40))
ax.scatter(x_i,mf_ij[:,(N_t-1)/20],marker='d',color='b',facecolor='none',s=markersize,label='{0:1.2f}s'.format(delt_t*(N_t-1)/20))
ax.scatter(x_i,mf_ij[:,(N_t-1)/10],marker='^',color='c',facecolor='none',s=markersize,label='{0:1.2f}s'.format(delt_t*(N_t-1)/10))
ax.scatter(x_i,mf_ij[:,(N_t-1)/6.7],marker='>',color='g',facecolor='none',s=markersize,label='{0:1.2f}s'.format(delt_t*(N_t-1)/6.7))

ax.legend(loc='upper right')
ax.set_xlim(0,L)
ax.set_xlabel('Distance[m]')
ax.set_ylabel('Mass fraction')

fig.savefig('../figs/1d_mass_species_demo.pdf',bbox_inches='tight')

plt.show()
