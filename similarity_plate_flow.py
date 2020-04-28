import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from pycse import bvp

def main():
    # The equation is going to be solved:
    # [@$\frac{1}{2} f f^{''} + f^{'''} = 0$ @]
    # Let 
    # [@$ f_1 = f $ @]
    # [@$f_2 = f_1^{'} = f^{'} $ @]
    # [@$f_3 = f_2^{'} = f^{''} $ @]
    # Then we have
    # [@$f_3^{'} = -\frac{1}{2} f_1 f_3$@]
    # subject to the boundary conditions 
    # [@$f_1(0) = 0$ @]
    # [@$f_2(0) = 0$ @]
    # [@$f_2(\infty) = 1$ @]
    
    ### Solve the Fluid Flow Equation

    eta_range = 6.5
    sol = plate_solver()
    eta = np.linspace(0.0,eta_range,100)
    sol.solve_f(eta)

    fig1,ax1 = plt.subplots(1,figsize=[6,4])
    ax1.plot(sol.eta,sol.f)
    ax1.set_xlim([0,eta_range])
    ax1.set_ylim([0,5])
    ax1.set_xlabel(r'$\eta$')
    ax1.set_ylabel(r'$f$')
    fig1.savefig('../figs/f_eta.pdf')
    #plt.show(fig1)

    fig2,ax2 = plt.subplots(1,figsize=[6,4])
    ax2.plot(sol.eta,sol.dfdeta)
    ax2.set_xlim([0,eta_range])
    ax2.set_ylim([0,1])
    ax2.set_xlabel(r'$\eta$')
    ax2.set_ylabel(r'$f\,^\prime$')
    fig2.savefig('../figs/f_prime_eta.pdf')
    #plt.show(fig2)

    fig3,ax3 = plt.subplots(1,figsize=[6,4])
    ax3.plot(sol.eta,sol.d2fdeta2)
    ax3.set_xlim([0,eta_range])
    ax3.set_ylim([0,.4])
    ax3.set_xlabel(r'$\eta$')
    ax3.set_ylabel(r'$f\,^{\prime \prime}$')
    fig3.savefig('../figs/f_double_prime_eta.pdf')
    #plt.show(fig3)

    ### Solve the Thermal Equation

    # [@$\frac{d^2\theta}{d\eta^2}+\frac{Pr}{2}\frac{d\theta}{d\eta}=0$@]
    # 
    # subject to the boundary conditions:
    # 
    # [@ $\theta(0)=1$ @]
    # 
    # [@$\theta(\infty)=0$@]
    
    eta_T_range=15.0
    eta_T = np.linspace(0.0,eta_T_range,100)
    Pr = 0.093
    sol.solve_theta(eta_T,Pr)

    fig1,ax1 = plt.subplots(1,figsize=[6,4])
    ax1.plot(sol.eta_T,sol.theta)
    ax1.set_xlabel(r'$\eta$')
    ax1.set_ylabel(r'$\theta$')
    ax1.set_xlim([0,eta_T_range])
    ax1.set_ylim([0,1])
    fig1.savefig('../figs/theta_eta.pdf')
    #plt.show(fig1)

    fig2,ax2 = plt.subplots(1,figsize=[6,4])
    ax2.plot(sol.eta_T,sol.dthetadeta)
    ax2.set_xlabel(r'$\eta$')
    ax2.set_ylabel(r'$\theta \, ^\prime$')
    ax2.set_xlim([0,eta_T_range])
    ax2.set_ylim([-0.14,0])
    fig2.savefig('../figs/theta_prime_eta.pdf')
    #plt.show(fig2)

    ### Visual Velocity Profile (set $U_{\infty}=0.06m/s$)

    x_grid = np.linspace(0.0001,0.005,50)
    y_grid = np.linspace(0,0.005,50)

    U_infty = 0.06
    rho = 4120.0
    vis = 0.003
    nu = vis/rho

    x_grid,y_grid = np.meshgrid(x_grid,y_grid)
    row,col = len(x_grid[:,0]),len(x_grid[0,:])
    U_x_grid = np.empty([row,col])

    fintp = interp1d(sol.eta,sol.dfdeta,kind='linear')

    for i in range(row):
        for j in range(col):
            eta_tmp = np.sqrt(U_infty/nu/x_grid[i,j])*y_grid[i,j]
            #print eta_tmp
            U_x_grid[i,j] = U_infty*fintp(eta_tmp)
    #print eta_tmp
    #print fintp(eta_tmp)
    fig,ax = plt.subplots(1)
    cf = ax.contourf(x_grid,y_grid,U_x_grid)
    plt.title('velocity')
    plt.colorbar(cf)

    ### Visual Temperature Profile

    x_grid = np.linspace(0.00001,0.005,500)
    y_grid = np.linspace(0,0.005,50)

    T_infty = 1770.0 # in C
    T_w = 1667.0
    k = 31.4
    U_infty = 0.06
    rho = 4120.0
    vis = 0.003
    nu = vis/rho

    x_grid,y_grid = np.meshgrid(x_grid,y_grid)
    row,col = len(x_grid[:,0]),len(x_grid[0,:])
    T_grid = np.empty([row,col])

    thetaintp = interp1d(sol.eta_T,sol.theta,kind='linear')

    for i in range(row):
        for j in range(col):
            eta_tmp = np.sqrt(U_infty/nu/x_grid[i,j])*y_grid[i,j]
            #print eta_tmp
            T_grid[i,j] = thetaintp(eta_tmp)*(T_w-T_infty)+T_infty

    fig,ax = plt.subplots(1)
    cf = ax.contourf(x_grid,y_grid,T_grid)
    plt.title('Temperature')
    plt.colorbar(cf)


    ### Calculate the Heat Transfer Coefficient

    h = k*(T_grid[0,:] - T_grid[1,:])/(y_grid[1,0]-y_grid[0,0])/(T_w-T_infty)

    fig,ax = plt.subplots(1)
    ax.plot(x_grid[0,:],h)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$h$')

    print "the average h is: {0}".format(sum(h)/len(h))


class plate_solver():
    def __init(self):
        pass
    
    def solve_f(self,eta,**kwargs):
        self.eta = eta
        
        def odefun(F, x):
            f1, f2, f3 = F
            return [f2,  
                    f3,             
                    -0.5 * f1 * f3] 
        
        def bcfun(Fa, Fb):
            return [Fa[0],        # f1(0) =  0
                    Fa[1],        # f2(0) = 0
                    1.0 - Fb[1]]  # f2(inf) = 1
        
        #initial guess (need to be carefully guessed)
        f1init = self.eta
        f2init = np.exp(self.eta)
        f3init = np.exp(-self.eta)
        
        #solve equation
        Finit = np.vstack([f1init, f2init, f3init])
        
        sol = bvp(odefun, bcfun, self.eta, Finit)
        
        self.f = sol[0]
        self.dfdeta = sol[1]
        self.d2fdeta2 = sol[2]
        
        #extrapolation
        self.eta[-1] = 200.0
        if kwargs.get('eta_extra'):
            self.eta[-1] = kwargs['eta_extra']
        self.dfdeta[-1] = 1.0
        self.f[-1] = (self.eta[-1]-self.eta[-2])*self.dfdeta[-1]
        
    def solve_theta(self,eta_T,Pr,**kwargs):
        self.Pr = Pr
        self.eta_T = eta_T
        
        f_intp = interp1d(self.eta,self.f,kind='linear')
        
        def odefun(F, x):
            f1, f2 = F
            df1dx = f2
            df2dx = -0.5*Pr*f_intp(x)*f2
            return [df1dx, df2dx]
        
        def bcfun(Fa, Fb):
            f1a, f2a = Fa
            f1b, f2b = Fb
            return [f1a-1.0,        # f1_a =  1
                    f1b-0.0]        # f1_b =  0
        
        f1init = np.exp(-self.eta_T)
        f2init = -np.exp(-self.eta_T)
        Finit = np.vstack([f1init, f2init])
        
        sol = bvp(odefun, bcfun, self.eta_T, Finit)
        self.theta = sol[0]
        self.dthetadeta = sol[1]
        #extrapolation    
        self.eta_T[-1] = 10000.0
        if kwargs.get('eta_extra'):
            self.eta_T[-1] = kwargs['eta_extra']
        self.dthetadeta[-1] = 0.0
        self.theta[-1] = (self.eta_T[-1]-self.eta_T[-2])*self.dthetadeta[-1]

main()
