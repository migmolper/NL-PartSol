#!/home/migmolper2/anaconda3/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FixedLocator, FormatStrFormatter

"""
  Shape functions based in :
  " Local maximum-entropy approximation schemes : a seamless 
  bridge between finite elements and meshfree methods "
  by M.Arroyo and M.Ortiz, 2006.

  Here we employ the same nomenclature as in the paper. With the single
  different of the "l" variable wich represents the distances between the
  evaluation point and the neighborhood nodes.

  List of functions :
  - LME_lambda_NR
  - LME_fa
  - LME_p
  - LME_r
  - LME_J
  - LME_dp
"""

def LME_fa(l_a, lambda_a, Beta_a):
    """
    """
    return - Beta_a * np.dot(l_a, l_a) + np.dot(l_a,lambda_a)

def LME_p(l, lambda_a, Beta_a):
    """
    Get the value of the shape function "pa" (1 x neighborhood) in the
    neighborhood nodes.
    
    Input parameters :
    -> l : Matrix with the distances to the
    neighborhood nodes (neiborghood x dim).
    -> lambda : Initial value of the lagrange multipliers (1 x dim).
    -> Beta : Tunning parameter (scalar).
    """
    
    # Definition of some parameters
    N_a = np.shape(l)[0]
    N_dim = 2
    Z = 0.0
    Z_m1 = 0.0

    # Vector with the values of the shape-function in the nodes
    p = np.zeros([N_a])

    # Get Z and the numerator
    for a in range(0,N_a):
        la = l[a,:]
        p[a] = np.exp(LME_fa(la,lambda_a,Beta_a))
        Z += p[a]

    # Get the inverse of Z
    Z_m1 = 1/Z

    # Divide by Z and get the final value
    for a in range(0,N_a):
        p[a] *= Z_m1
      
    # Return the value of the shape function
    return p


def LME_r(l,p):
    """
    Get the gradient "r" (dim x 1) of the function log(Z) = 0.
    Input parameters :
    -> l : Matrix with the distances to the 
    neighborhood nodes (neighborhood x dim).
    -> p : Shape function value in the
    neighborhood nodes (1 x neighborhood).
    """
    
    # Definition of some parameters
    N_a = np.shape(l)[0]
    N_dim = 2

    # Value of the gradient
    r = np.zeros([N_dim])
    
    # Fill ''r''
    for i in range(0,N_dim):
        for a in range(0,N_a):
            r[i] += p[a]*l[a][i]
            
    # Return the value of the gradient
    return r

def LME_J(l, p, r):
    """
    Get the Hessian "J" (dim x dim) of the function log(Z) = 0.
    Input parameters :
    -> l : Matrix with the distances to the
    neighborhood nodes (neighborhood x dim).
    -> p : Shape function value in the
    neighborhood nodes (neighborhood x 1).
    -> r : Gradient of log(Z) (dim x 1).
    """  
    # Definition of some parameters
    N_a = np.shape(l)[0]
    N_dim = 2
  
    # Allocate Hessian */
    J = np.zeros([N_dim,N_dim])

    # Fill the Hessian
    # Get the first component of the Hessian
    for a in range(0,N_a):
        J = J + p[a]*np.tensordot(l[a,:],l[a,:],axes = 0)
    # Get the second value of the Hessian
    J = J - np.tensordot(r,r,axes = 0)
    
    # Return the value of the Hessian
    return J


def LME_dp(l, p):
    """
    Value of the shape function gradient "dp" (dim x neighborhood) in 
    the neighborhood nodes.
    Input parameters :
    -> l : Matrix with the distances to the
    neighborhood nodes (neighborhood x dim).
    -> p : Shape function value in the
    neighborhood nodes (neighborhood x 1).
    """
    
    # Definition of some parameters */
    N_a = np.shape(l)[0]
    N_dim = 2
    dp = np.zeros([N_dim,N_a])
  
    # Get the Gradient and the Hessian of log(Z) 
    r = LME_r(l,p)
    J = LME_J(l,p,r)
  
    # Inverse of the Hessian */
    Jm1 = np.linalg.inv(J)
  
    # Fill the gradient for each node */
    for a in range(0,N_a):
        la = l[a,:] 
        Jm1_la = np.dot(Jm1, la)    
        for i in range(0,N_dim):
            dp[i][a] = - p[a]*Jm1_la[i]
  
    # Return the value of the shape function gradient
    return dp


def LME_lambda_NR(l, lambda_GP, Beta):
    """
    Get the lagrange multipliers "lambda" (1 x dim) for the LME 
    shape function. The numerical method for that is the Newton-Rapson.
    
    Input parameters :
    -> l : Matrix with the distances to the
    neighborhood nodes (neighborhood x dim).
    -> lambda : Initial value of the
    lagrange multipliers (1 x dim).
    -> Beta : Tunning parameter (scalar).
    -> h : Grid spacing (scalar).
    -> TOL_zero : Tolerance for Newton-Rapson.
    """
    
    # Definition of some parameters */
    MaxIter = 100
    Ndim = 2
    NumIter = 0 # Iterator counter
    norm_r = 10.0 # Value of the norm
    TOL_NR = 10e-6 # tolerance
    
    # Start with the Newton-Rapson method
    while norm_r > TOL_NR :
	
        # Get vector with the shape functions evaluated in the nodes
        p = LME_p(l,lambda_GP,Beta)

        # Get the gradient of log(Z) and its norm
        r = LME_r(l,p)
        norm_r = np.linalg.norm(r)

        # Get the Hessian of log(Z) and update it with +||r||*I 
        # according with Dennis M.Kochmann et al. 2019 (CMAME)
        J = LME_J(l,p,r)
        for i in range(0,Ndim):
            J[i][i] += norm_r

        # Check the conditioning number of the Hessian
        if np.linalg.cond(J) > 10 :
            print("Error in LME_lambda_NR. ","Iter :",NumIter,
    	        " The Hessian is near to singular matrix")
            exit()
    
        # Get the increment of lambda
        D_lambda = np.linalg.solve(J,r)

        # Update the value of lambda
        for i  in range(0,Ndim):
            lambda_GP[i] -= D_lambda[i]
    
        # Update the number of iterations
        NumIter += 1
        if NumIter >= MaxIter :
            print("Warning in LME_lambda.","No convergence in",
                  NumIter,MaxIter,"iterations")
            print("Error",norm_r)
            exit()
    
  
    # Once the stopping criteria is reached, 
    # return the lagrange multipliers value
    return lambda_GP

N_gp = 35
N_n = 5
L = 2.
DX = (2*L-0.2)/N_n
Gamma = 17.3 # 17.3
Gamma = ([17.3,10.0,7.0,5.0])
Node = 12 # 12
format_fig='png'
# Show
Show_grafs = False

xp = yp = np.linspace(-L+0.1, L-0.1, N_gp)
zp = np.zeros([1,N_gp])
Xp, Yp = np.meshgrid(xp, yp, indexing='ij', sparse=False)
Zp = np.zeros([N_gp,N_gp])
N_I = np.zeros_like(Zp)
dNdx_I = np.zeros_like(Zp)
dNdy_I = np.zeros_like(Zp)

xI = yI = np.linspace(-L, L, N_n)
zI = np.zeros([1,N_n])
XI, YI = np.meshgrid(xI, yI, indexing='ij', sparse=False)
ZI = np.zeros([N_n,N_n])

l = np.zeros([N_n*N_n,2])

lambda_GP = np.zeros([N_gp*N_gp,2])
Xp_Yp = np.zeros([N_gp*N_gp,2])

for i in range(0,N_gp):
    for j in range(0,N_gp):
        Xp_Yp[i*N_gp + j][0] = Xp[i][j]
        Xp_Yp[i*N_gp + j][1] = Yp[i][j]

for Gamma_i in Gamma:

    # Get the value of gamma
    Beta = Gamma_i/(L*L)
    
    # Distance to the GPa
    for a in range(0,N_gp*N_gp):
        for i in range(0,N_n):
            for j in range(0,N_n):
                l[i*N_n + j][0] = Xp_Yp[a][0] - XI[i][j]
                l[i*N_n + j][1] = Xp_Yp[a][1] - YI[i][j]

        lambda_GP[a,:] = LME_lambda_NR(l, np.zeros([2]), Beta);

        # GPs
        # ax.scatter(XI, YI, ZI, c='r', marker='o')

    for a in range(0,N_gp*N_gp):
        for i in range(0,N_n):
            for j in range(0,N_n):
                l[i*N_n + j][0] = Xp_Yp[a][0] - XI[i][j]
                l[i*N_n + j][1] = Xp_Yp[a][1] - YI[i][j]
            
        N_I_GP = LME_p(l, lambda_GP[a,:], Beta)
        dN_I_GP = LME_dp(l, N_I_GP)
        i = int(a/(N_gp))
        j = a % N_gp
        N_I[i][j] = N_I_GP[Node]
        dNdx_I[i][j] = dN_I_GP[0][Node]
        dNdy_I[i][j] = dN_I_GP[1][Node]


    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.grid(False)
    ax.xaxis.pane.set_edgecolor('black')
    ax.yaxis.pane.set_edgecolor('black')
    ax.zaxis.pane.set_edgecolor('black')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    for line in ax.xaxis.get_ticklines():
        line.set_visible(False)
    for line in ax.yaxis.get_ticklines():
        line.set_visible(False)
    for line in ax.zaxis.get_ticklines():
        line.set_visible(False)
    surf = ax.plot_surface(Xp, Yp, N_I, rstride=1, cstride=1, cmap=cm.jet,
                           linewidth=0, antialiased=False)
    ax.set_zlim3d(0.00, 1.00)
    ax.set_zticks([0, 0.5, 1])
    plt.tight_layout()
    plt.savefig('LME_%1.1f_Shape_Fun.%s'%(Gamma_i,format_fig))  
    if Show_grafs:
        plt.show()
    fig.clear()


    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.grid(False)
    ax.xaxis.pane.set_edgecolor('black')
    ax.yaxis.pane.set_edgecolor('black')
    ax.zaxis.pane.set_edgecolor('black')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    ax.axes.zaxis.set_ticklabels([])
    for line in ax.xaxis.get_ticklines():
        line.set_visible(False)
    for line in ax.yaxis.get_ticklines():
        line.set_visible(False)
    for line in ax.zaxis.get_ticklines():
        line.set_visible(False)    
    surf = ax.plot_surface(Xp, Yp, dNdx_I, rstride=1, cstride=1, cmap=cm.jet,
                           linewidth=0, antialiased=False)
    ax.set_zlim3d(-1.00, 1.00)
    plt.tight_layout()
    plt.savefig('LME_%1.1f_Shape_Fun_dx.%s'%(Gamma_i,format_fig))  
    if Show_grafs:
        plt.show()
    fig.clear()
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.grid(False)
    ax.xaxis.pane.set_edgecolor('black')
    ax.yaxis.pane.set_edgecolor('black')
    ax.zaxis.pane.set_edgecolor('black')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    ax.axes.zaxis.set_ticklabels([])
    for line in ax.xaxis.get_ticklines():
        line.set_visible(False)
    for line in ax.yaxis.get_ticklines():
        line.set_visible(False)
    for line in ax.zaxis.get_ticklines():
        line.set_visible(False)
    surf = ax.plot_surface(Xp, Yp, dNdy_I, rstride=1, cstride=1, cmap=cm.jet,
                        linewidth=0, antialiased=False)
    ax.set_zlim3d(-1.00, 1.00)
    plt.tight_layout()
    plt.savefig('LME_%1.1f_Shape_Fun_dy.%s'%(Gamma_i,format_fig))  
    if Show_grafs:
        plt.show()
    fig.clear()
