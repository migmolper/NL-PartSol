#!/opt/anaconda3/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FixedLocator, FormatStrFormatter

def MPM(L, Xp, Xi):
    """"""
    if (np.abs(Xp-Xi) >= L):
        Sip = 0
    elif (-L <= Xp-Xi <= 0):
        Sip = 1-np.abs(Xp-Xi)/L
    elif (L >= Xp-Xi >=0):
        Sip = 1-np.abs(Xp-Xi)/L
    return Sip

def d_MPM(L, Xp, Xi):
    """"""
    if (np.abs(Xp-Xi) >= L):
        Sip = 0
    elif (-L <= Xp-Xi <= 0):
        Sip = 1/L
    elif (L >= Xp-Xi >=0):
        Sip = -1/L
    return Sip

def print_dx_MPM():
    N = 50
    xp = yp = np.linspace(-2, 2, N)
    Xp, Yp = np.meshgrid(xp, yp, indexing='ij', sparse=False)
    Sip_dxdy = np.zeros([N,N])
    xi = yi = 0
    L = 1.
    format_fig = 'png'
    Show_grafs = False

    for i in range(0,N):
        for j in range(0,N):
            Sip_dxdy[i][j] = d_MPM(L,xp[i],xi) * MPM(L,yp[j],yi)
           
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
    
    surf = ax.plot_surface(Xp, Yp, Sip_dxdy, rstride=1,
                           cstride=1, cmap=cm.jet,
                           linewidth=0, antialiased=False)
    ax.set_zlim3d(-1.00, 1.00)
    plt.tight_layout()
    plt.savefig('MPM_Shape_Fun_dx.%s'%(format_fig))
    if Show_grafs:
        plt.show()
    fig.clear()

def print_dy_MPM():
    N = 50
    xp = yp = np.linspace(-2, 2, N)
    Xp, Yp = np.meshgrid(xp, yp, indexing='ij', sparse=False)
    Sip_dxdy = np.zeros([N,N])
    xi = yi = 0
    L = 1.
    lp = 1/4.
    format_fig = 'png'
    Show_grafs = False

    for i in range(0,N):
        for j in range(0,N):
            Sip_dxdy[i][j] = MPM(L,xp[i],xi) * d_MPM(L,yp[j],yi) 
            
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
        
    surf = ax.plot_surface(Xp, Yp, Sip_dxdy, rstride=1,
                           cstride=1, cmap=cm.jet,
                           linewidth=0, antialiased=False)
    ax.set_zlim3d(-1.00, 1.00)
    plt.tight_layout()
    plt.savefig('MPM_Shape_Fun_dy.%s'%(format_fig))
    if Show_grafs:
        plt.show()
    fig.clear()
    
def print_MPM():
    N = 50
    xp = yp = np.linspace(-2, 2, N)
    Xp, Yp = np.meshgrid(xp, yp, indexing='ij', sparse=False)
    Sip_xy = np.zeros([N,N])
    xi = yi = 0
    L = 1.
    lp = 1/4.
    format_fig = 'png'
    Show_grafs = False

    for i in range(0,N):
        for j in range(0,N):
            Sip_xy[i][j] = MPM(L,xp[i],xi) * MPM(L,yp[j],yi)
            

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

    surf = ax.plot_surface(Xp, Yp, Sip_xy, rstride=1,
                           cstride=1, cmap=cm.jet,
                           linewidth=0, antialiased=False)
    
    ax.set_zlim3d(0.00, 1.00)
    ax.set_zticks([0, 0.5, 1])
    plt.tight_layout()
    plt.savefig('MPM_Shape_Fun.%s'%(format_fig))
    if Show_grafs:
        plt.show()
    fig.clear()

print_MPM()
print_dx_MPM()
print_dy_MPM()
