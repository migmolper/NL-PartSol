#!/opt/anaconda3/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FixedLocator, FormatStrFormatter

def GIMP(L, lp, Xp, Xi):
    """"""
    if (np.abs(Xp-Xi) >= L+lp):
        Sip = 0
    elif (-L-lp < Xp-Xi <= -L+lp):
        Sip = (1/(4*L*lp))*(L+lp+Xp-Xi)**2
    elif (-L+lp < Xp-Xi <= -lp):
        Sip = 1 + (Xp-Xi)/(L)
    elif (-lp < Xp-Xi <= lp):
        Sip = 1 - ((Xp - Xi)**2 + (lp)**2)/(2*L*lp)
    elif (lp < Xp-Xi <= L-lp):
        Sip = 1 - (Xp-Xi)/(L)
    elif (L-lp < Xp-Xi <= L+lp):
        Sip = (1/(4*L*lp))*(L+lp-Xp+Xi)**2
    return Sip

def d_GIMP(L, lp, Xp, Xi):
    """"""
    if (np.abs(Xp-Xi) >= L+lp):
        d_Sip = 0
    elif (-L-lp < Xp-Xi <= -L+lp):
        d_Sip = (1/(2*L*lp))*(L+lp+Xp-Xi)
    elif (-L+lp < Xp-Xi <= -lp):
        d_Sip = 1/L
    elif (-lp < Xp-Xi <= lp):
        d_Sip = -(Xp - Xi)/(L*lp)
    elif (lp < Xp-Xi <= L-lp):
        d_Sip = -1/L
    elif (L-lp < Xp-Xi <= L+lp):
        d_Sip = -(1/(2*L*lp))*(L+lp-Xp+Xi)
    return d_Sip


def print_dx_GIMP():
    N = 50
    xp = yp = np.linspace(-2, 2, N)
    Xp, Yp = np.meshgrid(xp, yp, indexing='ij', sparse=False)
    Sip_xDy = np.zeros([N,N])
    xi = yi = 0
    L = 1.
    lp = 1/4.
    # Show
    Show_grafs = False
    format_fig = 'png'

    for i in range(0,N):
        for j in range(0,N):
            Sip_xDy[i][j] = d_GIMP(L,lp,xp[i],xi) * GIMP(L,lp,yp[j],yi) 
            
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
        
    surf = ax.plot_surface(Xp, Yp, Sip_xDy, rstride=1,
                           cstride=1, cmap=cm.jet,
                           linewidth=0, antialiased=False)
    ax.set_zlim3d(-1.00, 1.00)
    plt.tight_layout()
    plt.savefig('GIMP_Shape_Fun_dx.%s'%(format_fig))
    if Show_grafs:
        plt.show()
    fig.clear()

def print_dy_GIMP():
    N = 50
    xp = yp = np.linspace(-2, 2, N)
    Xp, Yp = np.meshgrid(xp, yp, indexing='ij', sparse=False)
    Sip_xDy = np.zeros([N,N])
    xi = yi = 0
    L = 1.
    lp = 1/4.
    # Show
    Show_grafs = False
    format_fig = 'png'

    for i in range(0,N):
        for j in range(0,N):
            Sip_xDy[i][j] = GIMP(L,lp,xp[i],xi) * d_GIMP(L,lp,yp[j],yi) 
            

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
            
    surf = ax.plot_surface(Xp, Yp, Sip_xDy, rstride=1,
                           cstride=1, cmap=cm.jet,
                           linewidth=0, antialiased=False)
    ax.set_zlim3d(-1.00, 1.00)
    plt.tight_layout()
    plt.savefig('GIMP_Shape_Fun_dy.%s'%(format_fig))  
    if Show_grafs:
        plt.show()
    fig.clear()
    
def print_GIMP():
    N = 50
    xp = yp = np.linspace(-2, 2, N)
    Xp, Yp = np.meshgrid(xp, yp, indexing='ij', sparse=False)
    Sip_xy = np.zeros([N,N])
    xi = yi = 0
    L = 1.
    lp = 1/4.
    # Show
    Show_grafs = False
    format_fig = 'png'

    for i in range(0,N):
        for j in range(0,N):
            Sip_xy[i][j] = GIMP(L,lp,xp[i],xi) * GIMP(L,lp,yp[j],yi)
            
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
    plt.savefig('GIMP_Shape_Fun.%s'%(format_fig))  
    if Show_grafs:
        plt.show()
    fig.clear()

print_GIMP()
print_dx_GIMP()
print_dy_GIMP()

