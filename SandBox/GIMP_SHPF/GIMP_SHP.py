#!/opt/anaconda3/bin/ python

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

N = 50
xp = yp = np.linspace(-1.2, 1.2, N)
Xp, Yp = np.meshgrid(xp, yp, indexing='ij', sparse=False)
Sip_xDy = np.zeros([N,N])
xi = yi = 0
L = 1.
lp = 1/4.

for i in range(0,N):
    for j in range(0,N):
        Sip_xDy[i][j] = GIMP(L,lp,xp[i],xi) * d_GIMP(L,lp,yp[j],yi) 

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(Xp, Yp, Sip_xDy, rstride=1, cstride=1, cmap=cm.jet,
        linewidth=0, antialiased=False)
ax.set_zlim3d(-1.00, 1.00)

# ax.w_zaxis.set_major_locator(LinearLocator(10))
# ax.w_zaxis.set_major_formatter(FormatStrFormatter('%.03f'))

plt.show()
