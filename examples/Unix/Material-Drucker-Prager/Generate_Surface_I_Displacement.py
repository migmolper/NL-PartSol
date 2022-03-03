import numpy as np
import matplotlib.pyplot as plt

uf = 0.1
tf = 12.0
tc = 6.0
CFL = 0.1
h = 0.833
CEL = 223.6
DT = CFL*h/CEL
Nsteps = (int)(tf/DT)

print(Nsteps)

Displacements = np.ones(Nsteps)*(-2.23533064266368e-05)
T = np.linspace(0,tf,Nsteps)

for i in range(0,Nsteps):
    Displacements[i] -= 2*DT*(uf/tf**2)*(T[i]-tc)*(T[i]>tc)

plt.plot(T,Displacements)
plt.show()

#np.savetxt('Curve_Surface_I.txt',Displacements,fmt='%2.15f',header='DAT_CURVE NUM#%i \n CUSTOM_CURVE'%(Nsteps) )

# DAT_CURVE NUM#632 
# CUSTOM_CURVE
