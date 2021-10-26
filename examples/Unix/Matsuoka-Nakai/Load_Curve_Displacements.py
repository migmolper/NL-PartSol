import numpy as np
import matplotlib.pyplot as plt

uf = 0.04
tf = 4.0
CFL = 0.1
h = 0.186484
CEL = 163.66
DT = CFL*h/CEL
Nsteps = (int)(tf/DT)

print(Nsteps)

Displacements = np.zeros(Nsteps)
T = np.linspace(0,tf,Nsteps)

for i in range(0,Nsteps):
    Displacements[i] = - 2*DT*(uf/tf**2)*T[i]

plt.plot(T,Displacements)
plt.show()

np.savetxt('Curve_Load.txt',Displacements,fmt='%2.15f',header='DAT_CURVE NUM#%i \n CUSTOM_CURVE'%(Nsteps) )

# DAT_CURVE NUM#632 
# CUSTOM_CURVE
