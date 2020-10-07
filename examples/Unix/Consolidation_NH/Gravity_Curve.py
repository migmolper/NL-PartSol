import numpy as np
import matplotlib.pyplot as plt

g = -9.81
E = 5E6 # Pa
rho = 6E3 # Kg/m^3 
CEL = np.sqrt(E/rho)
CFL = 0.1
DeltaX = 2 # m
DetalT = CFL*DeltaX/CEL # sg
Tend = 20 # sg
Nsteps = int(20/DetalT)

Gravity = np.zeros(Nsteps)
T = np.linspace(0,Tend,Nsteps)

for i in range(0,Nsteps):
    if T[i] < 10:
        Gravity[i] = 0.5*g*(np.sin(2*T[i]*np.pi/Tend - np.pi/2) + 1)
    else:
        Gravity[i] = g


plt.plot(T,Gravity)
plt.show()

np.savetxt('Gravity_Load.txt',Gravity,fmt='%2.5f',header='DAT_CURVE NUM#%i \n CUSTOM_CURVE'%(Nsteps) )

# DAT_CURVE NUM#632 
# CUSTOM_CURVE
