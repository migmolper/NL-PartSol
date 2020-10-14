import numpy as np
import matplotlib.pyplot as plt

g = -9.81
DetalT = 0.01 # sg
Nsteps = 1000
Tend = DetalT*Nsteps

Gravity = np.zeros(Nsteps)
T = np.linspace(0,Tend,Nsteps)

for i in range(0,Nsteps):
    if T[i] < 1.0:
        Gravity[i] = g*min(T[i],1)
    else:
        Gravity[i] = g


plt.plot(T,Gravity)
plt.show()

np.savetxt('Gravity_Load.txt',Gravity,fmt='%2.5f',header='DAT_CURVE NUM#%i \n CUSTOM_CURVE'%(Nsteps) )

# DAT_CURVE NUM#632 
# CUSTOM_CURVE
