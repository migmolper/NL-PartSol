import numpy as np
import math
import matplotlib.pyplot as plt

g = -10
DeltaT = 6.025086e-04 # sg
Tend = 1E3 
Nsteps = math.floor(Tend/DeltaT)

Gravity = np.zeros(Nsteps)
T = np.linspace(0,Tend,Nsteps)
Tc = 100

for i in range(0,Nsteps):
    Gravity[i] = g*min(1,T[i]/Tc)


plt.plot(T,Gravity)
plt.show()

np.savetxt('Gravity_Load.txt',Gravity,fmt='%2.5f',header='DAT_CURVE NUM#%i \n CUSTOM_CURVE'%(Nsteps) )

# DAT_CURVE NUM#632 
# CUSTOM_CURVE
