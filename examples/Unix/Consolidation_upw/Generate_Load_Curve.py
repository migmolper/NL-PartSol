import numpy as np
import matplotlib.pyplot as plt

w = -4e2
DetalT = 0.000022 # sg
Nsteps = 500
Tend = DetalT*Nsteps

Ramp = np.zeros(Nsteps)
T = np.linspace(0,Tend,Nsteps)
Tc = 0.01

for i in range(0,Nsteps):
    if T[i] < Tc:
        Ramp[i] = w*(min(T[i],Tc)/Tc)
    else:
        Ramp[i] = w


plt.plot(T,Ramp)
plt.show()

np.savetxt('CurveRamp.txt',Ramp,fmt='%2.5f',header='DAT_CURVE NUM#%i \n CUSTOM_CURVE'%(Nsteps) )

# DAT_CURVE NUM#632 
# CUSTOM_CURVE