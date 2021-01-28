import numpy as np
import matplotlib.pyplot as plt

w = -2e6
DetalT = 0.000759 # sg
Nsteps = 500
Tend = DetalT*Nsteps

Ramp = np.zeros(Nsteps)
T = np.linspace(0,Tend,Nsteps)

for i in range(0,Nsteps):
    if T[i] < 0.1:
        Ramp[i] = w*(min(T[i],0.1)/0.1)
    else:
        Ramp[i] = w


plt.plot(T,Ramp)
plt.show()

np.savetxt('CurveRamp.txt',Ramp,fmt='%2.5f',header='DAT_CURVE NUM#%i \n CUSTOM_CURVE'%(Nsteps) )

# DAT_CURVE NUM#632 
# CUSTOM_CURVE
