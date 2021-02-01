import numpy as np
import matplotlib.pyplot as plt

w = -40e3
DetalT = 0.001710 # sg
Nsteps = 500
Tend = DetalT*Nsteps

Ramp = np.zeros(Nsteps)
T = np.linspace(0,Tend,Nsteps)
T_c = 0.6

for i in range(0,Nsteps):
    if T[i] < T_c:
        Ramp[i] = w*(min(T[i],T_c)/T_c)
    else:
        Ramp[i] = w


plt.plot(T,Ramp)
plt.show()

np.savetxt('CurveRamp.txt',Ramp,fmt='%2.5f',header='DAT_CURVE NUM#%i \n CUSTOM_CURVE'%(Nsteps) )

# DAT_CURVE NUM#632 
# CUSTOM_CURVE
