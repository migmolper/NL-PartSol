import numpy as np
import matplotlib.pyplot as plt

g = - 10
tc = 6.0
tf = 8.0
CFL = 0.1
h = 0.833
CEL = 223.6
DT = CFL*h/CEL
Nsteps = (int)(tf/DT)

print(Nsteps)

Load = np.zeros(Nsteps)
T = np.linspace(0,tf,Nsteps)

for i in range(0,Nsteps):
        Load[i] = min(1,T[i]/tc)*g

plt.plot(T,Load)
plt.show()

np.savetxt('Curve_Gravity.txt',Load,fmt='%2.15f',header='DAT_CURVE NUM#%i \n CUSTOM_CURVE'%(Nsteps) )

# DAT_CURVE NUM#632 
# CUSTOM_CURVE
