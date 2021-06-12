#!/opt/anaconda3/bin/python

import pandas as pd                                                                                                                                                                           
import matplotlib.pyplot as plt                                                                                                                                                               
import numpy as np

data = pd.read_csv("Particle_evolution_idx_0.csv",header=None)                                                                                                                                

center = 0.5*(data[0] + data[3])
radius = np.sqrt((0.5*(data[0] - data[3]))**2 + data[1]**2)

sigma_I = center + radius
sigma_II = center - radius

#plt.plot(sigma_I,sigma_II)
#plt.scatter(sigma_I,sigma_II)

#plt.scatter(data[5],np.sqrt(0.5*((data[0] - data[3])**2 + (data[3])**2
#    + data[0]**2 + 6*(data[1]**2))))

#plt.scatter(data[5],sigma_I)

plt.scatter(data[5],data[1])

#plt.scatter(0.5*(sigma_I + sigma_II),0.5*(sigma_I-sigma_II))

plt.show()   
