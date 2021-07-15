#!/opt/anaconda3/bin/python3

import numpy as np
import matplotlib.pyplot as plt

strain_xy = np.zeros(4000) 
strain_yx = np.zeros(4000)
strain_xx = np.ones(4000)
strain_yy = np.ones(4000)

strain_xy[0:1000] = np.linspace(0.0, 0.01, num=1000)
strain_xy[1000:2000] = np.linspace(0.01, 0.0, num=1000)
strain_xy[2000:3000] = np.linspace(0.0, -0.01, num=1000)
strain_xy[3000:4000] = np.linspace(-0.01, 0.0, num=1000)

strain = np.array([strain_xx, strain_xy, strain_yx, strain_yy]).transpose()

plt.plot(strain_xy)
plt.show()

np.savetxt("Strains.txt", strain,
        header="NROWS=%i NCOLS=%i PARSER=%s"%(4000,4,"%d,%d,%d,%d"), delimiter=",")

