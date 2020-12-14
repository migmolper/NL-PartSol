#!/opt/anaconda3/bin/python3

import numpy as np

strain_xy = np.linspace(0, 1, num=4000) 
strain_yx = np.zeros(4000)
strain_xx = np.ones(4000)
strain_yy = np.ones(4000)

strain = np.array([strain_xx, strain_xy, strain_yx, strain_yy]).transpose()

np.savetxt("Strains.txt", strain,
        header="NROWS=%i NCOLS=%i PARSER=%s"%(4000,2,"%d,%d,%d,%d"), delimiter=",")

