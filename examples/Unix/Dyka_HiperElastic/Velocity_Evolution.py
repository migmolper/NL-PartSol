#!/opt/anaconda3/bin/python3

import pandas as pd
import matplotlib.pyplot as plt
import glob,os
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as VN


plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

control_particles_data = pd.read_csv("control_particles.txt",header=None)

Rute_Consistent = '%s/%s*.vtk' %('Resultados_consistent','Particles_')
Files_Consistent = sorted(glob.glob(Rute_Consistent),key=os.path.getmtime)
Num_Files_Consistent = len(Files_Consistent)

#T_step = 60
T_step = 100

reader = vtkUnstructuredGridReader()
reader.SetFileName(Files_Consistent[T_step])
reader.ReadAllTensorsOn()
reader.ReadAllVectorsOn()
reader.ReadAllScalarsOn()
reader.Update()
data = reader.GetOutput()
V_vtk = VN.vtk_to_numpy(data.GetCellData().GetArray('VELOCITY'))

plt.plot(V_vtk[control_particles_data[0],0])
plt.show()
