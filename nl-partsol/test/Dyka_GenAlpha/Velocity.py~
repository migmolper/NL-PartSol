#!/opt/anaconda3/bin/python

import glob,os
import numpy as np
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as VN
import matplotlib.pyplot as plt

CEL = 1
CFL = 0.7

fig, ax = plt.subplots()

Route_Results = '%s/%s/%s/%s*.vtk'%('test','Dyka_MPM','Resultados','MPM_MPM_VALUES_')
Files = sorted(glob.glob(Route_Results),key=os.path.getmtime)
Num_Files = len(Files)
DeltaX = 1
DeltaT = CFL/CEL*DeltaX;
t = np.linspace(0.0, Num_Files*DeltaT, Num_Files)
Velocity_GP = np.zeros((Num_Files,1))

# Read files
for i in range(0,Num_Files-1):
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(Files[i])
    reader.ReadAllTensorsOn()
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()
    V = VN.vtk_to_numpy(data.GetCellData().GetArray('VELOCITY'))
    i_GP = 9
    Velocity_GP[i] = V[i_GP,0]

ax.plot(Velocity_GP)    
plt.grid()
plt.show()
