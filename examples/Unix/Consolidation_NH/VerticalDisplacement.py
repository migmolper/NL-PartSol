import glob,os
import numpy as np
import pandas as pd
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as VN

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
# from matplotlib import rc
plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
plt.rc('text', usetex=True)

####################### Numerical parameters ##########################

rho = 6.00E+03
E = 5E6
CEL = np.sqrt(E/rho)
DeltaX = 2
T_END = 2E-4
RES = 'Resultados'
CFL = 0.1
TIS = ['FE','PCE']
I_GP = 97
i_RES = 1

######################## Compare velocity #############################

# Total plot
fig, ax = plt.subplots()
ax.set_xlabel('Time (sg)', fontsize=16)
ax.set_ylabel('Displacement (m)', fontsize=16)
# ax.set_xlim(0.00, 0.0002)
# ax.set_xticks([0, 0.00005 ,0.0001, 0.00015, 0.0002])
# ax.set_yticks([-5.0, -2.5, 0.0, 2.5 , 5.0])
ax.tick_params(axis='both', which='major', labelsize=15)
ax.grid(axis='both')


FEM_res = pd.read_csv("FEM/GeHoMadrid.txt","\s+",header=None)
T = FEM_res[0].values
D = FEM_res[1].values
ax.plot(T,D,label='FEM')

RES_MPM = 'Resultados_LME_05_PCE/%s*.vtk' %('MPM_MPM_VALUES_')
FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
NUM_FILES = len(FILES_MPM)
DELTA_T = DeltaX*CFL/CEL        
yDIS_MPM = np.zeros(NUM_FILES)
T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# Read files
for i in range(0,NUM_FILES):
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(FILES_MPM[i])
    reader.ReadAllVectorsOn()
    reader.Update()
    data = reader.GetOutput()
    Point_cordinates = data.GetPoints().GetData()
    D_vtk = VN.vtk_to_numpy(Point_cordinates)
    yDIS_MPM[i] = D_vtk[I_GP,1]    
ax.plot(T_MPM,yDIS_MPM-9.577350,label='LME-$\gamma=0.8$ PCE')


# RES_MPM = 'Resultados_LME_1_PCE/%s*.vtk' %('MPM_MPM_VALUES_')
# FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
# NUM_FILES = len(FILES_MPM)
# DELTA_T = DeltaX*CFL/CEL        
# yDIS_MPM = np.zeros(NUM_FILES)
# T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# # Read files
# for i in range(0,NUM_FILES):
#     reader = vtkUnstructuredGridReader()
#     reader.SetFileName(FILES_MPM[i])
#     reader.ReadAllVectorsOn()
#     reader.Update()
#     data = reader.GetOutput()
#     Point_cordinates = data.GetPoints().GetData()
#     D_vtk = VN.vtk_to_numpy(Point_cordinates)
#     yDIS_MPM[i] = D_vtk[I_GP,1]    
# ax.plot(T_MPM,yDIS_MPM-9.577350,label='LME-$\gamma=1$ PCE')

# RES_MPM = 'Resultados_LME_1_FE/%s*.vtk' %('MPM_MPM_VALUES_')
# FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
# NUM_FILES = len(FILES_MPM)
# DELTA_T = DeltaX*CFL/CEL        
# yDIS_MPM = np.zeros(NUM_FILES)
# T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# # Read files
# for i in range(0,NUM_FILES):
#     reader = vtkUnstructuredGridReader()
#     reader.SetFileName(FILES_MPM[i])
#     reader.ReadAllVectorsOn()
#     reader.Update()
#     data = reader.GetOutput()
#     Point_cordinates = data.GetPoints().GetData()
#     D_vtk = VN.vtk_to_numpy(Point_cordinates)
#     yDIS_MPM[i] = D_vtk[I_GP,1]    
# ax.plot(T_MPM,yDIS_MPM-9.577350,label='LME-$\gamma=1$ FE')

# RES_MPM = 'Resultados_LME_2_PCE/%s*.vtk' %('MPM_MPM_VALUES_')
# FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
# NUM_FILES = len(FILES_MPM)
# DELTA_T = DeltaX*CFL/CEL        
# yDIS_MPM = np.zeros(NUM_FILES)
# T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# # Read files
# for i in range(0,NUM_FILES):
#     reader = vtkUnstructuredGridReader()
#     reader.SetFileName(FILES_MPM[i])
#     reader.ReadAllVectorsOn()
#     reader.Update()
#     data = reader.GetOutput()
#     Point_cordinates = data.GetPoints().GetData()
#     D_vtk = VN.vtk_to_numpy(Point_cordinates)
#     yDIS_MPM[i] = D_vtk[I_GP,1]    
# ax.plot(T_MPM,yDIS_MPM-9.577350,label='LME-$\gamma=2$ PCE')

# RES_MPM = 'Resultados_LME_2_FE/%s*.vtk' %('MPM_MPM_VALUES_')
# FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
# NUM_FILES = len(FILES_MPM)
# DELTA_T = DeltaX*CFL/CEL        
# yDIS_MPM = np.zeros(NUM_FILES)
# T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# # Read files
# for i in range(0,NUM_FILES):
#     reader = vtkUnstructuredGridReader()
#     reader.SetFileName(FILES_MPM[i])
#     reader.ReadAllVectorsOn()
#     reader.Update()
#     data = reader.GetOutput()
#     Point_cordinates = data.GetPoints().GetData()
#     D_vtk = VN.vtk_to_numpy(Point_cordinates)
#     yDIS_MPM[i] = D_vtk[I_GP,1]    
# ax.plot(T_MPM,yDIS_MPM-9.577350,label='LME-$\gamma=2$ FE')


RES_MPM = 'Resultados_LME_3_PCE/%s*.vtk' %('MPM_MPM_VALUES_')
FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
NUM_FILES = len(FILES_MPM)
DELTA_T = DeltaX*CFL/CEL        
yDIS_MPM = np.zeros(NUM_FILES)
T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# Read files
for i in range(0,NUM_FILES):
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(FILES_MPM[i])
    reader.ReadAllVectorsOn()
    reader.Update()
    data = reader.GetOutput()
    Point_cordinates = data.GetPoints().GetData()
    D_vtk = VN.vtk_to_numpy(Point_cordinates)
    yDIS_MPM[i] = D_vtk[I_GP,1]    
ax.plot(T_MPM,yDIS_MPM-9.577350,label='LME-$\gamma=3.0$ PCE')

# RES_MPM = 'Resultados_LME_3_FE/%s*.vtk' %('MPM_MPM_VALUES_')
# FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
# NUM_FILES = len(FILES_MPM)
# DELTA_T = DeltaX*CFL/CEL        
# yDIS_MPM = np.zeros(NUM_FILES)
# T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# # Read files
# for i in range(0,NUM_FILES):
#     reader = vtkUnstructuredGridReader()
#     reader.SetFileName(FILES_MPM[i])
#     reader.ReadAllVectorsOn()
#     reader.Update()
#     data = reader.GetOutput()
#     Point_cordinates = data.GetPoints().GetData()
#     D_vtk = VN.vtk_to_numpy(Point_cordinates)
#     yDIS_MPM[i] = D_vtk[I_GP,1]    
# ax.plot(T_MPM,yDIS_MPM-9.577350,label='LME-$\gamma=3$ FE')



RES_MPM = 'Resultados_Q4_PCE/%s*.vtk' %('MPM_MPM_VALUES_')
FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
NUM_FILES = len(FILES_MPM)
DELTA_T = DeltaX*CFL/CEL        
yDIS_MPM = np.zeros(NUM_FILES)
T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# Read files
for i in range(0,NUM_FILES):
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(FILES_MPM[i])
    reader.ReadAllVectorsOn()
    reader.Update()
    data = reader.GetOutput()
    Point_cordinates = data.GetPoints().GetData()
    D_vtk = VN.vtk_to_numpy(Point_cordinates)
    yDIS_MPM[i] = D_vtk[I_GP,1]    
ax.plot(T_MPM,yDIS_MPM-9.577350,label='Q4 PCE')



# RES_MPM = 'Resultados_Q4_FE/%s*.vtk' %('MPM_MPM_VALUES_')
# FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
# NUM_FILES = len(FILES_MPM)
# DELTA_T = DeltaX*CFL/CEL        
# yDIS_MPM = np.zeros(NUM_FILES)
# T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# # Read files
# for i in range(0,NUM_FILES):
#     reader = vtkUnstructuredGridReader()
#     reader.SetFileName(FILES_MPM[i])
#     reader.ReadAllVectorsOn()
#     reader.Update()
#     data = reader.GetOutput()
#     Point_cordinates = data.GetPoints().GetData()
#     D_vtk = VN.vtk_to_numpy(Point_cordinates)
#     yDIS_MPM[i] = D_vtk[I_GP,1]    
# ax.plot(T_MPM,yDIS_MPM-9.577350,label='Q4 FE')



# RES_MPM = 'Resultados_GIMP_FE/%s*.vtk' %('MPM_MPM_VALUES_')
# FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
# NUM_FILES = len(FILES_MPM)
# DELTA_T = DeltaX*CFL/CEL        
# yDIS_MPM = np.zeros(NUM_FILES)
# T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# # Read files
# for i in range(0,NUM_FILES):
#     reader = vtkUnstructuredGridReader()
#     reader.SetFileName(FILES_MPM[i])
#     reader.ReadAllVectorsOn()
#     reader.Update()
#     data = reader.GetOutput()
#     Point_cordinates = data.GetPoints().GetData()
#     D_vtk = VN.vtk_to_numpy(Point_cordinates)
#     yDIS_MPM[i] = D_vtk[I_GP,1]    
# ax.plot(T_MPM,yDIS_MPM-9.5,label='uGIMP FE')

RES_MPM = 'Resultados_GIMP_PCE/%s*.vtk' %('MPM_MPM_VALUES_')
FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
NUM_FILES = len(FILES_MPM)
DELTA_T = DeltaX*CFL/CEL        
yDIS_MPM = np.zeros(NUM_FILES)
T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# Read files
for i in range(0,NUM_FILES):
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(FILES_MPM[i])
    reader.ReadAllVectorsOn()
    reader.Update()
    data = reader.GetOutput()
    Point_cordinates = data.GetPoints().GetData()
    D_vtk = VN.vtk_to_numpy(Point_cordinates)
    yDIS_MPM[i] = D_vtk[I_GP,1]    
ax.plot(T_MPM,yDIS_MPM-9.5,label='uGIMP PCE')

plt.legend(loc='upper right',markerscale=10)

# Zoom
axins = zoomed_inset_axes(ax, 1.8, 'upper center')
x1, x2, y1, y2 = 10, 15, -0.6, -0.5
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
# axins.set_yticks([-1.0, -0.5, 0.0, 0.5 , 1.0])
axins.tick_params(axis = 'both', which = 'both', labelsize = 10)
plt.yticks(visible=True)
plt.xticks(visible=False)
mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.1")
# axins.plot(T_MPM_PCE,VEL_MPM_FE/1000,markersize=5,label='FE')
# axins.plot(T_MPM_PCE,VEL_MPM_PCE/1000,markersize=5,label='PCE')
# axins.plot(T_MPM_PCE,VEL_ANALY_INT/1000,markersize=5,label='Analytical')


FEM_res = pd.read_csv("FEM/GeHoMadrid.txt","\s+",header=None)
T = FEM_res[0].values
D = FEM_res[1].values
axins.plot(T,D,label='FEM')

RES_MPM = 'Resultados_LME_05_PCE/%s*.vtk' %('MPM_MPM_VALUES_')
FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
NUM_FILES = len(FILES_MPM)
DELTA_T = DeltaX*CFL/CEL        
yDIS_MPM = np.zeros(NUM_FILES)
T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# Read files
for i in range(0,NUM_FILES):
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(FILES_MPM[i])
    reader.ReadAllVectorsOn()
    reader.Update()
    data = reader.GetOutput()
    Point_cordinates = data.GetPoints().GetData()
    D_vtk = VN.vtk_to_numpy(Point_cordinates)
    yDIS_MPM[i] = D_vtk[I_GP,1]    
axins.plot(T_MPM,yDIS_MPM-9.577350,label='LME-$\gamma=0.8$ PCE')


# RES_MPM = 'Resultados_LME_1_PCE/%s*.vtk' %('MPM_MPM_VALUES_')
# FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
# NUM_FILES = len(FILES_MPM)
# DELTA_T = DeltaX*CFL/CEL        
# yDIS_MPM = np.zeros(NUM_FILES)
# T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# # Read files
# for i in range(0,NUM_FILES):
#     reader = vtkUnstructuredGridReader()
#     reader.SetFileName(FILES_MPM[i])
#     reader.ReadAllVectorsOn()
#     reader.Update()
#     data = reader.GetOutput()
#     Point_cordinates = data.GetPoints().GetData()
#     D_vtk = VN.vtk_to_numpy(Point_cordinates)
#     yDIS_MPM[i] = D_vtk[I_GP,1]    
# axins.plot(T_MPM,yDIS_MPM-9.577350,label='LME-$\gamma=1$ PCE')

# RES_MPM = 'Resultados_LME_1_FE/%s*.vtk' %('MPM_MPM_VALUES_')
# FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
# NUM_FILES = len(FILES_MPM)
# DELTA_T = DeltaX*CFL/CEL        
# yDIS_MPM = np.zeros(NUM_FILES)
# T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# # Read files
# for i in range(0,NUM_FILES):
#     reader = vtkUnstructuredGridReader()
#     reader.SetFileName(FILES_MPM[i])
#     reader.ReadAllVectorsOn()
#     reader.Update()
#     data = reader.GetOutput()
#     Point_cordinates = data.GetPoints().GetData()
#     D_vtk = VN.vtk_to_numpy(Point_cordinates)
#     yDIS_MPM[i] = D_vtk[I_GP,1]    
# axins.plot(T_MPM,yDIS_MPM-9.577350,label='LME-$\gamma=1$ FE')

# RES_MPM = 'Resultados_LME_2_PCE/%s*.vtk' %('MPM_MPM_VALUES_')
# FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
# NUM_FILES = len(FILES_MPM)
# DELTA_T = DeltaX*CFL/CEL        
# yDIS_MPM = np.zeros(NUM_FILES)
# T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# # Read files
# for i in range(0,NUM_FILES):
#     reader = vtkUnstructuredGridReader()
#     reader.SetFileName(FILES_MPM[i])
#     reader.ReadAllVectorsOn()
#     reader.Update()
#     data = reader.GetOutput()
#     Point_cordinates = data.GetPoints().GetData()
#     D_vtk = VN.vtk_to_numpy(Point_cordinates)
#     yDIS_MPM[i] = D_vtk[I_GP,1]    
# axins.plot(T_MPM,yDIS_MPM-9.577350,label='LME-$\gamma=2$ PCE')

# RES_MPM = 'Resultados_LME_2_FE/%s*.vtk' %('MPM_MPM_VALUES_')
# FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
# NUM_FILES = len(FILES_MPM)
# DELTA_T = DeltaX*CFL/CEL        
# yDIS_MPM = np.zeros(NUM_FILES)
# T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# # Read files
# for i in range(0,NUM_FILES):
#     reader = vtkUnstructuredGridReader()
#     reader.SetFileName(FILES_MPM[i])
#     reader.ReadAllVectorsOn()
#     reader.Update()
#     data = reader.GetOutput()
#     Point_cordinates = data.GetPoints().GetData()
#     D_vtk = VN.vtk_to_numpy(Point_cordinates)
#     yDIS_MPM[i] = D_vtk[I_GP,1]    
# axins.plot(T_MPM,yDIS_MPM-9.577350,label='LME-$\gamma=2$ FE')


RES_MPM = 'Resultados_LME_3_PCE/%s*.vtk' %('MPM_MPM_VALUES_')
FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
NUM_FILES = len(FILES_MPM)
DELTA_T = DeltaX*CFL/CEL        
yDIS_MPM = np.zeros(NUM_FILES)
T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# Read files
for i in range(0,NUM_FILES):
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(FILES_MPM[i])
    reader.ReadAllVectorsOn()
    reader.Update()
    data = reader.GetOutput()
    Point_cordinates = data.GetPoints().GetData()
    D_vtk = VN.vtk_to_numpy(Point_cordinates)
    yDIS_MPM[i] = D_vtk[I_GP,1]    
axins.plot(T_MPM,yDIS_MPM-9.577350,label='LME-$\gamma=3.0$ PCE')

# RES_MPM = 'Resultados_LME_3_FE/%s*.vtk' %('MPM_MPM_VALUES_')
# FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
# NUM_FILES = len(FILES_MPM)
# DELTA_T = DeltaX*CFL/CEL        
# yDIS_MPM = np.zeros(NUM_FILES)
# T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# # Read files
# for i in range(0,NUM_FILES):
#     reader = vtkUnstructuredGridReader()
#     reader.SetFileName(FILES_MPM[i])
#     reader.ReadAllVectorsOn()
#     reader.Update()
#     data = reader.GetOutput()
#     Point_cordinates = data.GetPoints().GetData()
#     D_vtk = VN.vtk_to_numpy(Point_cordinates)
#     yDIS_MPM[i] = D_vtk[I_GP,1]    
# axins.plot(T_MPM,yDIS_MPM-9.577350,label='LME-$\gamma=3$ FE')



RES_MPM = 'Resultados_Q4_PCE/%s*.vtk' %('MPM_MPM_VALUES_')
FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
NUM_FILES = len(FILES_MPM)
DELTA_T = DeltaX*CFL/CEL        
yDIS_MPM = np.zeros(NUM_FILES)
T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# Read files
for i in range(0,NUM_FILES):
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(FILES_MPM[i])
    reader.ReadAllVectorsOn()
    reader.Update()
    data = reader.GetOutput()
    Point_cordinates = data.GetPoints().GetData()
    D_vtk = VN.vtk_to_numpy(Point_cordinates)
    yDIS_MPM[i] = D_vtk[I_GP,1]    
axins.plot(T_MPM,yDIS_MPM-9.577350,label='Q4 PCE')


# RES_MPM = 'Resultados_Q4_FE/%s*.vtk' %('MPM_MPM_VALUES_')
# FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
# NUM_FILES = len(FILES_MPM)
# DELTA_T = DeltaX*CFL/CEL        
# yDIS_MPM = np.zeros(NUM_FILES)
# T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# # Read files
# for i in range(0,NUM_FILES):
#     reader = vtkUnstructuredGridReader()
#     reader.SetFileName(FILES_MPM[i])
#     reader.ReadAllVectorsOn()
#     reader.Update()
#     data = reader.GetOutput()
#     Point_cordinates = data.GetPoints().GetData()
#     D_vtk = VN.vtk_to_numpy(Point_cordinates)
#     yDIS_MPM[i] = D_vtk[I_GP,1]    
# axins.plot(T_MPM,yDIS_MPM-9.577350,label='Q4 FE')



# RES_MPM = 'Resultados_GIMP_FE/%s*.vtk' %('MPM_MPM_VALUES_')
# FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
# NUM_FILES = len(FILES_MPM)
# DELTA_T = DeltaX*CFL/CEL        
# yDIS_MPM = np.zeros(NUM_FILES)
# T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# # Read files
# for i in range(0,NUM_FILES):
#     reader = vtkUnstructuredGridReader()
#     reader.SetFileName(FILES_MPM[i])
#     reader.ReadAllVectorsOn()
#     reader.Update()
#     data = reader.GetOutput()
#     Point_cordinates = data.GetPoints().GetData()
#     D_vtk = VN.vtk_to_numpy(Point_cordinates)
#     yDIS_MPM[i] = D_vtk[I_GP,1]    
# axins.plot(T_MPM,yDIS_MPM-9.5,label='uGIMP FE')

RES_MPM = 'Resultados_GIMP_PCE/%s*.vtk' %('MPM_MPM_VALUES_')
FILES_MPM = sorted(glob.glob(RES_MPM),key=os.path.getmtime)
NUM_FILES = len(FILES_MPM)
DELTA_T = DeltaX*CFL/CEL        
yDIS_MPM = np.zeros(NUM_FILES)
T_MPM = np.linspace(0.0,DELTA_T*NUM_FILES,NUM_FILES)
# Read files
for i in range(0,NUM_FILES):
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(FILES_MPM[i])
    reader.ReadAllVectorsOn()
    reader.Update()
    data = reader.GetOutput()
    Point_cordinates = data.GetPoints().GetData()
    D_vtk = VN.vtk_to_numpy(Point_cordinates)
    yDIS_MPM[i] = D_vtk[I_GP,1]    
axins.plot(T_MPM,yDIS_MPM-9.5,label='uGIMP PCE')

# plt.title('Left side')
plt.tight_layout()
plt.show()
