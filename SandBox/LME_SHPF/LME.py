import numpy as np
import matplotlib.pyplot as plt


def InOut(Coords_GP, Xmesh, Ymesh, Zmesh, h, gamma, TOL):
    """
    Prototipe of the search algorithm for the
    LME shape functions of M.Ortiz and M.Arroyo
    """
    DIST_SEARCH = h*np.sqrt(-np.log(TOL)/gamma)
    SHAPE = np.shape(Xmesh)
    for i in range(0,SHAPE[0]):
        for j in range(0,SHAPE[1]):
            dist = np.sqrt((Xmesh[i,j] - Coords_GP[0])**2 + (Ymesh[i,j] - Coords_GP[1])**2)
            if dist < DIST_SEARCH:
                Zmesh[i,j] = 0

    return DIST_SEARCH

x = np.linspace(0, 1, 10)
y = np.linspace(0, 1, 10)
h = 1/10. # Nodal spacing
TOL = 10**(-6) # Tolerance
gamma = 6 # Shape parameter
Xmesh, Ymesh = np.meshgrid(x, y)
Zmesh = np.ones_like(Xmesh)
Coords_GP = np.array([0.5,0.5])
   
DIST_SEARCH = InOut(Coords_GP, Xmesh, Ymesh, Zmesh, h, gamma, TOL)
circle = plt.Circle((0.5, 0.5), DIST_SEARCH, color='b', fill=False)

print(DIST_SEARCH)

fig, ax = plt.subplots()
ax.add_artist(circle)
plt.scatter(Xmesh,Ymesh,c = Zmesh)
ax.axis('equal')
plt.show()
