#!/opt/anaconda3/bin/python

import pandas as pd

# Nodes
data_left_wall   = pd.read_csv("left_GID.txt", header=None, sep="\s+")
data_right_wall  = pd.read_csv("right_GID.txt", header=None, sep="\s+")
data_bottom_wall = pd.read_csv("bottom_GID.txt", header=None, sep="\s+")
data_rollers     = pd.read_csv("rollers_GID.txt", header=None, sep="\s+")

(data_left_wall[0]-1).to_csv("../Left_wall.txt", header=None, index=False)
(data_right_wall[0]-1).to_csv("../Right_wall.txt", header=None, index=False)
(data_bottom_wall[0]-1).to_csv("../Bottom_wall.txt", header=None, index=False)
(data_rollers[0]-1).to_csv("../Rollers.txt", header=None, index=False)

# Particles
data_footing = pd.read_csv("footing_GID.txt", header=None, sep="\s+")
data_soil    = pd.read_csv("soil_GID.txt", header=None, sep="\s+")

(data_footing[0]-1).to_csv("../Footing.txt", header=None, index=False)
(data_soil[0]-1).to_csv("../Soil.txt", header=None, index=False)
