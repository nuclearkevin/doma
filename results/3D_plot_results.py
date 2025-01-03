import numpy as np
import matplotlib.pyplot as plt

from matplotlib import cm
import matplotlib.colors as colors

dim_x = 80
dim_y = 80
dim_z = 80
z_point = 40
#z_point = int(dim_z / 2)

diamond_raw = np.loadtxt('Kobayashi_1_diamond_flux.txt')
diamond_3D = diamond_raw.reshape((dim_x, dim_y, dim_z))
plt.pcolormesh(diamond_3D[z_point], norm=colors.LogNorm(vmin=diamond_raw.min(), vmax=diamond_raw.max()),
               cmap=cm.coolwarm, shading='auto')
plt.colorbar()
plt.savefig('diamond_difference.png')
plt.show()
