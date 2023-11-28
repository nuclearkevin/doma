import numpy as np
import matplotlib.pyplot as plt

from matplotlib import cm
import matplotlib.colors as colors

import matplotlib.animation as animation

dim_x = 80
dim_y = 80
dim_z = 80

blocks_raw = np.loadtxt('Kobayashi_1_diamond_blocks.txt')
blocks_3D = blocks_raw.reshape((dim_x, dim_y, dim_z))

diamond_raw = np.loadtxt('Kobayashi_1_diamond_flux.txt')
diamond_3D = diamond_raw.reshape((dim_x, dim_y, dim_z))
step_raw = np.loadtxt('Kobayashi_1_step_flux.txt')
step_3D = step_raw.reshape((dim_x, dim_y, dim_z))

fig, ax = plt.subplots()
artists = []
for z in range(dim_z):
  container = ax.imshow(blocks_3D[z], cmap=cm.Spectral, vmin=0.0, vmax=2.0)
  title = ax.text(2.0, 4.0, 'z = {}'.format(z))
  if z == 0:
    ax.imshow(blocks_3D[z], cmap=cm.Spectral, vmin=0.0, vmax=2.0)
  artists.append([container, title])

fig.colorbar(artists[40][0])
ani = animation.ArtistAnimation(fig, artists, interval=100, blit=True)
ani.save(filename="./varying_z/blocks_variable_z.gif", writer="pillow", dpi=500)
plt.show()

quit()

fig, ax = plt.subplots()
artists = []
for z in range(dim_z):
  container = ax.imshow(diamond_3D[z], norm=colors.LogNorm(vmin=np.abs(1e-5), vmax=diamond_raw.max()),
                        cmap=cm.coolwarm)
  title = ax.text(2.0, 4.0, 'z = {}'.format(z))
  if z == 0:
    ax.imshow(diamond_3D[z], norm=colors.LogNorm(vmin=1e-5, vmax=diamond_raw.max()),
              cmap=cm.coolwarm)
  artists.append([container, title])

ani = animation.ArtistAnimation(fig, artists, interval=100, blit=True)
ani.save(filename="./varying_z/dd_variable_z.gif", writer="pillow", dpi=500)
plt.show()

fig, ax = plt.subplots()
artists = []
for z in range(dim_z):
  container = ax.imshow(step_3D[z], norm=colors.LogNorm(vmin=step_3D.min(), vmax=step_3D.max()),
                        cmap=cm.coolwarm)
  title = ax.text(2.0, 4.0, 'z = {}'.format(z))
  if z == 0:
    ax.imshow(step_3D[z], norm=colors.LogNorm(vmin=np.abs(1e-5), vmax=diamond_raw.max()),
              cmap=cm.coolwarm)
  artists.append([container, title])

ani = animation.ArtistAnimation(fig, artists, interval=100, blit=True)
ani.save(filename="./varying_z/sc_variable_z.gif", writer="pillow", dpi=500)
plt.show()
