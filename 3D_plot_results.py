#!/usr/bin/env python3
import argparse as ap
import os

import numpy as np
import matplotlib.pyplot as plt

from matplotlib import cm
import matplotlib.colors as colors

def main():
  # Setup command line arguements.
  parser = ap.ArgumentParser(description='Plots the scalar flux generated by DOMA for 3D transport problems.')
  parser.add_argument('-i', type=str, dest='file_name', required=True, help='The path of the input file (assuming the outputs are in the same directory).')
  parser.add_argument('-z', type=int, dest='z_slice', required=True, help='The z-slice of the geometry that should be plotted.')
  parser.add_argument('-l', '--log-scale', action='store_true', dest='log_scale', help='Whether the plot is log scale or not..')
  parser.add_argument('-s', type=int, dest='step', required=False, default=-1, help='The timestep to plot. Negative timesteps indicate a steady state solve.')

  cli_args = parser.parse_args()
  dir_path = os.path.dirname(os.path.realpath(cli_args.file_name))
  input_name = os.path.splitext(os.path.basename(cli_args.file_name))[0]

  dim_file = open(dir_path + "/" + input_name + "_dims.txt", "r")
  dim_x = int(str(dim_file.readline()).replace("num_x: ", "").replace("\n", ""))
  dim_y = int(str(dim_file.readline()).replace("num_y: ", "").replace("\n", ""))
  dim_z = int(str(dim_file.readline()).replace("num_z: ", "").replace("\n", ""))
  grps = int(str(dim_file.readline()).replace("num_g: ", "").replace("\n", ""))

  x_vals = np.loadtxt(dir_path + "/" + input_name + "_meshx.txt")
  y_vals = np.loadtxt(dir_path + "/" + input_name + "_meshy.txt")
  x_3D = x_vals.reshape((dim_x, dim_y, dim_z))
  y_3D = y_vals.reshape((dim_x, dim_y, dim_z))

  s = int(cli_args.z_slice)
  if s < 0 or s >= dim_z:
    print(s, "is an invalid slice! Please ensure 0 <= s <", dim_z)
    quit()

  for grp in range(grps):
    raw_flux = np.array([1])
    if cli_args.step >= 0:
      raw_flux = np.loadtxt(dir_path + "/" + input_name + "_t" + str(cli_args.step) + "_g" + str(grp) + "_flux.txt")
    else:
      raw_flux = np.loadtxt(dir_path + "/" + input_name + "_g" + str(grp) + "_flux.txt")

    flux_3D = raw_flux.reshape((dim_x, dim_y, dim_z))
    if cli_args.log_scale == True:
      fig, ax = plt.subplots()
      mappable = ax.pcolor(x_3D[s], y_3D[s], flux_3D[s],
                          norm=colors.LogNorm(vmin=raw_flux.min(), vmax=raw_flux.max()),
                          cmap=cm.coolwarm, shading='auto')
      cbar = fig.colorbar(mappable)
      cbar.ax.set_ylabel('Group ' + str(grp) + ' Scalar Flux (s$^{-1}$ cm$^{-2}$)')
      ax.set_xlabel('x (cm)')
      ax.set_ylabel('y (cm)')

      if cli_args.step > 0:
        plt.savefig(dir_path + "/" + input_name + "_t" + str(cli_args.step) + "_g" + str(grp) + "_z" + str(s) + "_flux.png", format='png')
      else:
        plt.savefig(dir_path + "/" + input_name + "_g" + str(grp) + "_z" + str(s) + "_flux.png", format='png')
      plt.show()
      plt.close()
    else:
      fig, ax = plt.subplots()
      mappable = ax.pcolor(x_3D[s], y_3D[s], flux_3D[s],
                          cmap=cm.coolwarm, shading='auto')
      cbar = fig.colorbar(mappable)
      cbar.ax.set_ylabel('Group ' + str(grp) + ' Scalar Flux (s$^{-1}$ cm$^{-2}$)')
      ax.set_xlabel('x (cm)')
      ax.set_ylabel('y (cm)')

      if cli_args.step > 0:
        plt.savefig(dir_path + "/" + input_name + "_t" + str(cli_args.step) + "_g" + str(grp) + "_z" + str(s) + "_flux.png", format='png')
      else:
        plt.savefig(dir_path + "/" + input_name + "_g" + str(grp) + "_z" + str(s) + "_flux.png", format='png')
      plt.show()
      plt.close()

if __name__ == "__main__":
    main()