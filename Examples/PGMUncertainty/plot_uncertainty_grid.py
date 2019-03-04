import numpy as np
import matplotlib.pyplot as plt
import json
import scipy
from scipy.interpolate import griddata

def plot_uncertainty_grid(path):

    grid = np.loadtxt(path + "/grid.txt")
    uncertainty_percentage = np.loadtxt(path + "/uncertainty_over_reference_acc_percentage.txt")

    data = {}
    with open(path + "/input_file.json") as f:
        data = json.load(f)
    

    if data["PROJECTION_AXIS"] == 0:


        plt.imshow(uncertainty_percentage,interpolation = "bicubic")
        plt.colorbar()

        cont = plt.contour(uncertainty_percentage, levels = np.array([14,15,16,17]), origin='upper',cmap = plt.get_cmap("rainbow"))
        plt.gca().clabel(cont,cont.levels)
        # plt.xlabel("Y")
        # plt.ylabel("Z")
        # plt.title(r"Acceleration error $\left(\frac{\sqrt{\mathrm{trace}\left(P_{\mathbf{a}}\right)}}{\Vert \mathbf{a} \Vert }\cdot 100\ \%\right)$")
        plt.show()

       

  #       plt.figure()
        # CS = plt.contour(grid_y, grid_z, Z)
        # plt.clabel(CS, inline=1, fontsize=10)
        # plt.title('Simplest default with labels')



    elif data["PROJECTION_AXIS"] == 1:


        
        x_min = np.amin(grid[0,:])
        x_max = np.amax(grid[0,:])

        z_min = np.amin(grid[2,:])
        z_max = np.amax(grid[2,:])

        grid_x, grid_z = np.mgrid[x_min:x_max:1000j, z_min:z_max:1000j]
        grid_xz = griddata(grid[[0,2],:].T, uncertainty_percentage, (grid_x, grid_z), method='cubic')

        plt.imshow(grid_xz)


        plt.xlabel("X")
        plt.ylabel("Z")
        plt.title(r"Acceleration error $\left(\frac{\sqrt{\mathrm{trace}\left(P_{\mathbf{a}}\right)}}{\Vert \mathbf{a} \Vert }\cdot 100\ \%\right)$")
        plt.colorbar()
        plt.show()

       

    elif data["PROJECTION_AXIS"] == 2:
        x_min = np.amin(grid[0,:])
        x_max = np.amax(grid[0,:])

        y_min = np.amin(grid[2,:])
        y_max = np.amax(grid[2,:])

        grid_x, grid_y = np.mgrid[x_min:x_max:1000j, y_min:y_max:1000j]
        grid_xy = griddata(grid[[0,2],:].T, uncertainty_percentage, (grid_x, grid_y), method='cubic')

        plt.imshow(grid_xy)


        plt.xlabel("X")
        plt.ylabel("Y")
        plt.title(r"Acceleration error $\left(\frac{\sqrt{\mathrm{trace}\left(P_{\mathbf{a}}\right)}}{\Vert \mathbf{a} \Vert }\cdot 100\ \%\right)$")
        plt.colorbar()
        plt.show()

    else:
        raise(TypeError("Unrecognized case"))





plot_uncertainty_grid("/Users/bbercovici/GDrive/CUBoulder/Research/code/SBGAT/Examples/PGMUncertainty/output/PGMUncertainty_0")
# plot_uncertainty_grid("/Users/bbercovici/GDrive/CUBoulder/Research/code/SBGAT/Examples/PGMUncertainty/output/PGMUncertainty_1")
# plot_uncertainty_grid("/Users/bbercovici/GDrive/CUBoulder/Research/code/SBGAT/Examples/PGMUncertainty/output/PGMUncertainty_2")
