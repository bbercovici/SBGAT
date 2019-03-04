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
        

        y_min = np.amin(grid[1,:])
        y_max = np.amax(grid[1,:])
        
        z_min = np.amin(grid[2,:])
        z_max = np.amax(grid[2,:])

       
        grid_y, grid_z = np.mgrid[y_min:y_max:1000j, z_min:z_max:1000j]
        grid_yz = griddata(grid[[1,2],:].T, uncertainty_percentage, (grid_y, grid_z), method='cubic')

        plt.imshow(grid_yz)
        plt.show()

        plt.scatter(grid[1,:],grid[2,:],c = uncertainty_percentage)
        plt.colorbar()
        plt.show()


    elif data["PROJECTION_AXIS"] == 1:


        
        x_min = np.amin(grid[0,:])
        x_max = np.amax(grid[0,:])

        z_min = np.amin(grid[2,:])
        z_max = np.amax(grid[2,:])

        grid_x, grid_z = np.mgrid[x_min:x_max:1000j, z_min:z_max:1000j]
        grid_xz = griddata(grid[[0,2],:].T, uncertainty_percentage, (grid_x, grid_z), method='cubic')

        plt.imshow(grid_xz)
        plt.show()




        plt.scatter(grid[2,:],grid[0,:],c = uncertainty_percentage)
        plt.colorbar()
        plt.show()
    elif data["PROJECTION_AXIS"] == 2:
        plt.scatter(grid[0,:],grid[1,:],c = uncertainty_percentage)
        plt.colorbar()
        plt.show()
    else:
        raise(TypeError("Unrecognized case"))





# plot_uncertainty_grid("/Users/bbercovici/GDrive/CUBoulder/Research/code/SBGAT/Examples/PGMUncertainty/output/PGMUncertainty_0")
plot_uncertainty_grid("/Users/bbercovici/GDrive/CUBoulder/Research/code/SBGAT/Examples/PGMUncertainty/output/PGMUncertainty_1")
# plot_uncertainty_grid("/Users/bbercovici/GDrive/CUBoulder/Research/code/SBGAT/Examples/PGMUncertainty/output/PGMUncertainty_2")
