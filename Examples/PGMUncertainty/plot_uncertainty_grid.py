import numpy as np
import matplotlib.pyplot as plt
import json
import scipy
from scipy.interpolate import griddata
from matplotlib.colors import LogNorm

def plot_uncertainty_grid(path):

    grid = np.loadtxt(path + "/grid.txt")
    uncertainty_percentage = np.loadtxt(path + "/uncertainty_over_reference_acc_percentage.txt")
    inside_outside = np.loadtxt(path + "/inside_outside.txt")
    
    data = {}
    with open(path + "/input_file.json") as f:
        data = json.load(f)

    for i in range(uncertainty_percentage.shape[0]):
            for j in range(uncertainty_percentage.shape[1]):
                if inside_outside[i,j] != 0:
                    uncertainty_percentage[i,j] = np.nan
    uncertainty_percentage = uncertainty_percentage.T
    

    min_uncertainty = np.min(uncertainty_percentage[np.logical_not(np.isnan(uncertainty_percentage))])
    max_uncertainty = np.max(uncertainty_percentage[np.logical_not(np.isnan(uncertainty_percentage))])


    levels = np.array([min_uncertainty,12.5,13,14,15])

    if data["PROJECTION_AXIS"] == 0:


        plt.imshow(uncertainty_percentage,interpolation = "bicubic",origin = "lower")
        plt.colorbar()

        cont = plt.contour(uncertainty_percentage,origin = "lower", levels = levels,
        	cmap = plt.get_cmap("Set1"))
        plt.gca().clabel(cont,cont.levels)
        plt.xlabel("Y")
        plt.ylabel("Z")
        plt.title(r"Relative Acceleration Uncertainty $\left(\frac{\sqrt{\mathrm{trace}\left(P_{\mathbf{a}}\right)}}{\Vert \mathbf{a} \Vert }\cdot 100\ \%\right)$")
        plt.tight_layout()
        plt.show()

     

    elif data["PROJECTION_AXIS"] == 1:


        
        plt.imshow(uncertainty_percentage,interpolation = "bicubic",origin = "lower")
        plt.colorbar()

        cont = plt.contour(uncertainty_percentage,origin = "lower", levels = levels,cmap = plt.get_cmap("rainbow"))
        plt.gca().clabel(cont,cont.levels)
        plt.xlabel("X")
        plt.ylabel("Z")
        plt.title(r"Relative Acceleration Uncertainty $\left(\frac{\sqrt{\mathrm{trace}\left(P_{\mathbf{a}}\right)}}{\Vert \mathbf{a} \Vert }\cdot 100\ \%\right)$")
        plt.tight_layout()
        plt.show()

       

    elif data["PROJECTION_AXIS"] == 2:
        
        plt.imshow(uncertainty_percentage,interpolation = "bicubic",origin = "lower")
        plt.colorbar()

        cont = plt.contour(uncertainty_percentage,origin = "lower", levels = levels,cmap = plt.get_cmap("rainbow"))
        plt.gca().clabel(cont,cont.levels)
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.title(r"Relative Acceleration Uncertainty $\left(\frac{\sqrt{\mathrm{trace}\left(P_{\mathbf{a}}\right)}}{\Vert \mathbf{a} \Vert }\cdot 100\ \%\right)$")
        plt.tight_layout()
        plt.show()

    else:
        raise(TypeError("Unrecognized case"))





plot_uncertainty_grid("/Users/bbercovici/GDrive/CUBoulder/Research/code/SBGAT/Examples/PGMUncertainty/output/PGMUncertainty_0")
plot_uncertainty_grid("/Users/bbercovici/GDrive/CUBoulder/Research/code/SBGAT/Examples/PGMUncertainty/output/PGMUncertainty_1")
plot_uncertainty_grid("/Users/bbercovici/GDrive/CUBoulder/Research/code/SBGAT/Examples/PGMUncertainty/output/PGMUncertainty_2")
