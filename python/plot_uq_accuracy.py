import numpy as np
import matplotlib.pyplot as plt
import json
from draw_shape_slices import draw_slice
import glob

def plot_uq_accuracy(path):

    all_positions = np.loadtxt(path + "/all_positions_arma.txt")

    abs_value_cov_difference_analytical_vs_mc = np.loadtxt(path + "/abs_value_cov_difference_analytical_vs_mc.txt")
    rel_value_cov_difference_analytical_vs_mc = np.loadtxt(path + "/rel_value_cov_difference_analytical_vs_mc.txt")
    KL_divergence_analytical_vs_mc = np.loadtxt(path + "/KL_divergence_analytical_vs_mc.txt")

    data = {}
    with open(path + "/input_file.json") as f:
        data = json.load(f)

    if data["PROJECTION_AXIS"] == 0:

        slice_counter = len(glob.glob1(path,"slice_x_*"))
        draw_slice(0, [path + "/baseline_slice_x.txt"] ,delay_plot = True)
    
        plotted_points = [ np.abs(all_positions[0,i]) == 0 for i in range(all_positions.shape[1]) ]
        plt.scatter(all_positions[1,plotted_points],all_positions[2,plotted_points],c = np.log10(abs_value_cov_difference_analytical_vs_mc[plotted_points]))
        plt.colorbar()
        plt.ylabel("Z(m)")
        plt.xlabel("Y(m)")
        plt.title(r"$\mathrm{log}\left(\max \vert P_{\mathbf{a},mc} - P_{\mathbf{a}} \vert\right)$")
        plt.show()

        draw_slice(0, [path + "/baseline_slice_x.txt"] ,delay_plot = True)
        

        plt.scatter(all_positions[1,plotted_points],all_positions[2,plotted_points],c = rel_value_cov_difference_analytical_vs_mc[plotted_points])
        plt.colorbar()
        plt.ylabel("Z(m)")
        plt.xlabel("Y(m)")
        plt.title(r"$\frac{\sqrt{\Vert P_{\mathbf{a},mc} - P_{\mathbf{a}}\Vert}}{\mathrm{E}\left(\mathbf{a}_{mc}\right)}$")

        plt.show()
        draw_slice(0, [path + "/baseline_slice_x.txt"] ,delay_plot = True)
        
        plt.scatter(all_positions[1,plotted_points],all_positions[2,plotted_points],c = KL_divergence_analytical_vs_mc[plotted_points])
        plt.colorbar()
        plt.ylabel("Z(m)")
        plt.xlabel("Y(m)")
        plt.title(r"$\mathrm{KL}\left(\mathbf{a}_{ref},P_{\mathbf{a}}\Vert\mathbf{a}_{mc},P_{\mathbf{a}_{mc}}\right)$")

        plt.show()

     

    elif data["PROJECTION_AXIS"] == 1:

        slice_counter = len(glob.glob1(path,"slice_x_*"))
        draw_slice(1, [path + "/baseline_slice_y.txt"] ,delay_plot = True)
    
        plotted_points = [ np.abs(all_positions[1,i]) == 0 for i in range(all_positions.shape[1]) ]
        plt.scatter(all_positions[0,plotted_points],all_positions[2,plotted_points],c = np.log10(abs_value_cov_difference_analytical_vs_mc[plotted_points]))
        plt.colorbar()
        plt.ylabel("X(m)")
        plt.xlabel("Z(m)")
        plt.title(r"$\mathrm{log}\left(\max \vert P_{\mathbf{a},mc} - P_{\mathbf{a}} \vert\right)$")
        plt.show()

        draw_slice(1, [path + "/baseline_slice_y.txt"] ,delay_plot = True)
        
        plt.scatter(all_positions[0,plotted_points],all_positions[2,plotted_points],c = rel_value_cov_difference_analytical_vs_mc[plotted_points])
        plt.colorbar()
        plt.ylabel("X(m)")
        plt.xlabel("Z(m)")
        plt.title(r"$\frac{\sqrt{\Vert P_{\mathbf{a},mc} - P_{\mathbf{a}}\Vert}}{\mathrm{E}\left(\mathbf{a}_{mc}\right)}$")

        plt.show()
        draw_slice(1, [path + "/baseline_slice_y.txt"] ,delay_plot = True)
        
        plt.scatter(all_positions[0,plotted_points],all_positions[2,plotted_points],c = KL_divergence_analytical_vs_mc[plotted_points])
        plt.colorbar()
        plt.ylabel("X(m)")
        plt.xlabel("Z(m)")
        plt.title(r"$\mathrm{KL}\left(\mathbf{a}_{ref},P_{\mathbf{a}}\Vert\mathbf{a}_{mc},P_{\mathbf{a}_{mc}}\right)$")

        plt.show()

        

     
    elif data["PROJECTION_AXIS"] == 2:
        slice_counter = len(glob.glob1(path,"slice_x_*"))
        draw_slice(2, [path + "/baseline_slice_z.txt"] ,delay_plot = True)
    
        plotted_points = [ np.abs(all_positions[2,i]) == 0 for i in range(all_positions.shape[1]) ]
        plt.scatter(all_positions[0,plotted_points],all_positions[1,plotted_points],c = np.log10(abs_value_cov_difference_analytical_vs_mc[plotted_points]))
        plt.colorbar()
        plt.ylabel("Y(m)")
        plt.xlabel("Z(m)")
        plt.title(r"$\mathrm{log}\left(\max \vert P_{\mathbf{a},mc} - P_{\mathbf{a}} \vert\right)$")
        plt.show()

        draw_slice(2, [path + "/baseline_slice_z.txt"] ,delay_plot = True)
        

        plt.scatter(all_positions[0,plotted_points],all_positions[1,plotted_points],c = rel_value_cov_difference_analytical_vs_mc[plotted_points])
        plt.colorbar()
        plt.ylabel("Y(m)")
        plt.xlabel("Z(m)")
        plt.title(r"$\frac{\sqrt{\Vert P_{\mathbf{a},mc} - P_{\mathbf{a}}\Vert}}{\mathrm{E}\left(\mathbf{a}_{mc}\right)}$")

        plt.show()
        draw_slice(2, [path + "/baseline_slice_z.txt"] ,delay_plot = True)
        
        plt.scatter(all_positions[0,plotted_points],all_positions[1,plotted_points],c = KL_divergence_analytical_vs_mc[plotted_points])
        plt.colorbar()
        plt.ylabel("Y(m)")
        plt.xlabel("Z(m)")
        plt.title(r"$\mathrm{KL}\left(\mathbf{a}_{ref},P_{\mathbf{a}}\Vert\mathbf{a}_{mc},P_{\mathbf{a}_{mc}}\right)$")

        plt.show()

        

    else:
        raise(TypeError("Unrecognized case"))







plot_uq_accuracy("/Users/bbercovici/GDrive/CUBoulder/Research/code/SBGAT/Examples/PGMUncertainty/output/PGMUncertainty_0")
plot_uq_accuracy("/Users/bbercovici/GDrive/CUBoulder/Research/code/SBGAT/Examples/PGMUncertainty/output/PGMUncertainty_1")
plot_uq_accuracy("/Users/bbercovici/GDrive/CUBoulder/Research/code/SBGAT/Examples/PGMUncertainty/output/PGMUncertainty_2")
