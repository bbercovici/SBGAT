import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import glob

def draw_slices_in_dir(input_dir,output_dir = "",prefix = ""):

	slice_counter = len(glob.glob1(input_dir,"slice_x_*"))
	print("Printing slices from " + str(slice_counter) + " MC shapes\n")

	if(slice_counter>0):

		slices_x = [input_dir + "/slice_x_"+str(i)+ ".txt" for i in range(slice_counter)]
		slices_y = [input_dir + "/slice_y_"+str(i)+ ".txt" for i in range(slice_counter)]
		slices_z = [input_dir + "/slice_z_"+str(i)+ ".txt" for i in range(slice_counter)]

		slices_x = [input_dir + "/baseline_slice_x.txt"] + slices_x
		slices_y = [input_dir + "/baseline_slice_y.txt"] + slices_y
		slices_z = [input_dir + "/baseline_slice_z.txt"] + slices_z

		draw_slice(0,slices_x,output_dir = output_dir,prefix = prefix)
		draw_slice(1,slices_y,output_dir = output_dir,prefix = prefix)
		draw_slice(2,slices_z,output_dir = output_dir,prefix = prefix)


def draw_slice(axis,slices,delay_plot = False,output_dir = "",prefix = ""):

	cut_names = ["Y-Z","X-Z","X-Y"]

	cmap = plt.cm.get_cmap(plt.cm.viridis)

	indices = range(len(slices))[::-1]

	for s in indices:

		path = slices[s]
		lines_file = np.loadtxt(path)
		
		x_max = - float("inf")
		x_min = float("inf")

		y_max = - float("inf")
		y_min = float("inf")

		if s == 0:
			c = "black"
			alpha = 1
		else:
			c = "lightblue"
			alpha = 0.7

		for i in range(lines_file.shape[0]):
			x = [ lines_file[i][0], lines_file[i][2]]
			y = [ lines_file[i][1], lines_file[i][3]]
			x_max = max(max(x),x_max)
			y_max = max(max(y),y_max)
			x_min = min(min(x),x_min)
			y_min = min(min(y),y_min)

			plt.gca().add_line(mpl.lines.Line2D(x, y,color = c,alpha = alpha))


	if axis == 0:
		plt.xlabel("Y (m)")
		plt.ylabel("Z (m)")
	elif axis == 1:
		plt.xlabel("X (m)")
		plt.ylabel("Z (m)")
	elif axis == 2:
		plt.xlabel("X (m)")
		plt.ylabel("Y (m)")

	
	plt.xlim(1.5 * x_min, 1.5 * x_max)
	plt.ylim(1.5 * y_min, 1.5 * y_max)

	plt.axis("equal")
	if (delay_plot):
		return
	if (len(output_dir) > 0):
			plt.savefig(output_dir + "/" + prefix +"slice_" + str(axis) + ".pdf", bbox_inches='tight')
	else:
		plt.show()

	plt.cla()
	plt.clf()

# draw_slices_in_dir("/Users/bbercovici/GDrive/CUBoulder/Research/code/SBGAT/Examples/PGMUncertaintyMCSlopes/output/PGMUncertaintyMCSlopes_0")
