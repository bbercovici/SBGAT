import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Ellipse
import np_array_to_latex as nptlatex
from np_array_to_latex import np_array_to_latex
from pprint import pprint
import json
import os

def draw_slice(axis,slices,output_dir = None,zoom = False,prefix = ""):

	cut_names = ["Y-Z","X-Z","X-Y"]
	zoom_options = " with zoom ... \n" if zoom else " without zoom ... \n"

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
			x = [1e3 * lines_file[i][0],1e3 * lines_file[i][2]]
			y = [1e3 * lines_file[i][1],1e3 * lines_file[i][3]]
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

	plt.scatter(0,0,marker = ".",color = "black" )
	
	if zoom is False:
		plt.xlim(1.5 * x_min, 1.5 * x_max)
		plt.ylim(1.5 * y_min, 1.5 * y_max)

	plt.axis("equal")

	if (len(output_dir) > 0):
		if (zoom is True):
			plt.savefig(output_dir + "/" + prefix +"slice_zoom_" + str(axis) + ".pdf", bbox_inches='tight')
		else:
			plt.savefig(output_dir + "/" + prefix +"slice_" + str(axis) + ".pdf", bbox_inches='tight')
	else:
		plt.show()

	plt.cla()
	plt.clf()
