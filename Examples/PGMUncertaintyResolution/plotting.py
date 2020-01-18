import matplotlib.pyplot as plt
import numpy as np
import json



def load_haussdorff_distance(path):

	distance = np.loadtxt(path)
	distance = distance[:,-1]
	return distance


	distance = load_haussdorff_distance("/Users/bbercovici/Desktop/colored_mesh.ply")

	distance_sim = np.array(list(-distance) + list(distance))
	sd = np.std(distance_sim)
	print(sd)
	print(np.std(distance_sim))
	a = plt.hist(distance_sim,bins = 50)
	plt.plot(np.linspace(-20,20,1000),max(a[0]) * np.exp(-1./(sd ** 2 * 2) * np.linspace(-20,20,1000)* np.linspace(-20,20,1000)))

	plt.plot(np.linspace(-20,20,1000),max(a[0]) * np.exp(-1./(3 ** 2 * 2) * np.linspace(-20,20,1000)* np.linspace(-20,20,1000)))

	plt.show()



def plot_results(save = False):

	base_path = "output/PGMUncertaintyResolution_"
	for i in range(5):
		full_path = base_path + str(i)

		data = {}
		with open(full_path + "/input_file.json") as f:
			data = json.load(f)

		distance = np.loadtxt(full_path + "/all_positions_arma_normalized_acceleration_difference_sd.txt")
		plt.scatter(range(len(distance)),distance,label = "$l =" + str(data["CORRELATION_DISTANCE"]) + "\ \mathrm{m}$" )

	plt.ylabel(r"$\sqrt{\delta\mathbf{a}^T\left[P_{\mathbf{a}}\right]^{-1}\delta\mathbf{a}}$")
	plt.xlabel("Queried point index")
	plt.legend(loc = "upper center",ncol = 5,bbox_to_anchor = (0.5,1.1))
	plt.tight_layout()
	if save:
		# plt.savefig("/Users/bbercovici/GDrive/CUBoulder/Research/thesis/figs/pgm_uncertainty/pgm_uncertainty_resolution_correlation_distance.pdf")
		plt.savefig("/Users/bbercovici/GDrive/CUBoulder/Research/papers/pgm_uncertainty/R1/all_figures/pgm_uncertainty_resolution_correlation_distance.pdf")
		
	plt.show()

plot_results(True)
