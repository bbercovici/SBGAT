#include "SyntheticObservations.hpp"


using namespace SBGAT_CORE;



SyntheticObservations::SyntheticObservations(ShapeModel * shape_model, FrameGraph * frame_graph) {

	this -> shape_model = shape_model;
	this -> frame_graph = frame_graph;
}

void SyntheticObservations::compute_lightcurve_fixed_spin(
    arma::vec & spin_axis,
    double spin_rate,
    double radius,
    double t0,
    double tf,
    double initial_phase,
    double dt) {


	// Time vector
	arma::vec times ;
	if (dt > 0) {
		times = arma::linspace(t0, tf, (int)((tf - t0) / dt));
	}
	else {
		times = arma::linspace(t0, tf);
	}

	// Resetting data container
	this -> lightcurve_data.reset();
	this -> lightcurve_data = arma::mat(times.n_rows, 3);

	// Phase angle vector
	arma::vec phase_angle(times.n_rows);

	// Brightness vector
	arma::vec brightness = arma::vec(times.n_rows);

	// Unit vector
	arma::vec u = {1, 0, 0};

	std::cout << " Computing lightcurve " << std::endl;

	boost::progress_display progress(times.n_rows) ;

	for (unsigned int time_index = 0; time_index < times.n_rows; ++time_index) {

		// Time (years)
		double t = times(time_index);

		// Earth angular position (rad)
		double theta_earth = 2 * arma::datum::pi * (t - times(0)) ;

		// Body angular position (rad)
		double theta_body =  2 * arma::datum::pi * (t - times(0)) / (std::sqrt(radius * radius * radius));

		// Set frames origin
		arma::vec sun_to_earth_N = RBK::M3(theta_earth ).t() * u ;
		arma::vec sun_to_body_N = radius * RBK::M3(theta_body + initial_phase).t() * u ;

		this -> frame_graph -> set_transform_origin("N", "E", sun_to_earth_N);
		this -> frame_graph -> set_transform_origin("N", "T", sun_to_body_N) ;

		// Set body attitude
		double spin_angle = (t - times(0)) * spin_rate;
		arma::vec prv = spin_angle * spin_axis;
		arma::vec mrp = RBK::prv_to_mrp(prv);

		this -> frame_graph -> set_transform_mrp("N", "T", mrp);

		// Directions
		arma::vec earth_to_body_dir_N = arma::normalise(sun_to_body_N - sun_to_earth_N);
		arma::vec sun_to_body_dir_N = arma::normalise(sun_to_body_N);

		arma::vec earth_to_body_dir_T = this -> frame_graph -> convert(earth_to_body_dir_N, "N",
		                                "T", true);
		arma::vec sun_to_body_dir_T = this -> frame_graph -> convert(sun_to_body_dir_N, "N",
		                              "T", true);

		// Phase angle
		phase_angle(time_index) = std::asin(
		                              arma::norm(
		                                  arma::cross(sun_to_body_dir_N, arma::normalise(sun_to_earth_N)))
		                              * arma::norm(sun_to_earth_N) / arma::norm(sun_to_body_N - sun_to_earth_N));

		// Brightness
		brightness(time_index) = this -> collect_brightness(earth_to_body_dir_T, sun_to_body_dir_T);

		++progress;
	}

	this -> lightcurve_data.col(0) = times;
	this -> lightcurve_data.col(1) = phase_angle;
	this -> lightcurve_data.col(2) = brightness;

}


double SyntheticObservations::collect_brightness(arma::vec & earth_to_body_dir_T, arma::vec & sun_to_body_dir_T) {

	double brightness = 0;

	// The brightness reflected towards Earth is accumulated over each facet
	#pragma omp parallel for reduction(+:brightness)
	for (unsigned int facet_index = 0; facet_index < this -> shape_model -> get_NFacets(); ++facet_index) {

		Facet * facet = this -> shape_model -> get_facets() -> at(facet_index);

		arma::vec * n = facet -> get_facet_normal();

		double mu_0 = arma::dot(*n, - sun_to_body_dir_T);
		double mu = arma::dot(*n, - earth_to_body_dir_T);

		if (mu_0 > 0 && mu > 0) {
			brightness += facet -> get_albedo() * mu_0 * mu * ( 1. / (mu_0 + mu) + 0.1) * facet -> get_area();
		}

	}

	return brightness;
}

void SyntheticObservations::save_lightcurve(std::string filepath) const {
	this -> lightcurve_data.save(filepath, arma::raw_ascii);
}
