#ifndef HEADER_SYNTHETICOBSERVATIONS
#define HEADER_SYNTHETICOBSERVATIONS

#include "ShapeModel.hpp"
#include <armadillo>
#include <boost/progress.hpp>
#include "omp.h"
#include "OMP_flags.hpp"
#include "RigidBodyKinematics.hpp"

class SyntheticObservations {

public:
	/**
	Constructor
	Creates an instance of a SyntheticObservations object
	@param shape_model Pointer to shape model from which synthetic observations are generated
	@param frame_graph Pointer to the graph of reference frame involving the provided body and other
	entities
	*/
	SyntheticObservations(ShapeModel * shape_model, FrameGraph * frame_graph);

	/**
	Computes the lightcurve of the present body assuming a constant spin axis (in body-fixed frame), spin rate
	. Assumes Circular orbit of the earth with respect to the sun

	@param spin_axis Spin axis direction in body-fixed frame
	@param radius Sun/body distance (AU)
	@param t0 Initial time (years)
	@param tf Final time (years)
	@param initial_phase Phase offset between the earth and the body. If positive and less than 2pi, the body is initally leading the Earth
	@param dt Time step (years). Ignored if negative.
	*/
	void compute_lightcurve_fixed_spin(
	    arma::vec & spin_axis,
	    double spin_rate,
	    double radius,
	    double t0,
	    double tf,
	    double initial_phase,
	    double dt = -1);


	/**
	Save the lightcurve data to a prescribed file
	@param path Filepath
	*/
	void save_lightcurve(std::string filepath) const;


protected:
	ShapeModel * shape_model;
	FrameGraph * frame_graph;
	double collect_brightness(arma::vec & earth_to_body_dir_T, arma::vec & sun_to_body_dir_T);

	// Lightcurve data.
	// - First column : time (years)
	// - Second column : phase angle (rad)
	// - Third column : brightness
	arma::mat lightcurve_data;




};


#endif