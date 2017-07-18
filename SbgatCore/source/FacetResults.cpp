#include "FacetResults.hpp"


using namespace SBGAT_CORE;

FacetResults::FacetResults() {

}

double FacetResults::get_grav_potential() const {
	return this -> grav_potential;
}

double FacetResults::get_grav_slope() const {
	return this -> grav_slope;


}

void FacetResults::set_grav_potential(double grav_potential) {
	this -> grav_potential = grav_potential;

}

void FacetResults::set_grav_slope(double grav_slope) {
	this -> grav_slope = grav_slope;
}


arma::vec * FacetResults::get_grav_acceleration() {
	return &this -> grav_acceleration;

}

void FacetResults::set_grav_acceleration(arma::vec & gravity_acceleration) {
	this -> grav_acceleration = gravity_acceleration;

}



