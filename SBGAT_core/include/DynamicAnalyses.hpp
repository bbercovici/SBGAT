#ifndef HEADER_DYNAMICANALYSES
#define HEADER_DYNAMICANALYSES

#include "ShapeModel.hpp"
#include <armadillo>
#include <boost/progress.hpp>

class DynamicAnalyses {

public:
	DynamicAnalyses(ShapeModel * shape_model);
	void compute_pgm(bool return_pgm = false);
	arma::vec pgm_acceleration(arma::vec & facet_center ) const ;

protected:
	ShapeModel * shape_model;

};


#endif