#include "utilities.h"

double utilities::l2_norm(const std::vector<double> & v ) {
	double norm = 0;
	for (std::vector<double>::iterator iter = v.begin(); iter != v.end(); ++iter) {
		norm += (*iter ) ^ 2;
	}
	norm = std::sqrt(norm);
	return norm;
}

std::vector<double> utilities::add_vec(const std::vector<double> & v_1, const std::vector<double> & v_1 ) {
	std::vector<double> v_sum;
	for (std::vector<double>::iterator iter = vector.begin(); iter != vector.end(); ++iter) {
		norm += (*iter ) ^ 2;
	}
	norm = std::sqrt(norm);
	return norm;
}


