#include "../include/Vertex.hpp"

void Vertex::set_coordinates(std::shared_ptr<arma::vec> coordinates) {
	this -> coordinates = coordinates;
}


void Vertex::add_facet_ownership(Facet * facet) {
	
	bool already_present = false;

	for (unsigned int i = 0; i < this -> owning_facets.size(); ++i) {

		// If this vertex is alredy owned by $facet, nothing to do here

		if (this -> owning_facets[i] == facet) {
			already_present = true;
			break;
		}
	}
	if (already_present == false) {
		this -> owning_facets.push_back(facet);
	}


}


std::vector<Facet *>  Vertex::common_facets(std::shared_ptr<Vertex> vertex) const {

	std::vector<Facet *> common_facets;

	for (std::vector<Facet * >::const_iterator it = this -> owning_facets.begin();
	        it != this -> owning_facets.end(); ++it) {

		if (vertex -> is_owned_by(*it)) {
			common_facets.push_back(*it);
		}

	}

	return common_facets;

}




bool Vertex::is_owned_by(Facet * facet) const {
	for (std::vector<Facet * >::const_iterator it = this -> owning_facets.begin();
	        it != this -> owning_facets.end(); ++it) {

		if (*it == facet) {
			return true;
		}

	}
	return false;
}


arma::vec * Vertex::get_coordinates()  {
	return this -> coordinates.get();
}

unsigned int Vertex::get_number_of_owning_facets() const {
	return this -> owning_facets.size();
}
