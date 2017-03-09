#include <string>
#include <iostream>
#include <armadillo>


class ShapeModel : SBGAT_class{

public:

	int load_shape_model(std::string & filename);
	int save_shape_model(std::string & filename);


protected:

	void compute_normals();

	arma::mat vertex_coordinates;
	arma::mat cell_vertices;
	

};