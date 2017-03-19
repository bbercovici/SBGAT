#include "DynamicAnalyses.hpp"

DynamicAnalyses::DynamicAnalyses(ShapeModel * shape_model) {
	this -> shape_model = shape_model;
}

void DynamicAnalyses::compute_pgm(double density, bool return_pgm ) {

	arma::mat gravity_accelerations = arma::mat(3, this -> shape_model -> get_NFacets());

	std::cout << std::endl << " Computing PGM " << std::endl ;
	boost::progress_display progress(this -> shape_model -> get_NFacets()) ;

	for (unsigned int facet_index = 0; facet_index < this -> shape_model -> get_NFacets(); ++facet_index) {

		gravity_accelerations.col(facet_index) = this -> pgm_acceleration(
		            this -> shape_model -> get_facets() -> at(facet_index) -> get_facet_center() -> colptr(0) ,
		            density);

		++progress;
	}

}

arma::vec DynamicAnalyses::pgm_acceleration(double * point , double density) const {

	double ax = 0;
	double ay = 0;
	double az = 0;

	// Facet loop
	#pragma omp parallel for reduction(+:ax,ay,az) if (USE_OMP_DYNAMIC_ANALYSIS)
	for (unsigned int facet_index = 0; facet_index < this -> shape_model -> get_NFacets(); ++ facet_index) {

		std::vector<std::shared_ptr<Vertex > > * vertices = this -> shape_model -> get_facets() -> at(facet_index) -> get_vertices();

		const double * r1 =  vertices -> at(0) -> get_coordinates() -> colptr(0);
		const double * r2 =  vertices -> at(1) -> get_coordinates() -> colptr(0);
		const double * r3 =  vertices -> at(2) -> get_coordinates() -> colptr(0);

		double r1m[3];
		double r2m[3];
		double r3m[3];

		r1m[0] = r1[0] - point[0];
		r1m[1] = r1[1] - point[1];
		r1m[2] = r1[2] - point[2];

		r2m[0] = r2[0] - point[0];
		r2m[1] = r2[1] - point[1];
		r2m[2] = r2[2] - point[2];

		r3m[0] = r3[0] - point[0];
		r3m[1] = r3[1] - point[1];
		r3m[2] = r3[2] - point[2];



		double R1 = std::sqrt( r1m[0] * r1m[0]
		                       + r1m[1] * r1m[1]
		                       + r1m[2] * r1m[2]       );

		double R2 = std::sqrt( r2m[0] * r2m[0]
		                       + r2m[1] * r2m[1]
		                       + r2m[2] * r2m[2]      );


		double R3 = std::sqrt( r3m[0] * r3m[0]
		                       + r3m[1] * r3m[1]
		                       + r3m[2] * r3m[2]      );

		double r2_cross_r3_0 = r2m[1] * r3m[2] - r2m[2] * r3m[1];
		double r2_cross_r3_1 = r3m[0] * r2m[2] - r3m[2] * r2m[0];
		double r2_cross_r3_2 = r2m[0] * r3m[1] - r2m[1] * r3m[0];



		double wf = 2 * std::atan2(
		                r1m[0] * r2_cross_r3_0 + r1m[1] * r2_cross_r3_1 + r1m[2] * r2_cross_r3_2,

		                R1 * R2 * R3 + R1 * (r2m[0] * r3m[0] + r2m[1] * r3m[1]  + r2m[2] * r3m[2] )
		                + R2 * (r3m[0] * r1m[0] + r3m[1] * r1m[1] + r3m[2] * r1m[2])
		                + R3 * (r1m[0] * r2m[0] + r1m[1] * r2m[1] + r1m[2] * r2m[2]));


		arma::mat * Fdyad = this -> shape_model -> get_facets() -> at(facet_index) -> get_facet_dyad();

		double * F_col_0 = Fdyad -> colptr(0);
		double * F_col_1 = Fdyad -> colptr(1);
		double * F_col_2 = Fdyad -> colptr(2);

		ax += wf * (F_col_0[0] * r1m[0] + F_col_1[0] * r1m[1] +  F_col_2[0] * r1m[2]);
		ay += wf * (F_col_0[1] * r1m[0] + F_col_1[1] * r1m[1] +  F_col_2[1] * r1m[2]);
		az += wf * (F_col_0[2] * r1m[0] + F_col_1[2] * r1m[1] +  F_col_2[2] * r1m[2]);

	}


	// Edge loop
	#pragma omp parallel for reduction(-:ax,ay,az) if (USE_OMP_DYNAMIC_ANALYSIS)
	for (unsigned int edge_index = 0; edge_index < this -> shape_model -> get_NEdges(); ++ edge_index) {


		const double * r1 =  this -> shape_model -> get_edges() -> at(edge_index) -> get_v0() -> get_coordinates() -> colptr(0);
		const double * r2 =  this -> shape_model -> get_edges() -> at(edge_index) -> get_v1() -> get_coordinates() -> colptr(0);


		double r1m[3];
		double r2m[3];

		r1m[0] = r1[0] - point[0];
		r1m[1] = r1[1] - point[1];
		r1m[2] = r1[2] - point[2];

		r2m[0] = r2[0] - point[0];
		r2m[1] = r2[1] - point[1];
		r2m[2] = r2[2] - point[2];


		double R1 = std::sqrt( r1m[0] * r1m[0]
		                       + r1m[1] * r1m[1]
		                       + r1m[2] * r1m[2]       );

		double R2 = std::sqrt( r2m[0] * r2m[0]
		                       + r2m[1] * r2m[1]
		                       + r2m[2] * r2m[2]      );

		double Re = std::sqrt( (r2m[0] - r1m[0]) * (r2m[0] - r1m[0])
		                       + (r2m[1] - r1m[1]) * (r2m[1] - r1m[1])
		                       + (r2m[2] - r1m[2]) * (r2m[2] - r1m[2])      );


		double Le = std::log((R1 + R2 + Re) / (R1 + R2 - Re));


		arma::mat * Edyad = this -> shape_model -> get_edges() -> at(edge_index) -> get_edge_dyad();

		double * E_col_0 = Edyad -> colptr(0);
		double * E_col_1 = Edyad -> colptr(1);
		double * E_col_2 = Edyad -> colptr(2);


		ax -= Le * (E_col_0[0] * r1m[0] + E_col_1[0] * r1m[1] +  E_col_2[0] * r1m[2]);
		ay -= Le * (E_col_0[1] * r1m[0] + E_col_1[1] * r1m[1] +  E_col_2[1] * r1m[2]);
		az -= Le * (E_col_0[2] * r1m[0] + E_col_1[2] * r1m[1] +  E_col_2[2] * r1m[2]);


	}

	arma::vec acceleration = {ax, ay, az};
	acceleration = acceleration * arma::datum::G * density;

	return acceleration;

}
