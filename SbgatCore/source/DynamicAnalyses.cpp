#include "DynamicAnalyses.hpp"
#include "PolyhedralSphericalHarmo.hpp"


using namespace SBGAT_CORE;

DynamicAnalyses::DynamicAnalyses(ShapeModel * shape_model) {
	this -> shape_model = shape_model;
}

void DynamicAnalyses::compute_pgm_accelerations(double density ) {

	arma::vec acc(3);
	for (unsigned int facet_index = 0; facet_index < this -> shape_model -> get_NFacets(); ++facet_index) {
		Facet * facet = this -> shape_model -> get_facets() -> at(facet_index);

		acc = this -> pgm_acceleration(
			facet -> get_facet_center() -> colptr(0) ,
			density);
		facet -> get_facet_results() -> set_grav_acceleration(acc);

	}

}



void DynamicAnalyses::compute_pgm_potentials(double density) {

	for (unsigned int facet_index = 0; facet_index < this -> shape_model -> get_NFacets(); ++facet_index) {

		Facet * facet = this -> shape_model -> get_facets() -> at(facet_index);

		double potential = this -> pgm_potential(
			facet -> get_facet_center() -> colptr(0) ,
			density);
		facet -> get_facet_results() -> set_grav_potential(potential);

	}

}


arma::vec DynamicAnalyses::compute_gravity_slopes(
	arma::vec spin_axis,
	double spin_rate) {

	arma::vec omega_body_wrt_inertial = spin_rate * spin_axis;
	arma::vec stats;

	double min_slope = 181;
	double max_slope = -181;
	double mean_slope = 0;


	#pragma omp parallel for reduction(+:mean_slope), reduction(min:min_slope), reduction(max:max_slope)
	for (unsigned int facet_index = 0; facet_index < this -> shape_model -> get_NFacets();
		++facet_index) {

		Facet * facet = this -> shape_model -> get_facets() -> at(facet_index);

	arma::vec * n = facet -> get_facet_normal();
	double slope_d = std::acos(-arma::dot(
		*n,
		arma::normalise(*facet -> get_facet_results() -> get_grav_acceleration() - (
			arma::cross(
				omega_body_wrt_inertial,
				arma::cross(omega_body_wrt_inertial, *facet -> get_facet_center())))
		))) * 180 / arma::datum::pi ;

		// Statistics are accumulated
	if (slope_d > max_slope) {
		max_slope = slope_d;
	}

	if (slope_d < min_slope) {
		min_slope = slope_d;
	}
	mean_slope += slope_d / this -> shape_model -> get_NFacets();

		// The result for this facet is saved
	facet -> get_facet_results() -> set_grav_slope(slope_d);

}

stats = {min_slope, mean_slope, max_slope};
return stats;

}


double DynamicAnalyses::pgm_potential(double * point , double density) const {

	double potential = 0;

	// Facet loop
	#pragma omp parallel for reduction(+:potential) if (USE_OMP_DYNAMIC_ANALYSIS)
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

		double ax = wf * (F_col_0[0] * r1m[0] + F_col_1[0] * r1m[1] +  F_col_2[0] * r1m[2]);
		double ay = wf * (F_col_0[1] * r1m[0] + F_col_1[1] * r1m[1] +  F_col_2[1] * r1m[2]);
		double az = wf * (F_col_0[2] * r1m[0] + F_col_1[2] * r1m[1] +  F_col_2[2] * r1m[2]);

		potential += - (ax * r1m[0] + ay * r1m[1] + az * r1m[2]);

	}


	// Edge loop
	#pragma omp parallel for reduction(+:potential) if (USE_OMP_DYNAMIC_ANALYSIS)
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

		double ax = Le * (E_col_0[0] * r1m[0] + E_col_1[0] * r1m[1] +  E_col_2[0] * r1m[2]);
		double ay = Le * (E_col_0[1] * r1m[0] + E_col_1[1] * r1m[1] +  E_col_2[1] * r1m[2]);
		double az = Le * (E_col_0[2] * r1m[0] + E_col_1[2] * r1m[1] +  E_col_2[2] * r1m[2]);

		potential += (ax * r1m[0] + ay * r1m[1] + az * r1m[2]);

	}
	potential *= 0.5 * density * arma::datum::G;


	return potential;

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

void DynamicAnalyses::compute_exterior_sh_coefs(
	arma::mat & Cnm_total,
	arma::mat & Snm_total,
	unsigned int n_degree,
	double ref_radius,
	double density,
	bool normalized) {

	// Preallocation
	double total_mass = this -> shape_model -> get_volume() * density;

	// Normalized coefficients
	Cnm_total = arma::zeros<arma::mat>(n_degree + 1  , n_degree + 1);
	Snm_total = arma::zeros<arma::mat>(n_degree + 1  , n_degree + 1);

	// Main loop
	unsigned int num_facet = this -> shape_model -> get_NFacets();


	for (unsigned int facet_index = 0 ; facet_index <  num_facet; ++ facet_index) { // Keep adding differential C * Mass and S * Mass

		// Volume of tetrahedron associated with current facet
		std::vector<std::shared_ptr<SBGAT_CORE::Vertex > > * vertices = this -> shape_model -> get_facets() -> at(facet_index) -> get_vertices();

		arma::vec * r0 =  vertices -> at(0) -> get_coordinates();
		arma::vec * r1 =  vertices -> at(1) -> get_coordinates();
		arma::vec * r2 =  vertices -> at(2) -> get_coordinates();


		double dv = arma::dot(*r0, arma::cross(*r1 - *r0, *r2 - *r0)) / 6.;

		// Compute the CS for each polyhedron

		arma::mat Cnm2f ;
		arma::mat Snm2f ;

		ComputePolyhedralCS(
			Cnm2f,
			Snm2f,
			n_degree,
			ref_radius,
			dv * density,
			density,
			r0 -> colptr(0),
			r1 -> colptr(0),
			r2 -> colptr(0),
			total_mass,
			normalized
		);// % For normalized coeffs


		// Total \sum (C_i * Mass_i) for each density region

		Cnm_total += Cnm2f * dv * density;
		Snm_total += Snm2f * dv * density;

	}

	// Compute the overall coefficient values
	Cnm_total = Cnm_total / total_mass;
	Snm_total = Snm_total / total_mass;

}



void DynamicAnalyses::GetBnmNormalizedExterior(int n_degree,
	arma::mat & b_bar_real,
	arma::mat & b_bar_imag,
	const arma::vec & pos,
	double ref_radius){


	double r_sat = arma::norm(pos);
	double x_sat = pos(0);
	double y_sat = pos(1);
	double z_sat = pos(2);

	double delta_1_n;


    /*////////////////////////////////
    // -- Vertical Recurrences -- //
    ////////////////////////////////*/

	for (int mm = 0; mm <= n_degree+2; mm++) {


      // Pretty sure I don't need this

		double m = (double) mm;

		for (int nn = mm; nn <= n_degree+2; nn++) {

      // Pretty sure I don't need this

			double n = (double) nn;

            /*// Recursive Formulae*/

			if (mm == nn) {

				if (mm == 0) {

					b_bar_real(0,0) = ref_radius/r_sat;
					b_bar_imag(0,0) = 0.0;

				} else {

					if (nn == 1) {

						delta_1_n = 1.0;

					} else {

						delta_1_n = 0.0;

					} 

					b_bar_real(nn,nn) = sqrt( (1.0 + delta_1_n)*(2.0*n + 1.0)/(2.0*n) ) * (ref_radius/r_sat) * ( x_sat/r_sat*b_bar_real(nn-1,nn-1) - y_sat/r_sat*b_bar_imag(nn-1,nn-1) );
					b_bar_imag(nn,nn) = sqrt( (1.0 + delta_1_n)*(2.0*n + 1.0)/(2.0*n) ) * (ref_radius/r_sat) * ( y_sat/r_sat*b_bar_real(nn-1,nn-1) + x_sat/r_sat*b_bar_imag(nn-1,nn-1) );

				} 

			} 

			else {

				if ( nn >= 2 ) {

					b_bar_real(nn,mm) = sqrt( (4.0*n*n - 1.0)/(n*n - m*m) )*(ref_radius/r_sat)*(z_sat/r_sat)*b_bar_real(nn-1,mm) - sqrt( (2.0*n + 1.0)*( (n - 1.0)*(n - 1.0) - m*m )/( (2.0*n - 3.0)*(n*n - m*m) ) )*(ref_radius/r_sat)*(ref_radius/r_sat)*b_bar_real(nn-2,mm);
					b_bar_imag(nn,mm) = sqrt( (4.0*n*n - 1.0)/(n*n - m*m) )*(ref_radius/r_sat)*(z_sat/r_sat)*b_bar_imag(nn-1,mm) - sqrt( (2.0*n + 1.0)*( (n - 1.0)*(n - 1.0) - m*m )/( (2.0*n - 3.0)*(n*n - m*m) ) )*(ref_radius/r_sat)*(ref_radius/r_sat)*b_bar_imag(nn-2,mm);

				} else {

					b_bar_real(nn,mm) = sqrt( (4.0*n*n - 1.0)/(n*n - m*m) )*(ref_radius/r_sat)*(z_sat/r_sat)*b_bar_real(nn-1,mm);
					b_bar_imag(nn,mm) = sqrt( (4.0*n*n - 1.0)/(n*n - m*m) )*(ref_radius/r_sat)*(z_sat/r_sat)*b_bar_imag(nn-1,mm);

				} 

			} 

		} 

	} 

} 


arma::vec DynamicAnalyses::spherical_harmo_acc(const unsigned int n_degree,
	const double ref_radius,
	const double  mu,
	const arma::vec & pos, 
	const arma::mat & Cbar,
	const arma::mat & Sbar) {

	int n_max = 50;

	if (Cbar.n_rows < n_degree){
		throw std::runtime_error("Cbar has fewer rows than the queried spherical expansion degree");
	}

	if (Sbar.n_rows < n_degree){
		throw std::runtime_error("Sbar has fewer rows than the queried spherical expansion degree");
	}


	arma::mat b_bar_real = arma::zeros<arma::mat>(n_max + 3,n_max + 3);
	arma::mat b_bar_imag = arma::zeros<arma::mat>(n_max + 3,n_max + 3);

	GetBnmNormalizedExterior(n_degree,
		b_bar_real,
		b_bar_imag,
		pos,
		ref_radius);


	double K0 = 0.5 * mu / ref_radius / ref_radius;
	double x_ddot = 0;
	double y_ddot = 0;
	double z_ddot = 0;


	for (unsigned int nn = 0; nn<=n_degree; nn++){

		double n = (double) nn;

		for (unsigned int mm = 0; mm<=nn; mm++){

			double m = (double) mm;
			double delta_1_m;

			if (mm == 1){
				delta_1_m = 1.0;
			}
			else{
				delta_1_m = 0.0;
			} 

			double K1 = sqrt( (n+2.0) * (n+1.0) * (2.0*n+1.0) / 2.0 / (2.0*n+3.0) );
			double K2 = sqrt( (n+m+2.0) * (n+m+1.0) * (2.0*n+1.0) / (2.0*n+3.0) );
			double K3 = sqrt( 2.0 * (n-m+2.0) * (n-m+1.0) * (2.0*n+1.0) / (2.0 - delta_1_m) / (2.0*n+3.0) );

			if (mm == 0){

				x_ddot -= 2.0*K0 * ( Cbar(nn,mm)*K1*b_bar_real(nn+1,mm+1) );
				y_ddot -= 2.0*K0 * ( Cbar(nn,mm)*K1*b_bar_imag(nn+1,mm+1) );
				z_ddot -= 2.0*K0 * ( Cbar(nn,mm)*sqrt((n-m+1.0)*(n+m+1.0)*(2.0*n+1.0)/(2.0*n+3.0))*b_bar_real(nn+1,mm) );

			}
			else{

				x_ddot += K0 * ( -Cbar(nn,mm)*K2*b_bar_real(nn+1,mm+1) -Sbar(nn,mm)*K2*b_bar_imag(nn+1,mm+1) +Cbar(nn,mm)*K3*b_bar_real(nn+1,mm-1) +Sbar(nn,mm)*K3*b_bar_imag(nn+1,mm-1));
				y_ddot += K0 * ( -Cbar(nn,mm)*K2*b_bar_imag(nn+1,mm+1) +Sbar(nn,mm)*K2*b_bar_real(nn+1,mm+1) -Cbar(nn,mm)*K3*b_bar_imag(nn+1,mm-1) +Sbar(nn,mm)*K3*b_bar_real(nn+1,mm-1));
				z_ddot -= 2.0*K0 * ( Cbar(nn,mm)*sqrt((n-m+1.0)*(n+m+1.0)*(2.0*n+1.0)/(2.0*n+3.0))*b_bar_real(nn+1,mm) +Sbar(nn,mm)*sqrt((n-m+1.0)*(n+m+1.0)*(2*n+1.0)/(2.0*n+3.0))*b_bar_imag(nn+1,mm) );

			} 

		} 

	} 

	arma::vec acceleration = {x_ddot,y_ddot,z_ddot};

	return acceleration;

} 




