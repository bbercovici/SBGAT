#include <SBGATPolyhedronGravityModelUQ.hpp>


double SBGATPolyhedronGravityModelUQ::GetPotentialVariance(double const * point) const{

}

double SBGATPolyhedronGravityModelUQ::GetPotentialVariance(const arma::vec::fixed<3> & point) const{

}

void SBGATPolyhedronGravityModelUQ::GetPotentialVarianceAccelerationCovariance(double const * point,double & potential_var, 
	arma::mat::fixed<3,3> & acc_cov) const{

}

void SBGATPolyhedronGravityModelUQ::GetPotentialVarianceAccelerationCovariance(const arma::vec::fixed<3> & point,double & potential_var, 
	arma::mat::fixed<3,3> & acc_cov) const{

}



arma::rowvec::fixed<10> SBGATPolyhedronGravityModelUQ::PartialUePartialXe(const arma::vec::fixed<3> & pos,const int & e) const{


	arma::rowvec::fixed<10> partial;

	double Le = this -> pgm_model -> GetLe(pos,e);
	arma::vec::fixed<3> r_ei_0 = this -> pgm_model -> GetRe(pos,e); 
	arma::mat::fixed<3,6> R_ei_0 = {
		{r_ei_0[0],0,0,r_ei_0[1],r_ei_0[2],0},
		{0,r_ei_0[1],0,r_ei_0[0],0,r_ei_0[2]},
		{0,0,r_ei_0[2],0,r_ei_0[0],r_ei_0[1]}
	};

	arma::vec::fixed<3> Ee_times_r_ei_0 = R_ei_0 * this -> pgm_model -> GetEeParam(e);


	partial(0) = arma::dot(r_ei_0,Ee_times_r_ei_0);
	partial.subvec(1,3) = 2 * Le * Ee_times_r_ei_0.t();
	partial.subvec(4,9) = Le * r_ei_0.t() * R_ei_0;

	return partial;



}


arma::rowvec::fixed<10> SBGATPolyhedronGravityModelUQ::PartialUfPartialXf(const arma::vec::fixed<3> & pos,
	const int & f) const{


	arma::rowvec::fixed<10> partial;

	double omega_f = this -> pgm_model -> GetOmegaf(pos,f);
	arma::vec::fixed<3> r_fi_0 = this -> pgm_model -> GetRf(pos,f); 
	arma::mat::fixed<3,6> R_fi_0 = {
		{r_fi_0[0],0,0,r_fi_0[1],r_fi_0[2],0},
		{0,r_fi_0[1],0,r_fi_0[0],0,r_fi_0[2]},
		{0,0,r_fi_0[2],0,r_fi_0[0],r_fi_0[1]}
	};

	arma::vec::fixed<3> F_times_r_fi_0 = R_fi_0 * this -> pgm_model -> GetFfParam(f);


	partial(0) = arma::dot(r_fi_0,F_times_r_fi_0);
	partial.subvec(1,3) = 2 * omega_f * F_times_r_fi_0.t();
	partial.subvec(4,9) = omega_f * r_fi_0.t() * R_fi_0;

	return partial;





}


arma::mat::fixed<10,9> SBGATPolyhedronGravityModelUQ::PartialXfPartialTf(const arma::vec::fixed<3> & pos, const int & f) const{

}


arma::rowvec::fixed<9> SBGATPolyhedronGravityModelUQ::PartialOmegafPartialTf(const arma::vec::fixed<3> & pos,const int & f) const{

}



arma::mat::fixed<6,9> SBGATPolyhedronGravityModelUQ::PartialFfPartialTf(const arma::vec::fixed<3> & pos,const int & f) const{

}



arma::mat::fixed<3,3> SBGATPolyhedronGravityModelUQ::PartialNormalizedVPartialNonNormalizedV(const arma::vec::fixed<3> & non_normalized_V){




	return (arma::eye<arma::mat>(3,3) / arma::norm(non_normalized_V) - non_normalized_V * non_normalized_V.t() / std::pow(arma::norm(non_normalized_V),3));



}


arma::mat::fixed<3,3> SBGATPolyhedronGravityModelUQ::PartialNfPartialTf(const arma::vec::fixed<3> & pos,const int & f) const{

}


arma::mat::fixed<3,3> SBGATPolyhedronGravityModelUQ::PartialFfPartialnf(const int & f) const{

}


arma::rowvec::fixed<6> SBGATPolyhedronGravityModelUQ::PartialLePartialAe(const arma::vec::fixed<3> & pos,const int & e) const{

}




arma::mat::fixed<3,6> SBGATPolyhedronGravityModelUQ::PartialRadiusEePartialAe(const arma::vec::fixed<3> & pos,const int & e) const{

}



arma::mat::fixed<6,6> SBGATPolyhedronGravityModelUQ::PartialEePartialAe(const arma::vec::fixed<3> & pos,const int & e) const{

}



arma::mat::fixed<6,9> SBGATPolyhedronGravityModelUQ::PartialEePartialTf(const arma::vec::fixed<3> & pos,const int & f) const{

}


arma::mat::fixed<10,24> SBGATPolyhedronGravityModelUQ::PartialXePartialBe(const arma::vec::fixed<3> & pos,const int & e) const{

}



arma::rowvec::fixed<6> SBGATPolyhedronGravityModelUQ::PartialEdgeLengthPartialAe(const arma::vec::fixed<3> & pos,const int & e) const{

}

arma::rowvec::fixed<24> SBGATPolyhedronGravityModelUQ::PartialEqrPartialBe(const arma::vec::fixed<3> & pos,const int & e,const int & q,const int & r) const{

}

arma::mat::fixed<6,24> SBGATPolyhedronGravityModelUQ::PartialEPartialBe(const arma::vec::fixed<3> & pos,const int & e,const int & q,const int & r) const{

}

void SBGATPolyhedronGravityModelUQ::TestPartials(double tol) const{


	this -> TestPartialUePartialXe(tol);
	this -> TestPartialUfPartialXf(tol);
	this -> TestPartialXfPartialTf(tol);
	this -> TestPartialOmegafPartialTf(tol);
	this -> TestPartialFfPartialTf(tol);
	this -> TestPartialNormalizedVPartialNonNormalizedV(tol);
	this -> TestPartialNfPartialTf(tol);
	this -> TestPartialFfPartialnf(tol);
	this -> TestPartialLePartialAe(tol);
	this -> TestPartialRadiusEePartialAe(tol);
	this -> TestPartialEePartialAe(tol);
	this -> TestPartialEePartialTf(tol);
	this -> TestPartialXePartialBe(tol);
	this -> TestPartialEdgeLengthPartialAe(tol);
	this -> TestPartialEqrPartialBe(tol);
	this -> TestPartialEPartialBe(tol);

}


void SBGATPolyhedronGravityModelUQ::TestPartialUePartialXe(double tol) const{

  	// Test point 
	arma::vec::fixed<3> pos = {3,2,3};

  	// Pick e
	int e = 1;

  	// Edge point
	arma::vec::fixed<3> rei_0 = this -> pgm_model -> GetRe(pos,e);

  	// Xe before
	arma::vec::fixed<10> Xe = this -> pgm_model -> GetXe(pos,e);

  	// Ue before
	double Ue = this -> pgm_model -> GetUe(Xe);

  	// Perturbation to Xe 
	arma::vec::fixed<10> dXe = 1e-3 * Xe;

  	// Xe after
	arma::vec::fixed<10> Xe_p = Xe + dXe;

  	// Ue after
	double Ue_p = this -> pgm_model -> GetUe(Xe_p);

  	// Non-linear difference
	double dUe = Ue_p - Ue;

  	// Partial 

	arma::rowvec::fixed<10> partial_Ue_partial_Xe = this -> PartialUePartialXe(pos,e);

  	// Linear difference
	double dUe_lin = arma::dot(partial_Ue_partial_Xe.t(), dXe);

	assert(std::abs(dUe - dUe_lin) / std::abs(dUe_lin) < tol);



}
void SBGATPolyhedronGravityModelUQ::TestPartialUfPartialXf(double tol) const{


	// Test point 
	arma::vec::fixed<3> pos = {3,2,3};

  	// Pick f
	int f = 1;

  	// Edge point
	arma::vec::fixed<3> rfi_0 = this -> pgm_model -> GetRf(pos,f);

  	// Xe before
	arma::vec::fixed<10> Xf = this -> pgm_model -> GetXf(pos,f);

  	// Ue before
	double Uf = this -> pgm_model -> GetUf(Xf);

  	// Perturbation to Xe 
	arma::vec::fixed<10> dXf = 1e-3 * Xf;

  	// Xe after
	arma::vec::fixed<10> Xf_p = Xf + dXf;

  	// Ue after
	double Uf_p = this -> pgm_model -> GetUf(Xf_p);

  	// Non-linear difference
	double dUf = Uf_p - Uf;

  	// Partial 

	arma::rowvec::fixed<10> partial_Uf_partial_Xf = this -> PartialUfPartialXf(pos,f);

  	// Linear difference
	double dUf_lin = arma::dot(partial_Uf_partial_Xf.t(), dXf);

	assert(std::abs(dUf - dUf_lin) / std::abs(dUf_lin) < tol);








}
void SBGATPolyhedronGravityModelUQ::TestPartialXfPartialTf(double tol) const{



}
void SBGATPolyhedronGravityModelUQ::TestPartialOmegafPartialTf(double tol) const{



}
void SBGATPolyhedronGravityModelUQ::TestPartialFfPartialTf(double tol) const{



}
void SBGATPolyhedronGravityModelUQ::TestPartialNormalizedVPartialNonNormalizedV(double tol) const{

  	// non-normalized vector
	arma::vec::fixed<3> X = {1,2,3};


  	// normalized vector
	arma::vec::fixed<3> x = arma::normalise(X);


  	// disturbance 
	arma::vec::fixed<3> dX = {1e-2,2e-2,-3e-2};


  	// normalized vector, perturbed
	arma::vec::fixed<3> x_p = arma::normalise(X + dX);

  	// difference, non-linear
	arma::vec::fixed<3> dx = x_p - x;

  	// difference, linear
	arma::vec::fixed<3> dx_lin = this -> PartialNormalizedVPartialNonNormalizedV(X) * dX;


	assert(arma::norm(x_p - x - dx_lin)/arma::norm(dx_lin) < tol);



}
void SBGATPolyhedronGravityModelUQ::TestPartialNfPartialTf(double tol) const{



}
void SBGATPolyhedronGravityModelUQ::TestPartialFfPartialnf(double tol) const{



}
void SBGATPolyhedronGravityModelUQ::TestPartialLePartialAe(double tol) const{



}
void SBGATPolyhedronGravityModelUQ::TestPartialRadiusEePartialAe(double tol) const{



}
void SBGATPolyhedronGravityModelUQ::TestPartialEePartialAe(double tol) const{



}
void SBGATPolyhedronGravityModelUQ::TestPartialEePartialTf(double tol) const{



}
void SBGATPolyhedronGravityModelUQ::TestPartialXePartialBe(double tol) const{



}
void SBGATPolyhedronGravityModelUQ::TestPartialEdgeLengthPartialAe(double tol) const{



}
void SBGATPolyhedronGravityModelUQ::TestPartialEqrPartialBe(double tol) const{



}
void SBGATPolyhedronGravityModelUQ::TestPartialEPartialBe(double tol) const{



}






