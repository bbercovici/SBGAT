#include <SBGATPolyhedronGravityModelUQ.hpp>
#include <RigidBodyKinematics.hpp>
#include <vtkOBJReader.h>
#include <vtkCleanPolyData.h>


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




	arma::mat::fixed<10,9> partial;

	partial.row(0) = this -> PartialOmegafPartialTf(pos,f);

	partial.rows(1,3) = this -> PartialRadiusFfPartialTf();

	partial.rows(4,9) = this -> PartialFfPartialTf(f);

	return partial;







}


arma::mat::fixed<3,6> SBGATPolyhedronGravityModelUQ::PartialRadiusEePartialAe() const{

	arma::mat::fixed<3,6> partial = arma::zeros<arma::mat>(3,6);
	partial.submat(0,0,2,2) = arma::eye<arma::mat>(3,3);
	return partial;
}

arma::mat::fixed<3,9> SBGATPolyhedronGravityModelUQ::PartialRadiusFfPartialTf() const{

	arma::mat::fixed<3,9> partial = arma::zeros<arma::mat>(3,9);
	partial.submat(0,0,2,2) = arma::eye<arma::mat>(3,3);
	return partial;
}





arma::rowvec::fixed<9> SBGATPolyhedronGravityModelUQ::PartialOmegafPartialTf(const arma::vec::fixed<3> & pos,const int & f) const{


	double r0[3],r1[3],r2[3];

	this -> pgm_model -> GetVerticesInFacet(f,r0,r1,r2);

	arma::vec::fixed<3> R0 = {r0[0] - pos[0],r0[1] - pos[1],r0[2] - pos[2]};
	arma::vec::fixed<3> R1 = {r1[0] - pos[0],r1[1] - pos[1],r1[2] - pos[2]};
	arma::vec::fixed<3> R2 = {r2[0] - pos[0],r2[1] - pos[1],r2[2] - pos[2]};

	arma::vec::fixed<3> r0_hat = arma::normalise(R0);
	arma::vec::fixed<3> r1_hat = arma::normalise(R1);
	arma::vec::fixed<3> r2_hat = arma::normalise(R2);

	arma::vec::fixed<9> r_hat = {
		r0_hat(0),r0_hat(1),r0_hat(2),
		r1_hat(0),r1_hat(1),r1_hat(2),
		r2_hat(0),r2_hat(1),r2_hat(2)
	};

	double alpha = 1 + arma::dot(r0_hat,r1_hat) + arma::dot(r0_hat,r2_hat) + arma::dot(r1_hat,r2_hat);
	double gamma = arma::dot(r0_hat,arma::cross(r1_hat,r2_hat));

	arma::vec::fixed<2> Zf = {alpha,gamma};

	arma::mat::fixed<9,9> partial_r_normalized_partial_R = arma::zeros<arma::mat>(9,9);

	partial_r_normalized_partial_R.submat(0,0,2,2) = SBGATPolyhedronGravityModelUQ::PartialNormalizedVPartialNonNormalizedV(R0);
	partial_r_normalized_partial_R.submat(3,3,5,5) = SBGATPolyhedronGravityModelUQ::PartialNormalizedVPartialNonNormalizedV(R1);
	partial_r_normalized_partial_R.submat(6,6,8,8) = SBGATPolyhedronGravityModelUQ::PartialNormalizedVPartialNonNormalizedV(R2);

	return ( 2 * SBGATPolyhedronGravityModelUQ::PartialAtan2PartialZf(Zf) * SBGATPolyhedronGravityModelUQ::PartialZfPartialUnitRf(r_hat) * partial_r_normalized_partial_R);

}



arma::mat::fixed<2,9> SBGATPolyhedronGravityModelUQ::PartialZfPartialUnitRf(const arma::vec::fixed<9> & UnitRf){


	arma::mat::fixed<9,2> partial;

	const arma::vec::fixed<3> & r0_hat = UnitRf.subvec(0,2);
	const arma::vec::fixed<3> & r1_hat = UnitRf.subvec(3,5);
	const arma::vec::fixed<3> & r2_hat = UnitRf.subvec(6,8);

	double dot_prod = arma::dot(r0_hat,arma::cross(r1_hat,r2_hat));

	partial.submat(0,0,2,0) = r1_hat + r2_hat;
	partial.submat(3,0,5,0) = r0_hat + r2_hat;
	partial.submat(6,0,8,0) = r0_hat + r1_hat;


	partial.submat(0,1,2,1) = RBK::tilde(r1_hat) * r2_hat;
	partial.submat(3,1,5,1) = RBK::tilde(r2_hat) * r0_hat;
	partial.submat(6,1,8,1) = RBK::tilde(r0_hat) * r1_hat;


	return partial.t();



}


arma::rowvec::fixed<2> SBGATPolyhedronGravityModelUQ::PartialAtan2PartialZf(const arma::vec::fixed<2> & Zf){

	arma::vec::fixed<2> e1 = {1,0};
	arma::vec::fixed<2> e2 = {0,1};

	arma::rowvec::fixed<2> partial = 2 * (
		((arma::norm(Zf) + arma::dot(e1,Zf)) * e2.t() - arma::dot(e2,Zf) * (e1.t() + Zf.t()/arma::norm(Zf)))
		/(std::pow(arma::norm(Zf) + arma::dot(e1,Zf),2) + std::pow(arma::dot(e2,Zf),2)));
	return partial;

}








arma::mat::fixed<6,9> SBGATPolyhedronGravityModelUQ::PartialFfPartialTf(const int & f) const{




	double r0[3],r1[3],r2[3];

	this -> pgm_model -> GetVerticesInFacet(f,r0,r1,r2);


	arma::vec::fixed<3> R0 = {r0[0],r0[1],r0[2]};
	arma::vec::fixed<3> R1 = {r1[0],r1[1],r1[2]};
	arma::vec::fixed<3> R2 = {r2[0],r2[1],r2[2]};


	return this -> PartialFfPartialnf(f) * SBGATPolyhedronGravityModelUQ::PartialNormalizedVPartialNonNormalizedV(arma::cross(R1 - R0,R2 - R0)) * this -> PartialNfPartialTf(f);



}



arma::mat::fixed<3,3> SBGATPolyhedronGravityModelUQ::PartialNormalizedVPartialNonNormalizedV(const arma::vec::fixed<3> & non_normalized_V){




	return (arma::eye<arma::mat>(3,3) / arma::norm(non_normalized_V) - non_normalized_V * non_normalized_V.t() / std::pow(arma::norm(non_normalized_V),3));



}


arma::mat::fixed<3,9> SBGATPolyhedronGravityModelUQ::PartialNfPartialTf(const int & f) const{


	arma::mat::fixed<3,9> partial;

	double r0[3], r1[3], r2[3];

	this -> pgm_model -> GetVerticesInFacet(f,r0,r1,r2);

	arma::vec::fixed<3> r0_arma = {r0[0],r0[1],r0[2]};
	arma::vec::fixed<3> r1_arma = {r1[0],r1[1],r1[2]};
	arma::vec::fixed<3> r2_arma = {r2[0],r2[1],r2[2]};

	partial.cols(0,2) = RBK::tilde(r2_arma - r1_arma);
	partial.cols(3,5) = RBK::tilde(r0_arma - r2_arma);
	partial.cols(6,8) = RBK::tilde(r1_arma - r0_arma);

	return partial;


}


arma::mat::fixed<6,3> SBGATPolyhedronGravityModelUQ::PartialFfPartialnf(const int & f) const{



	double n[3];
	this -> pgm_model -> GetFacetNormal(f,n);
	arma::vec::fixed<3> e0 = {1,0,0};
	arma::vec::fixed<3> e1 = {0,1,0};
	arma::vec::fixed<3> e2 = {0,0,1};


	arma::vec n_arma = {n[0],n[1],n[2]};

	arma::mat::fixed<6,3> partial;

	partial.row(0) = n_arma.t() * e0*e0.t();
	partial.row(1) = n_arma.t() * e1*e1.t();
	partial.row(2) = n_arma.t() * e2*e2.t();
	partial.row(3) = n_arma.t() * e0*e1.t();
	partial.row(4) = n_arma.t() * e0*e2.t();
	partial.row(5) = n_arma.t() * e1*e2.t();

	return partial;


}


arma::rowvec::fixed<6> SBGATPolyhedronGravityModelUQ::PartialLePartialAe(const arma::vec::fixed<3> & pos,const int & e) const{

	double r0[3],r1[3];

	this -> pgm_model -> GetVerticesOnEdge(e,r0,r1);

	arma::vec::fixed<3> r0_arma = {r0[0],r0[1],r0[2]};
	arma::vec::fixed<3> r1_arma = {r1[0],r1[1],r1[2]};


	r0_arma -= pos;
	r1_arma -= pos;

	double re_0 = arma::norm(r0_arma);
	double re_1 = arma::norm(r1_arma);
	double le = arma::norm(r1_arma - r0_arma);

	double beta_plus = re_0 + re_1 + le;
	double beta_minus = re_0 + re_1 - le;

	arma::vec::fixed<3> beta_vectors = {
		1./beta_plus - 1./beta_minus,
		1./beta_plus + 1./beta_minus,
		1./beta_plus + 1./beta_minus
	};


	arma::mat::fixed<3,6> partial_re_0_partial_Ae = arma::zeros<arma::mat>(3,6);
	arma::mat::fixed<3,6> partial_re_1_partial_Ae = arma::zeros<arma::mat>(3,6);

	partial_re_0_partial_Ae.cols(0,2) = arma::eye<arma::mat>(3,3);
	partial_re_1_partial_Ae.cols(3,5) = arma::eye<arma::mat>(3,3);


	arma::mat::fixed<3,6> partial_mat;
	partial_mat.row(0) = this -> PartialEdgeLengthPartialAe(e);
	partial_mat.row(1) = r0_arma.t() / re_0 * partial_re_0_partial_Ae;
	partial_mat.row(2) = r1_arma.t() / re_1 * partial_re_1_partial_Ae;


	return beta_vectors.t() * partial_mat;

}




arma::mat::fixed<10,24> SBGATPolyhedronGravityModelUQ::PartialXePartialBe(const arma::vec::fixed<3> & pos,const int & e) const{


	arma::mat::fixed<10,24> partial = arma::zeros<arma::mat>(10,24);

	partial.submat(0,0,0,5) = this -> PartialLePartialAe(pos,e);
	partial.submat(1,0,3,5) = this -> PartialRadiusEePartialAe();

	partial.submat(4,0,9,23) = this -> PartialEPartialBe(e);


	return partial;



}



arma::rowvec::fixed<6> SBGATPolyhedronGravityModelUQ::PartialEdgeLengthPartialAe(const int & e) const{


	arma::rowvec::fixed<6> partial;

	double r0[3], r1[3];
	this -> pgm_model -> GetVerticesOnEdge(e,r0,r1);

	arma::vec::fixed<3> r0_arma = {r0[0],r0[1],r0[2]};
	arma::vec::fixed<3> r1_arma = {r1[0],r1[1],r1[2]};

	double le = arma::norm(r1_arma - r0_arma);

	arma::mat::fixed<3,6> mat;
	mat.cols(0,2) = arma::eye<arma::mat>(3,3);
	mat.cols(3,5) = - arma::eye<arma::mat>(3,3);

	return (r0_arma - r1_arma).t()/le * mat;


}

arma::rowvec::fixed<24> SBGATPolyhedronGravityModelUQ::PartialEqrPartialBe(const int & e,const int & q,const int & r) const{

	arma::rowvec::fixed<24> partial;
	arma::vec::fixed<3> e_q = arma::zeros<arma::vec>(3);
	arma::vec::fixed<3> e_r = arma::zeros<arma::vec>(3);
	e_q(q) = 1;
	e_r(r) = 1;

	double r0[3],r1[3];


	this -> pgm_model -> GetVerticesOnEdge(e,r0,r1);

	arma::vec::fixed<3> r0_arma = {r0[0],r0[1],r0[2]};
	arma::vec::fixed<3> r1_arma = {r1[0],r1[1],r1[2]};

	int f0,f1;

	this -> pgm_model -> GetIndicesOfAdjacentFacets(e,f0,f1);

	double r0_e_f0[3],r1_e_f0[3],r2_e_f0[3];
	double r0_e_f1[3],r1_e_f1[3],r2_e_f1[3];

	this -> pgm_model -> GetVerticesInFacet(f0,r0_e_f0,r1_e_f0,r2_e_f0);
	this -> pgm_model -> GetVerticesInFacet(f1,r0_e_f1,r1_e_f1,r2_e_f1);

	arma::vec::fixed<3> r0_e_f0_arma = {r0_e_f0[0],r0_e_f0[1],r0_e_f0[2]};
	arma::vec::fixed<3> r1_e_f0_arma = {r1_e_f0[0],r1_e_f0[1],r1_e_f0[2]};
	arma::vec::fixed<3> r2_e_f0_arma = {r2_e_f0[0],r2_e_f0[1],r2_e_f0[2]};

	arma::vec::fixed<3> r0_e_f1_arma = {r0_e_f1[0],r0_e_f1[1],r0_e_f1[2]};
	arma::vec::fixed<3> r1_e_f1_arma = {r1_e_f1[0],r1_e_f1[1],r1_e_f1[2]};
	arma::vec::fixed<3> r2_e_f1_arma = {r2_e_f1[0],r2_e_f1[1],r2_e_f1[2]};


	arma::vec::fixed<3> N_e_f0 = arma::cross(r1_e_f0_arma - r0_e_f0_arma,r2_e_f0_arma - r0_e_f0_arma);
	arma::vec::fixed<3> N_e_f1 = arma::cross(r1_e_f1_arma - r0_e_f0_arma,r2_e_f1_arma - r0_e_f1_arma);

	arma::vec::fixed<3> n_e_f0 = arma::normalise(N_e_f0);
	arma::vec::fixed<3> n_e_f1 = arma::normalise(N_e_f1);

	double E_qr = arma::dot(e_q,(n_e_f0 * n_e_f0.t() - n_e_f1 * n_e_f1.t()) * RBK::tilde(r0_arma - r1_arma) * e_r);


	double le = arma::norm(r0_arma - r1_arma);



	partial.subvec(0,5) = 1./le * (E_qr * this -> PartialEdgeLengthPartialAe(e)
		- e_q.t() * (n_e_f0 * n_e_f0.t() - n_e_f1 * n_e_f1.t()) * RBK::tilde(e_r));

	partial.subvec(6,14) = 1./le * (n_e_f0.t() * RBK::tilde(r0_arma - r1_arma) * e_r * e_q.t()
		- arma::dot(e_q, n_e_f0) * e_r.t() * RBK::tilde(r0_arma - r1_arma)) * this -> PartialNormalizedVPartialNonNormalizedV(N_e_f0) * this -> PartialNfPartialTf(f0) ;

	partial.subvec(15,23) = 1./le * (n_e_f1.t() * RBK::tilde(r0_arma - r1_arma) * e_r * e_q.t()
		- arma::dot(e_q, n_e_f1) * e_r.t() * RBK::tilde(r0_arma - r1_arma))* this -> PartialNormalizedVPartialNonNormalizedV(N_e_f1) * this -> PartialNfPartialTf(f1) ;



	return partial;

}

arma::mat::fixed<6,24> SBGATPolyhedronGravityModelUQ::PartialEPartialBe(const int & e) const{

	arma::mat::fixed<6,24> partial;

	partial.row(0) = this -> PartialEqrPartialBe(e,0,0);
	partial.row(1) = this -> PartialEqrPartialBe(e,1,1);
	partial.row(2) = this -> PartialEqrPartialBe(e,2,2);
	partial.row(3) = this -> PartialEqrPartialBe(e,0,1);
	partial.row(4) = this -> PartialEqrPartialBe(e,0,2);
	partial.row(5) = this -> PartialEqrPartialBe(e,1,2);


	return partial;


}

void SBGATPolyhedronGravityModelUQ::TestPartials(double tol){

	SBGATPolyhedronGravityModelUQ::TestPartialAtan2PartialZf(tol);
	SBGATPolyhedronGravityModelUQ::TestPartialZfPartialUnitRf(tol);
	SBGATPolyhedronGravityModelUQ::TestPartialOmegafPartialTf(tol);
	SBGATPolyhedronGravityModelUQ::TestPartialEdgeLengthPartialAe(tol);
	SBGATPolyhedronGravityModelUQ::TestPartialNfPartialTf(tol);
	SBGATPolyhedronGravityModelUQ::TestPartialUePartialXe(tol);
	SBGATPolyhedronGravityModelUQ::TestPartialUfPartialXf(tol);
	SBGATPolyhedronGravityModelUQ::TestPartialXfPartialTf(tol);
	// SBGATPolyhedronGravityModelUQ::TestPartialFfPartialTf(tol);
	// SBGATPolyhedronGravityModelUQ::TestPartialNormalizedVPartialNonNormalizedV(tol);
	// SBGATPolyhedronGravityModelUQ::TestPartialNfPartialTf(tol);
	// SBGATPolyhedronGravityModelUQ::TestPartialFfPartialnf(tol);
	// SBGATPolyhedronGravityModelUQ::TestPartialLePartialAe(tol);
	// SBGATPolyhedronGravityModelUQ::TestPartialEePartialAe(tol);
	// SBGATPolyhedronGravityModelUQ::TestPartialEePartialTf(tol);
	// SBGATPolyhedronGravityModelUQ::TestPartialXePartialBe(tol);
	// SBGATPolyhedronGravityModelUQ::TestPartialEqrPartialBe(tol);
	// SBGATPolyhedronGravityModelUQ::TestPartialEPartialBe(tol);

}


void SBGATPolyhedronGravityModelUQ::TestPartialUePartialXe(double tol){

	std::cout << "\t In TestPartialUePartialXe ...\n";

	std::string filename  = "../input/cube.obj";

	// Reading
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName(filename.c_str());
	reader -> Update(); 

	vtkSmartPointer<vtkCleanPolyData> cleaner =
	vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner -> SetInputConnection (reader -> GetOutputPort());
	cleaner -> Update();

	// Creating the PGM dyads
	vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
	pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
	pgm_filter -> SetDensity(1970); // density in kg/m^3
	pgm_filter -> SetScaleMeters();
	pgm_filter -> Update();

	SBGATPolyhedronGravityModelUQ shape_uq;
	shape_uq.SetPGM(pgm_filter);


  	// Test point 
	arma::vec::fixed<3> pos = {3,2,3};

  	// Pick e
	int e = 1;

  	// Edge point
	arma::vec::fixed<3> rei_0 = shape_uq.GetPGMModel() -> GetRe(pos,e);

  	// Xe before
	arma::vec::fixed<10> Xe = shape_uq.GetPGMModel() -> GetXe(pos,e);

  	// Ue before
	double Ue = shape_uq.GetPGMModel() -> GetUe(Xe);

  	// Perturbation to Xe 
	arma::vec::fixed<10> dXe = 1e-3 * Xe;

  	// Xe after
	arma::vec::fixed<10> Xe_p = Xe + dXe;

  	// Ue after
	double Ue_p = shape_uq.GetPGMModel() -> GetUe(Xe_p);

  	// Non-linear difference
	double dUe = Ue_p - Ue;

  	// Partial 

	arma::rowvec::fixed<10> partial_Ue_partial_Xe = shape_uq.PartialUePartialXe(pos,e);

  	// Linear difference
	double dUe_lin = arma::dot(partial_Ue_partial_Xe.t(), dXe);

	assert(std::abs(dUe - dUe_lin) / std::abs(dUe_lin) < tol);

	std::cout << "\t Passed TestPartialUePartialXe ...\n";


}
void SBGATPolyhedronGravityModelUQ::TestPartialUfPartialXf(double tol){


	std::cout << "\t In TestPartialUfPartialXf ...\n";

	std::string filename  = "../input/cube.obj";

	// Reading
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName(filename.c_str());
	reader -> Update(); 

	// Creating the PGM dyads
	vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
	pgm_filter -> SetInputConnection(reader -> GetOutputPort());
	pgm_filter -> SetDensity(1970); // density in kg/m^3
	pgm_filter -> SetScaleMeters();
	pgm_filter -> Update();

	SBGATPolyhedronGravityModelUQ shape_uq;
	shape_uq.SetPGM(pgm_filter);

	// Test point 
	arma::vec::fixed<3> pos = {3,2,3};

  	// Pick f
	int f = 1;

  	// Edge point
	arma::vec::fixed<3> rfi_0 = shape_uq.GetPGMModel() -> GetRf(pos,f);

  	// Xe before
	arma::vec::fixed<10> Xf = shape_uq.GetPGMModel() -> GetXf(pos,f);

  	// Ue before
	double Uf = shape_uq.GetPGMModel() -> GetUf(Xf);

  	// Perturbation to Xe 
	arma::vec::fixed<10> dXf = 1e-3 * Xf;

  	// Xe after
	arma::vec::fixed<10> Xf_p = Xf + dXf;

  	// Ue after
	double Uf_p = shape_uq.GetPGMModel() -> GetUf(Xf_p);

  	// Non-linear difference
	double dUf = Uf_p - Uf;

  	// Partial 

	arma::rowvec::fixed<10> partial_Uf_partial_Xf = shape_uq. PartialUfPartialXf(pos,f);

  	// Linear difference
	double dUf_lin = arma::dot(partial_Uf_partial_Xf.t(), dXf);

	assert(std::abs(dUf - dUf_lin) / std::abs(dUf_lin) < tol);

	std::cout << "\t Passed TestPartialUfPartialXf ...\n";


}
void SBGATPolyhedronGravityModelUQ::TestPartialXfPartialTf(double tol){


	std::cout << "\t In TestPartialXfPartialTf ...\n";

	std::string filename  = "../input/cube.obj";

	// Reading
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName(filename.c_str());
	reader -> Update(); 

	// Creating the PGM dyads
	vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
	pgm_filter -> SetInputConnection(reader -> GetOutputPort());
	pgm_filter -> SetDensity(1970); // density in kg/m^3
	pgm_filter -> SetScaleMeters();
	pgm_filter -> Update();

	SBGATPolyhedronGravityModelUQ shape_uq;
	shape_uq.SetPGM(pgm_filter);

	int f = 1;
	arma::vec pos = {2,3,4};

	// Nominal Xf
	arma::vec::fixed<10> Xf;
	Xf(0) = pgm_filter -> GetOmegaf(pos,f);
	Xf.subvec(1,3) = pgm_filter -> GetRf(pos,f);
	Xf.subvec(4,9) = pgm_filter -> GetFfParam(f);

	// Deviation
	arma::vec::fixed<9> delta_Tf = {1e-2,-1e-2,1e-2,2e-2,-1e-2,1e-2,1e-2,1e-2,-1e-2};

	// Apply Tf deviation
	shape_uq. ApplyTfDeviation(delta_Tf,f);

	// Perturbed Xf
	arma::vec::fixed<10> Xf_p;
	Xf_p(0) = pgm_filter -> GetOmegaf(pos,f);
	Xf_p.subvec(1,3) = pgm_filter -> GetRf(pos,f);
	Xf_p.subvec(4,9) = pgm_filter -> GetFfParam(f);

	// Non-linear dXf
	arma::vec::fixed<10> dXf = Xf_p - Xf;

	// Linear dXf
	arma::vec::fixed<10> dXf_lin = shape_uq.PartialXfPartialTf(pos,f) * delta_Tf;

	assert(arma::norm(dXf - dXf_lin) / arma::norm(dXf_lin) < 1e-2);

















	std::cout << "\t Passed TestPartialXfPartialTf ...\n";

}
void SBGATPolyhedronGravityModelUQ::TestPartialOmegafPartialTf(double tol){




	std::cout << "\t In TestPartialOmegafPartialTf ...\n";

	std::string filename  = "../input/cube.obj";

	// Reading
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName(filename.c_str());
	reader -> Update(); 

	// Creating the PGM dyads
	vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
	pgm_filter -> SetInputConnection(reader -> GetOutputPort());
	pgm_filter -> SetDensity(1970); // density in kg/m^3
	pgm_filter -> SetScaleMeters();
	pgm_filter -> Update();

	SBGATPolyhedronGravityModelUQ shape_uq;
	shape_uq.SetPGM(pgm_filter);

	int f = 1;
	arma::vec pos = {2,3,4};

	// Nominal omega_f
	double omega_f = pgm_filter -> GetOmegaf(pos,f);

	// Deviation
	arma::vec::fixed<9> delta_Tf = {1e-2,1e-2,1e-2,1e-2,-1e-2,1e-2,1e-2,1e-2,-1e-2};

	// Apply Tf deviation
	shape_uq. ApplyTfDeviation(delta_Tf,f);

	// Perturbed omega_f
	double omega_f_p = pgm_filter -> GetOmegaf(pos,f);
	
	// Non-linear domega_f
	double domega_f = omega_f_p - omega_f;

	// Linear dXf
	double domega_f_lin = arma::dot(shape_uq.PartialOmegafPartialTf(pos,f), delta_Tf);

	std::cout << "linear : " << domega_f_lin << std::endl;
	std::cout << "non-linear : " << domega_f << std::endl;


	assert(std::abs(domega_f - domega_f_lin) / std::abs(domega_f_lin) < 1e-2);

	std::cout << "\t Passed TestPartialOmegafPartialTf ...\n";





}
void SBGATPolyhedronGravityModelUQ::TestPartialFfPartialTf(double tol){



}


void SBGATPolyhedronGravityModelUQ::TestPartialNormalizedVPartialNonNormalizedV(double tol){

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
	arma::vec::fixed<3> dx_lin = SBGATPolyhedronGravityModelUQ::PartialNormalizedVPartialNonNormalizedV(X) * dX;


	assert(arma::norm(x_p - x - dx_lin)/arma::norm(dx_lin) < tol);



}
void SBGATPolyhedronGravityModelUQ::TestPartialNfPartialTf(double tol){

	std::cout << "\t In TestPartialNfPartialTf ...\n";

	std::string filename  = "../input/cube.obj";

	// Reading
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName(filename.c_str());
	reader -> Update(); 

	// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cleaner =
	vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner -> SetInputConnection (reader -> GetOutputPort());
	cleaner -> Update();

	// Creating the PGM dyads
	vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
	pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
	pgm_filter -> SetDensity(1970); // density in kg/m^3
	pgm_filter -> SetScaleMeters();
	pgm_filter -> Update();


	SBGATPolyhedronGravityModelUQ shape_uq;
	shape_uq.SetPGM(pgm_filter);


	int f = 1;
	arma::vec::fixed<9> delta_Tf = {1e-2,-1e-2,1e-2,1e-2,-2e-2,3e-2,1e-2,-1e-2,1e-2};

	double r0[3],r1[3],r2[3];

	shape_uq.GetPGMModel() -> GetVerticesInFacet(f,r0,r1,r2);

	arma::vec::fixed<3> r0_arma = {r0[0],r0[1],r0[2]};
	arma::vec::fixed<3> r1_arma = {r1[0],r1[1],r1[2]};
	arma::vec::fixed<3> r2_arma = {r2[0],r2[1],r2[2]};

	arma::vec::fixed<3> N = arma::cross(r1_arma - r0_arma,r2_arma - r0_arma);

	// Apply Tf deviation
	shape_uq. ApplyTfDeviation(delta_Tf,f);

	// Perturbed normal
	shape_uq.GetPGMModel() -> GetVerticesInFacet(f,r0,r1,r2);

	arma::vec::fixed<3> r0_arma_p = {r0[0],r0[1],r0[2]};
	arma::vec::fixed<3> r1_arma_p = {r1[0],r1[1],r1[2]};
	arma::vec::fixed<3> r2_arma_p = {r2[0],r2[1],r2[2]};

	assert(arma::norm(r0_arma_p - r0_arma - delta_Tf.subvec(0,2))/arma::abs(delta_Tf.subvec(0,2)).max() < 1e-4);
	assert(arma::norm(r1_arma_p - r1_arma - delta_Tf.subvec(3,5))/arma::abs(delta_Tf.subvec(3,5)).max() < 1e-4);
	assert(arma::norm(r2_arma_p - r2_arma - delta_Tf.subvec(6,8))/arma::abs(delta_Tf.subvec(6,8)).max() < 1e-4);

	arma::vec::fixed<3> N_p = arma::cross(r1_arma_p - r0_arma_p,r2_arma_p - r0_arma_p);

	// Non-linear difference
	arma::vec::fixed<3> dN = N_p - N;

	// Linear difference
	arma::vec::fixed<3> dN_lin = shape_uq . PartialNfPartialTf(f) * delta_Tf;


	assert(arma::abs(dN - dN_lin).max()/arma::abs(dN_lin).max() < tol);
	

	std::cout << "\t Passed TestPartialNfPartialTf ...\n";



}
void SBGATPolyhedronGravityModelUQ::TestPartialFfPartialnf(double tol){



}
void SBGATPolyhedronGravityModelUQ::TestPartialLePartialAe(double tol){



}

void SBGATPolyhedronGravityModelUQ::TestPartialEePartialAe(double tol){



}
void SBGATPolyhedronGravityModelUQ::TestPartialEePartialTf(double tol){



}
void SBGATPolyhedronGravityModelUQ::TestPartialXePartialBe(double tol){


}

void SBGATPolyhedronGravityModelUQ::TestPartialEdgeLengthPartialAe(double tol){

	std::cout << "\t In TestPartialEdgeLengthPartialAe ...\n";
	std::string filename  = "../input/cube.obj";

	// Reading
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName(filename.c_str());
	reader -> Update(); 

	// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cleaner =
	vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner -> SetInputConnection (reader -> GetOutputPort());
	cleaner -> Update();

	// Creating the PGM dyads
	vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
	pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
	pgm_filter -> SetDensity(1970); // density in kg/m^3
	pgm_filter -> SetScaleMeters();
	pgm_filter -> Update();

	SBGATPolyhedronGravityModelUQ shape_uq;
	shape_uq.SetPGM(pgm_filter);

	int e = 1;
	arma::vec::fixed<6> delta_Ae = {1e-2,-1e-3,1e-3,-1e-2,1e-3,-1e-2};

	double r0[3],r1[3];

	shape_uq.GetPGMModel() -> GetVerticesOnEdge(e,r0,r1);
	
	double le = std::sqrt(vtkMath::Distance2BetweenPoints(r0,r1));

	arma::vec::fixed<3> r0_arma= {r0[0],r0[1],r0[2]};
	arma::vec::fixed<3> r1_arma= {r1[0],r1[1],r1[2]};


	// Apply Ae deviation
	shape_uq.ApplyAeDeviation(delta_Ae,e);

	shape_uq.GetPGMModel() -> GetVerticesOnEdge(e,r0,r1);

	arma::vec::fixed<3> r0_arma_p = {r0[0],r0[1],r0[2]};
	arma::vec::fixed<3> r1_arma_p = {r1[0],r1[1],r1[2]};
	assert(arma::norm(r0_arma_p - r0_arma - delta_Ae.subvec(0,2))/arma::abs(delta_Ae.subvec(0,2)).max()<1e-3);
	assert(arma::norm(r1_arma_p - r1_arma - delta_Ae.subvec(3,5))/arma::abs(delta_Ae.subvec(3,5)).max()<1e-3);

	
	// Perturbed length
	double le_p = std::sqrt(vtkMath::Distance2BetweenPoints(r0,r1));


	// Non-linear difference
	double dle = le_p - le;

	// Linear difference
	double dle_lin = arma::dot(shape_uq.PartialEdgeLengthPartialAe(e), delta_Ae);


	assert(std::abs(dle - dle_lin)/std::abs(dle_lin) < tol);
	
	std::cout << "\t Passed TestPartialEdgeLengthPartialAe ...\n";

}
void SBGATPolyhedronGravityModelUQ::TestPartialEqrPartialBe(double tol){



}
void SBGATPolyhedronGravityModelUQ::TestPartialEPartialBe(double tol){



}




void SBGATPolyhedronGravityModelUQ::ApplyAeDeviation(arma::vec::fixed<6> delta_Ae,const int & e){
	
	int v0_index,v1_index;

	this -> pgm_model -> GetIndicesVerticesOnEdge(e,v0_index,v1_index);

	double r0[3],r1[3];
	vtkPolyData * polydata = vtkPolyData::SafeDownCast(this -> pgm_model -> GetInput());

	polydata -> GetPoint(v0_index,r0);
	polydata -> GetPoint(v1_index,r1);


	r0[0] += delta_Ae(0);
	r0[1] += delta_Ae(1);
	r0[2] += delta_Ae(2);

	r1[0] += delta_Ae(3);
	r1[1] += delta_Ae(4);
	r1[2] += delta_Ae(5);

	
	polydata -> GetPoints() -> SetPoint(v0_index,r0);
	polydata -> GetPoints() -> SetPoint(v1_index,r1);


	polydata -> GetPoints() -> Modified();

	polydata -> Modified();


	this -> pgm_model -> Modified();

	this -> pgm_model -> Update();


}


void SBGATPolyhedronGravityModelUQ::ApplyTfDeviation(arma::vec::fixed<9> delta_Tf,const int & f){
	
	int v0_index,v1_index,v2_index;
	this -> pgm_model -> GetIndicesVerticesInFacet(f,v0_index,v1_index,v2_index);

	double r0[3],r1[3],r2[3];
	vtkPolyData * polydata = vtkPolyData::SafeDownCast(this -> pgm_model -> GetInput());

	polydata -> GetPoint(v0_index,r0);
	polydata -> GetPoint(v1_index,r1);
	polydata -> GetPoint(v2_index,r2);


	r0[0] += delta_Tf(0);
	r0[1] += delta_Tf(1);
	r0[2] += delta_Tf(2);

	r1[0] += delta_Tf(3);
	r1[1] += delta_Tf(4);
	r1[2] += delta_Tf(5);

	r2[0] += delta_Tf(6);
	r2[1] += delta_Tf(7);
	r2[2] += delta_Tf(8);

	polydata -> GetPoints() -> SetPoint(v0_index,r0);
	polydata -> GetPoints() -> SetPoint(v1_index,r1);
	polydata -> GetPoints() -> SetPoint(v2_index,r2);

	polydata -> GetPoints() -> Modified();

	polydata -> Modified();
	this -> pgm_model -> Modified();

	this -> pgm_model -> Update();

}


void SBGATPolyhedronGravityModelUQ::TestPartialAtan2PartialZf(double tol){

	std::cout << "\t In TestPartialAtan2PartialZf ...\n";

	arma::vec::fixed<2> e1 = {1,0};
	arma::vec::fixed<2> e2 = {0,1};

	arma::vec::fixed<2> Zf = {0.1,-0.2};

	double at2 = 2 * std::atan(arma::dot(Zf,e2) / (arma::norm(Zf) + arma::dot(Zf,e1)));


	arma::vec::fixed<2> dZf = {1e-3,-1e-3};

	double at2_p = 2 * std::atan(arma::dot(Zf + dZf,e2) / (arma::norm(Zf + dZf) + arma::dot(Zf + dZf,e1)));


	double dat2 = at2_p - at2;

	double dat2_lin = arma::dot(SBGATPolyhedronGravityModelUQ::PartialAtan2PartialZf(Zf),dZf);

	assert(std::abs(dat2 - dat2_lin)/std::abs(dat2_lin) < tol);
	std::cout << "\t Passed TestPartialAtan2PartialZf ...\n";

}


void SBGATPolyhedronGravityModelUQ::TestPartialZfPartialUnitRf(double tol){

	std::cout << "\t In TestPartialZfPartialUnitRf ...\n";

	arma::vec::fixed<9> Rf = {0,1,2,3,4,5,6,7,8};
	Rf.subvec(0,2) = arma::normalise(Rf.subvec(0,2));
	Rf.subvec(3,5) = arma::normalise(Rf.subvec(3,5));
	Rf.subvec(6,8) = arma::normalise(Rf.subvec(6,8));



	arma::vec::fixed<3> r0_hat = Rf.subvec(0,2);
	arma::vec::fixed<3> r1_hat = Rf.subvec(3,5);
	arma::vec::fixed<3> r2_hat = Rf.subvec(6,8);

	double alpha = 1 + arma::dot(r0_hat,r1_hat) + arma::dot(r0_hat,r2_hat) + arma::dot(r1_hat,r2_hat);
	double gamma = arma::dot(r0_hat,arma::cross(r1_hat,r2_hat));

	arma::vec::fixed<2> Zf = {alpha,gamma};


	arma::vec::fixed<9> dRf = {1e-3,1e-3,-1e-3,2e-3,-1e-3,2e-3,1e-3,1e-3,-1e-3};

	arma::vec::fixed<9> Rf_p = Rf + dRf;

	arma::vec::fixed<3> r0_hat_p = Rf_p.subvec(0,2);
	arma::vec::fixed<3> r1_hat_p = Rf_p.subvec(3,5);
	arma::vec::fixed<3> r2_hat_p = Rf_p.subvec(6,8);

	double alpha_p = 1 + arma::dot(r0_hat_p,r1_hat_p) + arma::dot(r0_hat_p,r2_hat_p) + arma::dot(r1_hat_p,r2_hat_p);
	double gamma_p = arma::dot(r0_hat_p,arma::cross(r1_hat_p,r2_hat_p));

	arma::vec::fixed<2> Zf_p = {alpha_p,gamma_p};


	arma::vec::fixed<2> dZf = Zf_p - Zf;
	arma::vec::fixed<2> dZf_lin = SBGATPolyhedronGravityModelUQ::PartialZfPartialUnitRf(Rf) * dRf;

	assert(arma::norm(dZf - dZf_lin)/arma::norm(dZf_lin) < tol);
	std::cout << "\t Passed TestPartialZfPartialUnitRf ...\n";

}


