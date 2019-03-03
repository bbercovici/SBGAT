#include <SBGATPolyhedronGravityModelUQ.hpp>
#include <SBGATMassProperties.hpp>
#include <RigidBodyKinematics.hpp>
#include <vtkOBJReader.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>


#include <json.hpp>

#pragma omp declare reduction( + : arma::rowvec : omp_out += omp_in ) \
initializer( omp_priv = arma::zeros<arma::rowvec>(omp_orig.n_cols))


#pragma omp declare reduction( + : arma::mat : omp_out += omp_in ) \
initializer( omp_priv = arma::zeros<arma::mat>(omp_orig.n_rows,omp_orig.n_cols))




double SBGATPolyhedronGravityModelUQ::GetVariancePotential(double const * point) const{

	arma::rowvec partial = this -> GetPartialUPartialC(point);

	return arma::dot(partial, this -> P_CC * partial.t()); 

}

double SBGATPolyhedronGravityModelUQ::GetVariancePotential(const arma::vec::fixed<3> & point) const{

	return this -> GetVariancePotential(point.colptr(0));

}


arma::mat::fixed<3,3> SBGATPolyhedronGravityModelUQ::GetCovarianceAcceleration(double const * point) const{

	arma::mat partial = this -> GetPartialAPartialC(point);
	return partial *  this -> P_CC * partial.t();

}

arma::mat::fixed<3,3> SBGATPolyhedronGravityModelUQ::GetCovarianceAcceleration(const arma::vec::fixed<3> & point) const{
	return this -> GetCovarianceAcceleration(point.colptr(0));
}


void SBGATPolyhedronGravityModelUQ::GetVariancePotentialAccelerationCovariance(const arma::vec::fixed<3> & point,double & potential_var, 
	arma::mat::fixed<3,3> & acc_cov) const{
	throw(std::runtime_error("SBGATPolyhedronGravityModelUQ::GetVariancePotentialAccelerationCovariance is not implemented yet"));

}



arma::rowvec::fixed<10> SBGATPolyhedronGravityModelUQ::PartialUePartialXe(const arma::vec::fixed<3> & pos,const int & e) const{


	arma::rowvec::fixed<10> partial;

	SBGATPolyhedronGravityModel * pgm_model = SBGATPolyhedronGravityModel::SafeDownCast(this -> model);

	double Le = pgm_model -> GetLe(pos,e);
	arma::vec::fixed<3> r_ei_0 = pgm_model -> GetRe(pos,e); 
	arma::mat::fixed<3,6> R_ei_0 = {
		{r_ei_0[0],0,0,r_ei_0[1],r_ei_0[2],0},
		{0,r_ei_0[1],0,r_ei_0[0],0,r_ei_0[2]},
		{0,0,r_ei_0[2],0,r_ei_0[0],r_ei_0[1]}
	};

	arma::vec::fixed<3> Ee_times_r_ei_0 = R_ei_0 * pgm_model -> GetEeParam(e);

	partial(0) = arma::dot(r_ei_0,Ee_times_r_ei_0);
	partial.subvec(1,3) = 2 * Le * Ee_times_r_ei_0.t();
	partial.subvec(4,9) = Le * r_ei_0.t() * R_ei_0;

	return partial;

}


arma::rowvec::fixed<10> SBGATPolyhedronGravityModelUQ::PartialUfPartialXf(const arma::vec::fixed<3> & pos,
	const int & f) const{


	arma::rowvec::fixed<10> partial;

	SBGATPolyhedronGravityModel * pgm_model = SBGATPolyhedronGravityModel::SafeDownCast(this -> model);


	double omega_f = pgm_model -> GetOmegaf(pos,f);
	arma::vec::fixed<3> r_fi_0 = pgm_model -> GetRf(pos,f); 
	arma::mat::fixed<3,6> R_fi_0 = {
		{r_fi_0[0],0,0,r_fi_0[1],r_fi_0[2],0},
		{0,r_fi_0[1],0,r_fi_0[0],0,r_fi_0[2]},
		{0,0,r_fi_0[2],0,r_fi_0[0],r_fi_0[1]}
	};

	arma::vec::fixed<3> F_times_r_fi_0 = R_fi_0 * pgm_model -> GetFfParam(f);


	partial(0) = arma::dot(r_fi_0,F_times_r_fi_0);
	partial.subvec(1,3) = 2 * omega_f * F_times_r_fi_0.t();
	partial.subvec(4,9) = omega_f * r_fi_0.t() * R_fi_0;

	return - partial;

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

	this -> model -> GetVerticesInFacet(f,r0,r1,r2);

	arma::vec::fixed<3> R0 = arma::vec({r0[0],r0[1],r0[2]}) * this -> model -> GetScaleFactor();
	arma::vec::fixed<3> R1 = arma::vec({r1[0],r1[1],r1[2]}) * this -> model -> GetScaleFactor();
	arma::vec::fixed<3> R2 = arma::vec({r2[0],r2[1],r2[2]}) * this -> model -> GetScaleFactor();

	R0 -= pos;
	R1 -= pos;
	R2 -= pos;

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

	arma::vec::fixed<3> Nf = this -> model -> GetNonNormalizedFacetNormal(f);

	

	return SBGATPolyhedronGravityModelUQ::PartialFfPartialnf(arma::normalise(Nf)) * SBGATPolyhedronGravityModelUQ::PartialNormalizedVPartialNonNormalizedV(Nf) * this -> PartialNfPartialTf(f);

}



arma::mat::fixed<3,3> SBGATPolyhedronGravityModelUQ::PartialNormalizedVPartialNonNormalizedV(const arma::vec::fixed<3> & non_normalized_V){

	return (arma::eye<arma::mat>(3,3) / arma::norm(non_normalized_V) - non_normalized_V * non_normalized_V.t() / std::pow(arma::norm(non_normalized_V),3));

}



arma::mat::fixed<6,3> SBGATPolyhedronGravityModelUQ::PartialFfPartialnf(const arma::vec::fixed<3> & nf){

	arma::vec::fixed<3> e0 = {1,0,0};
	arma::vec::fixed<3> e1 = {0,1,0};
	arma::vec::fixed<3> e2 = {0,0,1};

	arma::mat::fixed<6,3> partial;

	partial.row(0) = 2 * nf.t() * e0 * e0.t();
	partial.row(1) = 2 * nf.t() * e1 * e1.t();
	partial.row(2) = 2 * nf.t() * e2 * e2.t();
	partial.row(3) = nf.t() *( e0 * e1.t() + e1 * e0.t());
	partial.row(4) = nf.t() *( e0 * e2.t() + e2 * e0.t());
	partial.row(5) = nf.t() *( e1 * e2.t() + e2 * e1.t());


	return partial;


}


arma::rowvec::fixed<6> SBGATPolyhedronGravityModelUQ::PartialLePartialAe(const arma::vec::fixed<3> & pos,const int & e) const{

	double r0[3],r1[3];

	this -> model -> GetVerticesOnEdge(e,r0,r1);

	arma::vec::fixed<3> r0_arma = arma::vec({r0[0],r0[1],r0[2]}) * this -> model -> GetScaleFactor() - pos;
	arma::vec::fixed<3> r1_arma = arma::vec({r1[0],r1[1],r1[2]}) * this -> model -> GetScaleFactor() - pos;


	double re_0 = arma::norm(r0_arma);
	double re_1 = arma::norm(r1_arma);
	double le = arma::norm(r1_arma - r0_arma);
	double beta_plus = re_0 + re_1 + le;
	double beta_minus = re_0 + re_1 - le;

	arma::vec::fixed<3> beta_vectors = {
		1./beta_plus + 1./beta_minus,
		1./beta_plus - 1./beta_minus,
		1./beta_plus - 1./beta_minus
	};

	arma::mat::fixed<3,6> partial_re_0_partial_Ae = arma::zeros<arma::mat>(3,6);
	arma::mat::fixed<3,6> partial_re_1_partial_Ae = arma::zeros<arma::mat>(3,6);

	partial_re_0_partial_Ae.cols(0,2) = arma::eye<arma::mat>(3,3);
	partial_re_1_partial_Ae.cols(3,5) = arma::eye<arma::mat>(3,3);


	arma::mat::fixed<3,6> partial_mat;
	partial_mat.row(0) = this -> PartialEdgeLengthPartialAe(e);
	partial_mat.row(1) = r0_arma.t() / re_0 * partial_re_0_partial_Ae;
	partial_mat.row(2) = r1_arma.t() / re_1 * partial_re_1_partial_Ae;


	return beta_vectors.t() * partial_mat ;

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
	this -> model -> GetVerticesOnEdge(e,r0,r1);

	arma::vec::fixed<3> r0_arma = {r0[0],r0[1],r0[2]};
	arma::vec::fixed<3> r1_arma = {r1[0],r1[1],r1[2]};

	arma::mat::fixed<3,6> mat;
	mat.cols(0,2) = arma::eye<arma::mat>(3,3);
	mat.cols(3,5) = - arma::eye<arma::mat>(3,3);

	return arma::normalise(r0_arma - r1_arma).t() * mat;


}

arma::rowvec::fixed<24> SBGATPolyhedronGravityModelUQ::PartialEqrPartialBe(const int & e,const int & q,const int & r) const{

	arma::vec::fixed<3> e_q = arma::zeros<arma::vec>(3);
	arma::vec::fixed<3> e_r = arma::zeros<arma::vec>(3);
	e_q(q) = 1;
	e_r(r) = 1;

	double r0[3],r1[3];

	this -> model -> GetVerticesOnEdge(e,r0,r1);

	arma::vec::fixed<3> r0_arma = this -> model -> GetScaleFactor() * arma::vec({r0[0],r0[1],r0[2]});
	arma::vec::fixed<3> r1_arma = this -> model -> GetScaleFactor() * arma::vec({r1[0],r1[1],r1[2]});

	int f0,f1;

	this -> model -> GetIndicesOfAdjacentFacets(e,f0,f1);

	arma::vec::fixed<3> N_e_f0 = this -> model -> GetNonNormalizedFacetNormal(f0) ;
	arma::vec::fixed<3> N_e_f1 = this -> model -> GetNonNormalizedFacetNormal(f1) ;

	arma::vec::fixed<3> n_e_f0 = arma::normalise(N_e_f0);
	arma::vec::fixed<3> n_e_f1 = arma::normalise(N_e_f1);

	arma::vec::fixed<3> Vr = RBK::tilde(r1_arma - r0_arma) * e_r;

	double le = arma::norm(r0_arma - r1_arma);

	double E_qr = arma::dot(e_q,(n_e_f1 * n_e_f1.t() - n_e_f0 * n_e_f0.t()) * Vr) / le;

	arma::mat::fixed<3,6> M = arma::zeros<arma::mat>(3,6);
	M.cols(0,2) = - arma::eye<arma::mat>(3,3);
	M.cols(3,5) = arma::eye<arma::mat>(3,3);

	arma::vec::fixed<13> first_partial;

	first_partial(0) = - E_qr;
	first_partial.subvec(1,6) = M.t() * RBK::tilde(e_r) * (arma::dot(n_e_f1,e_q) * n_e_f1 - arma::dot(n_e_f0,e_q) * n_e_f0);
	first_partial.subvec(7,9) = - (e_q * n_e_f0.t() + arma::dot(n_e_f0,e_q) * arma::eye<arma::mat>(3,3)) * Vr;
	first_partial.subvec(10,12) = (e_q * n_e_f1.t() + arma::dot(n_e_f1,e_q) * arma::eye<arma::mat>(3,3)) * Vr;
	
	first_partial *= 1/le;
	
	arma::mat::fixed<13,24> second_partial = arma::zeros<arma::mat>(13,24);
	second_partial.submat(0,0,0,5) = this -> PartialEdgeLengthPartialAe(e);	
	second_partial.submat(1,0,6,5) = arma::eye<arma::mat>(6,6);

	second_partial.submat(7,6,9,14) = SBGATPolyhedronGravityModelUQ::PartialNormalizedVPartialNonNormalizedV(N_e_f0) * this -> PartialNfPartialTf(f0);
	
	second_partial.submat(10,15,12,23) = SBGATPolyhedronGravityModelUQ::PartialNormalizedVPartialNonNormalizedV(N_e_f1) * this -> PartialNfPartialTf(f1);

	return first_partial.t() * second_partial;

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

void SBGATPolyhedronGravityModelUQ::TestPartials(std::string input,double tol,bool shape_in_meters){

	std::cout << "\tRunning SBGATPolyhedronGravityModelUQ::TestPartials on " << input << std::endl;
	
	SBGATPolyhedronGravityModelUQ::TestPartialNfPartialTf(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestPartialNormalizedVPartialNonNormalizedV(input,tol,shape_in_meters);	
	SBGATPolyhedronGravityModelUQ::TestPartialAtan2PartialZf(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestPartialZfPartialUnitRf(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestPartialOmegafPartialTf(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestPartialFfPartialnf(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestPartialFfPartialNonNormalizedNf(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestPartialFfPartialTf(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestPartialUfPartialXf(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestPartialUfPartialTf(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestPartialEdgeLengthPartialAe(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestPartialEPartialBe(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestPartialLePartialAe(input,tol,shape_in_meters);

	SBGATPolyhedronGravityModelUQ::TestPartialXePartialBe(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestPartialXfPartialTf(input,tol,shape_in_meters);

	SBGATPolyhedronGravityModelUQ::TestPartialBePartialC(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestPartialUePartialXe(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestPartialUePartialBe(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestPartialUfPartialC(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestPartialUePartialC(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestAddPartialSumUfPartialC(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestAddPartialSumUePartialC(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestAddPartialSumAccfPartialC(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestAddPartialSumAccePartialC(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestPartialUPartialC(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestGetPartialAPartialC(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestPartialOmegaPartialwC(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestPartialBodyFixedAccelerationfPartialC(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestPartialBodyFixedAccelerationfPartialwC(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestPartialSlopeArgumentPartialOmegaC(input,tol,shape_in_meters);
	SBGATPolyhedronGravityModelUQ::TestGetPartialSlopePartialwPartialC(input,tol,shape_in_meters);

}

void SBGATPolyhedronGravityModelUQ::TestPartialUfPartialTf(std::string filename,double tol,bool shape_in_meters){


	std::cout << "\t In TestPartialUfPartialTf ... ";
	int successes = 0;
	arma::arma_rng::set_seed(0);
	int N = 100;
	#pragma omp parallel for reduction(+:successes)
	
	for (int i = 0; i < N ; ++i){

		
		// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner =
		vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();


	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);

		int N_facets = pgm_filter -> GetN_facets();
		arma::ivec f_vec = arma::randi<arma::ivec>(1,arma::distr_param(0,N_facets - 1));
		int f = f_vec(0);
		arma::vec pos = {2,3,4};

	// Nominal Uf
		double Uf = pgm_filter -> GetUf(pgm_filter -> GetXf(pos,f));

	// Deviation
		arma::vec::fixed<9> delta_Tf = 1e-2 * arma::randn<arma::vec>(9)/ pgm_filter -> GetScaleFactor();

	// Linear dUf
		double dUf_lin = arma::dot(shape_uq.PartialUfPartialXf(pos,f) * shape_uq.PartialXfPartialTf(pos,f),  pgm_filter -> GetScaleFactor() * delta_Tf);

	// Apply Tf deviation
		shape_uq. ApplyTfDeviation(delta_Tf,f);

	// Perturbed Uf
		double Uf_p = pgm_filter -> GetUf(pgm_filter -> GetXf(pos,f));
		

	// Non-linear dUf
		double dUf = Uf_p - Uf;

		if(std::abs(dUf - dUf_lin) / std::abs(dUf_lin) < tol){
			++successes;
		}



	}

	std::cout << "\t Passed TestPartialUfPartialTf with " << double(successes)/N * 100 << " \% of sucesses. \n";


}
void SBGATPolyhedronGravityModelUQ::TestPartialUePartialBe(std::string filename,double tol,bool shape_in_meters){



	std::cout << "\t In TestPartialUePartialBe ... ";
	
	int successes = 0;
	arma::arma_rng::set_seed(0);
	int N = 100;

	
	// Reading
	vtkSmartPointer<vtkOBJReader> r = vtkSmartPointer<vtkOBJReader>::New();
	r -> SetFileName(filename.c_str());
	r -> Update(); 

	// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cl =
	vtkSmartPointer<vtkCleanPolyData>::New();
	cl -> SetInputConnection (r -> GetOutputPort());
	cl -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
	cl -> Update();	


	vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
	
	mass_prop -> SetInputConnection(cl -> GetOutputPort());
	if(shape_in_meters){
		mass_prop -> SetScaleMeters();
	}
	else{
		mass_prop -> SetScaleKiloMeters();
	}

	mass_prop -> Update();
	double xmin, xmax, ymin, ymax, zmin, zmax;
	mass_prop -> GetBoundingBox( xmin, xmax, ymin, ymax, zmin, zmax);

	arma::vec::fixed<3> pos = 1.5 * arma::vec({xmax,ymax,zmax});


	#pragma omp parallel for reduction(+:successes)
	
	for(int i = 0; i < N; ++i){

		// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner =
		vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 

		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}

		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);

		int N_edges = pgm_filter -> GetN_edges();
		arma::ivec e_vec = arma::randi<arma::ivec>(1,arma::distr_param(0,N_edges - 1));
		int e = e_vec(0);



		arma::vec::fixed<6> Ee_param = pgm_filter -> GetEeParam(e);

		arma::vec deviation = 1e-2 * arma::randn<arma::vec>(3 * pgm_filter -> GetN_vertices()) / pgm_filter -> GetScaleFactor();

		arma::rowvec::fixed<24> partial = shape_uq.PartialUePartialXe(pos,e) * shape_uq.PartialXePartialBe(pos,e);

		double Ue = pgm_filter -> GetUe( pgm_filter -> GetXe(pos,e));

		// Apply global deviation and get all dBes deviation
		arma::vec dBe = shape_uq.ApplyAndGetBeDeviation(deviation);

		double dUe = pgm_filter -> GetUe( pgm_filter -> GetXe(pos,e)) - Ue;
		double dUe_lin = arma::dot(partial, dBe.subvec(24 * e, 24 * e + 23) * pgm_filter -> GetScaleFactor());
		
		if(std::abs(dUe - dUe_lin)/std::abs(dUe_lin) < tol){
			++successes;
		}


	}
	std::cout << "\t Passed TestPartialUePartialBe with " << double(successes)/N * 100 << " \% of successes .\n";






}


void SBGATPolyhedronGravityModelUQ::TestPartialUePartialXe(std::string filename,double tol,bool shape_in_meters){

	std::cout << "\t In TestPartialUePartialXe ... ";

	int successes = 0;
	arma::arma_rng::set_seed(0);
	int N = 100;

	


	#pragma omp parallel for reduction(+:successes)

	for (int i = 0; i < N; ++i){

			// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner =
		vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);


  	// Test point 
		arma::vec::fixed<3> pos;
		if (shape_in_meters){
			pos = {1,3,4};
		}
		else{
			pos = {1e3,3e3,4e3};

		}

  	// Pick e
		int e = 1;

  	// Edge point
		arma::vec::fixed<3> rei_0 = pgm_filter -> GetRe(pos,e);

  	// Xe before
		arma::vec::fixed<10> Xe = pgm_filter -> GetXe(pos,e);

  	// Ue before
		double Ue = pgm_filter -> GetUe(Xe);

  	// Perturbation to Xe 
		arma::vec::fixed<10> dXe = 1e-3 * arma::randn<arma::vec>(10);

  	// Xe after
		arma::vec::fixed<10> Xe_p = Xe + dXe;

  	// Ue after
		double Ue_p = pgm_filter -> GetUe(Xe_p);

  	// Non-linear difference
		double dUe = Ue_p - Ue;

  	// Partial 

		arma::rowvec::fixed<10> partial_Ue_partial_Xe = shape_uq.PartialUePartialXe(pos,e);

  	// Linear difference
		double dUe_lin = arma::dot(partial_Ue_partial_Xe.t(), dXe);

		if(std::abs(dUe - dUe_lin) / std::abs(dUe_lin) < tol){
			++successes;
		}
	}

	std::cout << "\t Passed TestPartialUePartialXe with " << double(successes) / N * 100 <<  " \% of successes .\n";


}
void SBGATPolyhedronGravityModelUQ::TestPartialUfPartialXf(std::string filename,double tol,bool shape_in_meters){


	std::cout << "\t In TestPartialUfPartialXf ... ";

	int successes = 0;
	arma::arma_rng::set_seed(0);
	#pragma omp parallel for reduction(+:successes)

	for (int i =0; i < 100; ++i){
		
	// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(reader -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);

	// Test point 
		arma::vec::fixed<3> pos;
		if (shape_in_meters){
			pos = {1,3,4};
		}
		else{
			pos = {1e3,3e3,4e3};

		}

  	// Pick f
		int N_facets = pgm_filter -> GetN_facets();

		arma::ivec f_vec = arma::randi<arma::ivec>(1,arma::distr_param(0,N_facets - 1));
		int f = f_vec(0);
  	// Edge point
		arma::vec::fixed<3> rfi_0 = pgm_filter -> GetRf(pos,f);

  	// Xe before
		arma::vec::fixed<10> Xf = pgm_filter -> GetXf(pos,f);

  	// Ue before
		double Uf = pgm_filter -> GetUf(Xf);

  	// Perturbation to Xe 
		arma::vec::fixed<10> dXf = 1e-2 * arma::randn<arma::vec>(10);

  	// Xe after
		arma::vec::fixed<10> Xf_p = Xf + dXf;

  	// Ue after
		double Uf_p = pgm_filter -> GetUf(Xf_p);

  	// Non-linear difference
		double dUf = Uf_p - Uf;

  	// Partial 

		arma::rowvec::fixed<10> partial_Uf_partial_Xf = shape_uq. PartialUfPartialXf(pos,f);

  	// Linear difference
		double dUf_lin = arma::dot(partial_Uf_partial_Xf.t(), dXf);

		if(std::abs(dUf - dUf_lin) / std::abs(dUf_lin) < tol){
			++successes;
		}
	}


	std::cout << "\t Passed TestPartialUfPartialXf with " << successes<< "\% of successes. \n";


}

void SBGATPolyhedronGravityModelUQ::TestPartialXfPartialTf(std::string filename,double tol,bool shape_in_meters){


	std::cout << "\t In TestPartialXfPartialTf ... ";
	int successes = 0;
	arma::arma_rng::set_seed(0);
	


	
	// Reading
	vtkSmartPointer<vtkOBJReader> r = vtkSmartPointer<vtkOBJReader>::New();
	r -> SetFileName(filename.c_str());
	r -> Update(); 

	// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cl =
	vtkSmartPointer<vtkCleanPolyData>::New();
	cl -> SetInputConnection (r -> GetOutputPort());
	cl -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
	cl -> Update();	


	vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
	
	mass_prop -> SetInputConnection(cl -> GetOutputPort());
	if(shape_in_meters){
		mass_prop -> SetScaleMeters();
	}
	else{
		mass_prop -> SetScaleKiloMeters();
	}

	mass_prop -> Update();
	double xmin, xmax, ymin, ymax, zmin, zmax;
	mass_prop -> GetBoundingBox( xmin, xmax, ymin, ymax, zmin, zmax);

	arma::vec::fixed<3> pos = 1.5 * arma::vec({xmax,ymax,zmax});



	#pragma omp parallel for reduction(+:successes)
	for (int i = 0; i < 100 ; ++i){

		
	// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner =
		vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();



	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);
		int N_facets = pgm_filter -> GetN_facets();

		arma::ivec f_vec = arma::randi<arma::ivec>(1,arma::distr_param(0,N_facets - 1));
		int f = f_vec(0);		

	// Nominal Xf
		arma::vec::fixed<10> Xf = pgm_filter -> GetXf(pos,f);

	// Deviation
		arma::vec::fixed<9> delta_Tf = 1e-2 * arma::randn<arma::vec>(9)/ pgm_filter -> GetScaleFactor();

	// Linear dXf
		arma::vec::fixed<10> dXf_lin = shape_uq.PartialXfPartialTf(pos,f) *  pgm_filter -> GetScaleFactor() * delta_Tf;

	// Apply Tf deviation
		shape_uq. ApplyTfDeviation(delta_Tf,f);

	// Perturbed Xf
		arma::vec::fixed<10> Xf_p = pgm_filter -> GetXf(pos,f);

	// Non-linear dXf
		arma::vec::fixed<10> dXf = Xf_p - Xf;

		if(arma::norm(dXf - dXf_lin) / arma::norm(dXf_lin) < tol){
			++successes;
		}



	}

	std::cout << "\t Passed TestPartialXfPartialTf with " << successes << " \% of sucesses. \n";

}
void SBGATPolyhedronGravityModelUQ::TestPartialOmegafPartialTf(std::string filename,double tol,bool shape_in_meters){


	std::cout << "\t In TestPartialOmegafPartialTf ... ";



	
	// Reading
	vtkSmartPointer<vtkOBJReader> r = vtkSmartPointer<vtkOBJReader>::New();
	r -> SetFileName(filename.c_str());
	r -> Update(); 

	// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cl =
	vtkSmartPointer<vtkCleanPolyData>::New();
	cl -> SetInputConnection (r -> GetOutputPort());
	cl -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
	cl -> Update();	


	vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
	
	mass_prop -> SetInputConnection(cl -> GetOutputPort());
	if(shape_in_meters){
		mass_prop -> SetScaleMeters();
	}
	else{
		mass_prop -> SetScaleKiloMeters();
	}

	mass_prop -> Update();
	double xmin, xmax, ymin, ymax, zmin, zmax;
	mass_prop -> GetBoundingBox( xmin, xmax, ymin, ymax, zmin, zmax);

	arma::vec::fixed<3> pos = 1.5 * arma::vec({xmax,ymax,zmax});



	int successes = 0;
	arma::arma_rng::set_seed(0);
	#pragma omp parallel for reduction(+:successes)

	for (int i = 0; i < 100; ++i){
		
	// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner =
		vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

		vtkSmartPointer<vtkPolyData> polydata = cleaner -> GetOutput();
		
	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputData(polydata);
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);
		int N_facets = pgm_filter -> GetN_facets();

		arma::ivec f_vec = arma::randi<arma::ivec>(1,arma::distr_param(0,N_facets - 1));
		int f = f_vec(0);		

	// Nominal omega_f
		double omega_f = pgm_filter -> GetOmegaf(pos,f);

	// Deviation
		arma::vec::fixed<9> delta_Tf = 1e-2 * arma::randn<arma::vec>(9)/ pgm_filter -> GetScaleFactor();

	// Linear dXf
		double domega_f_lin = arma::dot(shape_uq.PartialOmegafPartialTf(pos,f), pgm_filter -> GetScaleFactor() * delta_Tf);

	// Apply Tf deviation
		shape_uq. ApplyTfDeviation(delta_Tf,f);

	// Perturbed omega_f
		double omega_f_p = pgm_filter -> GetOmegaf(pos,f);


	// Non-linear domega_f
		double domega_f = omega_f_p - omega_f;

		
		if(std::abs(domega_f - domega_f_lin) / std::abs(domega_f_lin) < tol){
			++successes;
		}

	}
	std::cout << "\t Passed TestPartialOmegafPartialTf with " << successes  << " \% of successes. \n";





}
void SBGATPolyhedronGravityModelUQ::TestPartialFfPartialTf(std::string filename,double tol,bool shape_in_meters){


	std::cout << "\t In TestPartialFfPartialTf ... ";
	int successes = 0;
	arma::arma_rng::set_seed(0);
	#pragma omp parallel for reduction(+:successes)

	for (int i = 0; i < 100; ++i){
		
	// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner =
		vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);
		int N_facets = pgm_filter -> GetN_facets();

		arma::ivec f_vec = arma::randi<arma::ivec>(1,arma::distr_param(0,N_facets - 1));
		int f = f_vec(0);

	// Nominal dyad
		arma::vec::fixed<6> Ff = pgm_filter -> GetFfParam(f);

	// Deviation
		arma::vec::fixed<9> delta_Tf = 1e-3 * arma::randn<arma::vec>(9)/ pgm_filter -> GetScaleFactor();

	// Linear difference
		arma::vec::fixed<6> dFf_lin = shape_uq.PartialFfPartialTf(f) * pgm_filter -> GetScaleFactor() * delta_Tf;

	// Apply Tf deviation
		shape_uq. ApplyTfDeviation(delta_Tf,f);

	// Perturbed dyad
		arma::vec::fixed<6> Ff_p = pgm_filter -> GetFfParam(f);

	// Non-linear difference
		arma::vec::fixed<6> dFf = Ff_p - Ff;


		if(arma::norm(dFf_lin - dFf) / arma::norm(dFf_lin) < 1e-2){
			++successes;
		}
	}
	std::cout << "\t Passed TestPartialFfPartialTf with " << successes<<  " \% of successes. \n";




}




arma::mat SBGATPolyhedronGravityModelUQ::PartialOmegaPartialwC(const arma::vec::fixed<3> & Omega) const{
	

	int N_C = this -> model -> GetN_vertices();


	SBGATPolyhedronGravityModel * pgm_model = SBGATPolyhedronGravityModel::SafeDownCast(this -> model);

	arma::mat partial = arma::zeros<arma::mat>(3, 1 + 3 * N_C);
	
	arma::mat::fixed<3,3> PB = pgm_model -> GetPrincipalAxes();
	arma::vec::fixed<3> rotation_axis_principal_frame = PB * arma::normalise(Omega);

	partial.col(0) = arma::normalise(Omega);
	partial.cols(1,3 * N_C) = - 4 * arma::norm(Omega) * PB.t() * RBK::tilde(rotation_axis_principal_frame) * this -> GetPartialSigmaPartialC();

	return partial;

}


arma::mat SBGATPolyhedronGravityModelUQ::PartialBodyFixedAccelerationfPartialOmegaC(const int & f,const arma::vec::fixed<3> & Omega) const{

	int N_C = this -> model -> GetN_vertices();

	arma::mat partial = arma::zeros<arma::mat>(3, 3 + 3 * N_C );

	partial.cols(0,2) = this -> PartialBodyFixedAccelerationfPartialOmega(f,Omega);
	partial.cols(3,3 * N_C + 3 - 1 ) = this -> PartialBodyFixedAccelerationfPartialC(f,Omega);

	return partial;

}

arma::mat SBGATPolyhedronGravityModelUQ::PartialBodyFixedAccelerationfPartialC(const int & f,const arma::vec::fixed<3> & Omega) const{

	SBGATPolyhedronGravityModel * pgm_model = SBGATPolyhedronGravityModel::SafeDownCast(this -> model);

	arma::mat::fixed<3,9> mat;
	mat.cols(0,2) = arma::eye<arma::mat>(3,3);
	mat.cols(3,5) = arma::eye<arma::mat>(3,3);
	mat.cols(6,8) = arma::eye<arma::mat>(3,3);

	arma::vec::fixed<3> facet_center = this -> model  -> GetFacetCenter(f);

	return (this -> GetPartialAPartialC(facet_center)
		+ pgm_model -> GetGravityGradient(facet_center) * 1./3 * mat * this -> PartialTfPartialC(f)
		+ RBK::tilde(Omega) * RBK::tilde(Omega) * (this -> GetPartialComPartialC() - 
			1./3 * mat * this -> PartialTfPartialC(f)));

}

double SBGATPolyhedronGravityModelUQ::PartialSlopePartialSlopeArgument(const double & u){

	return 1./std::sqrt(1 - u * u);
}


arma::rowvec SBGATPolyhedronGravityModelUQ::GetPartialSlopePartialwPartialC(const int & f,const arma::vec::fixed<3> & Omega) const{

	SBGATPolyhedronGravityModel * pgm_model = SBGATPolyhedronGravityModel::SafeDownCast(this -> model);

	arma::vec::fixed<3> body_fixed_acc = pgm_model -> GetBodyFixedAccelerationf(f,Omega);

	double slope = std::acos(- arma::dot(arma::normalise(body_fixed_acc),arma::normalise(pgm_model -> GetNonNormalizedFacetNormal(f))));

	double u = - std::cos(slope);


	return (this -> PartialSlopePartialSlopeArgument(u) 
		* this -> PartialSlopeArgumentPartialOmegaC(f,Omega,body_fixed_acc) 
		* this -> PartialOmegaCPartialwC(Omega));


}


arma::sp_mat SBGATPolyhedronGravityModelUQ::PartialOmegaCPartialwC(const arma::vec::fixed<3> & Omega) const{

	int N_C = this -> model -> GetN_vertices();

	arma::sp_mat partial(3 + 3 * N_C, 1 + 3 * N_C);


	partial.rows(0,2) = this -> PartialOmegaPartialwC(Omega);


	partial.submat(3,1,3 + 3 * N_C - 1, 3 * N_C) = arma::eye<arma::mat>(3 * N_C,3 * N_C);


	return partial;

}

arma::rowvec SBGATPolyhedronGravityModelUQ::PartialSlopeArgumentPartialOmegaC(const int & f,const arma::vec::fixed<3> & Omega, const arma::vec::fixed<3> & body_fixed_acc) const{

	int N_C = this -> model -> GetN_vertices();

	arma::vec::fixed<6> unit_vectors;
	arma::vec::fixed<3> Nf = this -> model -> GetNonNormalizedFacetNormal(f);
	unit_vectors.subvec(0,2) = arma::normalise(Nf);
	unit_vectors.subvec(3,5) = arma::normalise(body_fixed_acc);

	arma::mat partial_mat = arma::zeros<arma::mat>(6, 3 * N_C + 3);

	partial_mat.rows(0,2) = this -> PartialNormalizedVPartialNonNormalizedV(body_fixed_acc) * this -> PartialBodyFixedAccelerationfPartialOmegaC(f,Omega);
	partial_mat.submat(3,3,5,3 * N_C + 2) = this -> PartialNormalizedVPartialNonNormalizedV(Nf) * this -> PartialNfPartialTf(f) * this -> PartialTfPartialC(f);


	return (unit_vectors.t() * partial_mat);

}





arma::mat::fixed<3,3> SBGATPolyhedronGravityModelUQ::PartialBodyFixedAccelerationfPartialOmega(const int & f,const arma::vec::fixed<3> & Omega) const{


	arma::vec::fixed<3> Pf = this -> model  -> GetFacetCenter(f);
	const arma::vec::fixed<3> G = SBGATMassProperties::SafeDownCast(this -> model)  -> GetCenterOfMass();

	return RBK::tilde(arma::cross(Omega,Pf - G)) + RBK::tilde(Omega) * RBK::tilde(Pf - G);


}




void SBGATPolyhedronGravityModelUQ::TestPartialNormalizedVPartialNonNormalizedV(std::string filename,double tol,bool shape_in_meters){

	std::cout << "\t In TestPartialNormalizedVPartialNonNormalizedV ... ";

	int N = 100;
	int successes = 0;
	arma::arma_rng::set_seed(0);
	#pragma omp parallel for reduction(+:successes)

	for (int i = 0; i < N ;++i){
  	// non-normalized vector
		arma::vec::fixed<3> X = arma::randn<arma::vec>(3);

  	// normalized vector
		arma::vec::fixed<3> x = arma::normalise(X);


  	// disturbance 
		arma::vec::fixed<3> dX = 1e-2 * arma::randn<arma::vec>(3);

  	// normalized vector, perturbed
		arma::vec::fixed<3> x_p = arma::normalise(X + dX);

  	// difference, non-linear
		arma::vec::fixed<3> dx = x_p - x;

  	// difference, linear
		arma::vec::fixed<3> dx_lin = SBGATPolyhedronGravityModelUQ::PartialNormalizedVPartialNonNormalizedV(X) * dX;


		if(arma::norm(dx - dx_lin)/arma::norm(dx_lin) < tol){
			++successes;
		}

	}
	std::cout << "\t Passed TestPartialNormalizedVPartialNonNormalizedV with " << double (successes) /N * 100<< " \% of successes.\n";

}
void SBGATPolyhedronGravityModelUQ::TestPartialNfPartialTf(std::string filename,double tol,bool shape_in_meters){

	std::cout << "\t In TestPartialNfPartialTf ... ";

	
	int successes = 0;
	arma::arma_rng::set_seed(0);
	#pragma omp parallel for reduction(+:successes)
	for (int i = 0; i < 100; ++i){
	// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner =
		vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();


		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);

		int N_facets = pgm_filter -> GetN_facets();

		arma::ivec f_vec = arma::randi<arma::ivec>(1,arma::distr_param(0,N_facets - 1));
		int f = f_vec(0);		

		arma::vec::fixed<9> delta_Tf = 1e-2 * arma::randn<arma::vec>(9) / pgm_filter -> GetScaleFactor();

		arma::vec::fixed<3> N = pgm_filter -> GetNonNormalizedFacetNormal(f);

	// Linear difference
		arma::vec::fixed<3> dN_lin = shape_uq . PartialNfPartialTf(f) * delta_Tf * pgm_filter -> GetScaleFactor();

	// Apply Tf deviation
		shape_uq. ApplyTfDeviation(delta_Tf,f);

	// Perturbed normal
		arma::vec::fixed<3> N_p = pgm_filter -> GetNonNormalizedFacetNormal(f);
		
	// Non-linear difference
		arma::vec::fixed<3> dN = N_p - N;

		if(arma::norm(dN - dN_lin)/arma::norm(dN_lin) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestPartialNfPartialTf with " << successes <<" \% of successes. \n";



}
void SBGATPolyhedronGravityModelUQ::TestPartialFfPartialnf(std::string filename,double tol,bool shape_in_meters){

	std::cout << "\t In TestPartialFfPartialnf ... ";

	int N = 100;
	int successes = 0;
	arma::arma_rng::set_seed(0);

	#pragma omp parallel for reduction(+:successes)
	for (int i = 0; i < N; ++i){

		arma::vec::fixed<3> nf = {1,2,-3};
		nf = arma::normalise(nf);


	// Dyad
		arma::mat::fixed<3,3> Ff = nf * nf.t();

	// Dyad parametrization
		arma::vec::fixed<6> Ff_vec = {Ff(0,0),Ff(1,1),Ff(2,2),Ff(0,1),Ff(0,2),Ff(1,2)};

	// Perturbation

		arma::vec::fixed<3> dnf = 1e-2 * arma::randn<arma::vec>(3);

		arma::vec::fixed<3> nf_p = dnf + nf;


	// Pertyrbed dyad
		arma::mat::fixed<3,3> Ff_p = nf_p * nf_p.t();

	// Perturbed dyad parametrization
		arma::vec::fixed<6> Ff_p_vec = {Ff_p(0,0),Ff_p(1,1),Ff_p(2,2),Ff_p(0,1),Ff_p(0,2),Ff_p(1,2)};


	// Differences
		arma::vec::fixed<6> dFf = Ff_p_vec - Ff_vec;
		arma::vec::fixed<6> dFf_linear = SBGATPolyhedronGravityModelUQ::PartialFfPartialnf(nf) * dnf;

		if(arma::norm(dFf_linear - dFf ) / arma::norm(dFf) < tol){
			++successes;
		}
	}
	std::cout << "\t Passed TestPartialFfPartialnf with " << double(successes)/N * 100<< " \% of successes.\n";



}

void SBGATPolyhedronGravityModelUQ::TestPartialFfPartialNonNormalizedNf(std::string filename,double tol,bool shape_in_meters){

	std::cout << "\t In TestPartialFfPartialNonNormalizedNf ... ";

	int N = 100;
	int successes = 0;
	arma::arma_rng::set_seed(0);

	#pragma omp parallel for reduction(+:successes)
	for (int i =0; i < N; ++i){
		arma::vec::fixed<3> Nf = {1,2,3};
		arma::vec::fixed<3> nf = arma::normalise(Nf);


	// Dyad
		arma::mat::fixed<3,3> Ff = nf * nf.t();

	// Dyad parametrization
		arma::vec::fixed<6> Ff_vec = {Ff(0,0),Ff(1,1),Ff(2,2),Ff(0,1),Ff(0,2),Ff(1,2)};

	// Perturbation

		arma::vec::fixed<3> dNf = 1e-2 * arma::randn<arma::vec>(3);

		arma::vec::fixed<3> Nf_p = dNf + Nf;
		arma::vec::fixed<3> nf_p = arma::normalise(Nf_p);

	// Pertyrbed dyad
		arma::mat::fixed<3,3> Ff_p = nf_p * nf_p.t();

	// Perturbed dyad parametrization
		arma::vec::fixed<6> Ff_p_vec = {Ff_p(0,0),Ff_p(1,1),Ff_p(2,2),Ff_p(0,1),Ff_p(0,2),Ff_p(1,2)};

	// Differences
		arma::vec::fixed<6> dFf = Ff_p_vec - Ff_vec;
		arma::vec::fixed<6> dFf_linear = SBGATPolyhedronGravityModelUQ::PartialFfPartialnf(nf) * SBGATPolyhedronGravityModelUQ::PartialNormalizedVPartialNonNormalizedV(Nf) * dNf;


		if(arma::norm(dFf_linear - dFf ) / arma::norm(dFf) < tol){
			++successes;
		}
	}


	std::cout << "\t Passed TestPartialFfPartialNonNormalizedNf with " <<double(successes) / N * 100 << " \% of successes. \n";

}







void SBGATPolyhedronGravityModelUQ::TestPartialLePartialAe(std::string filename,double tol,bool shape_in_meters){


	std::cout << "\t In TestPartialLePartialAe ... ";
	
	int successes = 0;
	arma::arma_rng::set_seed(0);
	arma::vec pos = {1,2,3};
	#pragma omp parallel for reduction(+:successes)

	for (int i = 0; i < 100; ++i){
	// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner =
		vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);

		int e = 1;
		arma::vec::fixed<6> delta_Ae = 1e-2 * arma::randn<arma::vec>(6) / pgm_filter -> GetScaleFactor();

		double Le = pgm_filter -> GetLe(pos,e);
		

	// Apply Ae deviation
		shape_uq.ApplyAeDeviation(delta_Ae,e);

	// Perturbed length
		double Le_p = pgm_filter -> GetLe(pos,e);

	// Non-linear difference
		double dLe = Le_p - Le;

	// Linear difference
		double dLe_lin = arma::dot(shape_uq.PartialLePartialAe(pos,e), pgm_filter -> GetScaleFactor() * delta_Ae);


		if(std::abs(dLe - dLe_lin)/std::abs(dLe_lin) < tol){
			++successes;
		}
	}
	std::cout << "\t Passed TestPartialLePartialAe with " << successes<< " \% of successes. \n";


}


void SBGATPolyhedronGravityModelUQ::TestPartialEdgeLengthPartialAe(std::string filename,double tol,bool shape_in_meters){

	std::cout << "\t In TestPartialEdgeLengthPartialAe ... ";
	
	int successes = 0;
	arma::arma_rng::set_seed(0);
	#pragma omp parallel for reduction(+:successes)
	
	for (int i = 0; i < 100; ++i){
	// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner =
		vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);

		int e = 1;
		arma::vec::fixed<6> delta_Ae = 1e-2 * arma::randn<arma::vec>(6) / pgm_filter -> GetScaleFactor() ;

		
	// Nominal length
		double le = pgm_filter -> GetEdgeLength(e);

	// Apply Ae deviation
		shape_uq.ApplyAeDeviation(delta_Ae,e);


	// Perturbed length
		double le_p = pgm_filter -> GetEdgeLength(e);


	// Non-linear difference
		double dle = le_p - le;

	// Linear difference
		double dle_lin = arma::dot(shape_uq.PartialEdgeLengthPartialAe(e), pgm_filter -> GetScaleFactor() * delta_Ae);


		if(std::abs(dle - dle_lin)/std::abs(dle_lin) < tol){
			++successes;
		}
	}
	std::cout << "\t Passed TestPartialEdgeLengthPartialAe with " << successes<< " \% of successes. \n";

}


void SBGATPolyhedronGravityModelUQ::TestPartialEPartialBe(std::string filename,double tol,bool shape_in_meters){


	std::cout << "\t In TestPartialEPartialBe ... ";


	
	int successes = 0;
	arma::arma_rng::set_seed(0);
	int N = 100;
	#pragma omp parallel for reduction(+:successes)
	for(int i = 0; i < N; ++i){

		// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner =
		vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);

		int N_edges = pgm_filter -> GetN_edges();
		arma::ivec e_vec = arma::randi<arma::ivec>(1,arma::distr_param(0,N_edges - 1));
		int e = e_vec(0);

		arma::vec::fixed<6> Ee_param = pgm_filter -> GetEeParam(e);

		arma::vec deviation = 1e-2 * arma::randn<arma::vec>(3 * pgm_filter -> GetN_vertices()) / pgm_filter -> GetScaleFactor();

		auto partial = shape_uq.PartialEPartialBe(e);

		// Apply global deviation and get all dBes deviation
		arma::vec dBe = shape_uq.ApplyAndGetBeDeviation(deviation);

		arma::vec::fixed<6> Ee_param_p = pgm_filter -> GetEeParam(e);


		arma::vec::fixed<6> dEe_param = Ee_param_p - Ee_param;
		arma::vec::fixed<6> dEe_param_lin = partial * dBe.subvec(24 * e,24 * e + 23) * pgm_filter -> GetScaleFactor();


		if(arma::norm(dEe_param - dEe_param_lin)/arma::norm(dEe_param_lin) < tol){
			++successes;
		}



	}
	std::cout << "\t Passed TestPartialEPartialBe with " << double(successes)/N * 100 << " \% of successes .\n";



}




arma::sp_mat  SBGATPolyhedronGravityModelUQ::PartialBePartialC(const int & e) const{

	arma::sp_mat table(24, 3 * vtkPolyData::SafeDownCast(this -> model -> GetInput()) -> GetNumberOfPoints());

	// Ae
	int v0_e,v1_e;
	this -> model -> GetIndicesVerticesOnEdge(e,v0_e,v1_e);
	table.submat(0,3 * v0_e, 2,3 * v0_e + 2) = arma::eye<arma::mat>(3,3);
	table.submat(3,3 * v1_e, 5,3 * v1_e + 2) = arma::eye<arma::mat>(3,3);

	// Ts
	int f0_e,f1_e;
	this -> model -> GetIndicesOfAdjacentFacets(e,f0_e,f1_e);

	// T0
	int v0_f0_e,v1_f0_e,v2_f0_e;
	this -> model -> GetIndicesVerticesInFacet(f0_e, v0_f0_e,v1_f0_e,v2_f0_e);

	table.submat(6,3 * v0_f0_e, 8,3 * v0_f0_e + 2) = arma::eye<arma::mat>(3,3);
	table.submat(9,3 * v1_f0_e, 11,3 * v1_f0_e + 2) = arma::eye<arma::mat>(3,3);
	table.submat(12,3 * v2_f0_e, 14,3 * v2_f0_e + 2) = arma::eye<arma::mat>(3,3);

	// T1
	int v0_f1_e,v1_f1_e,v2_f1_e;
	this -> model -> GetIndicesVerticesInFacet(f1_e, v0_f1_e,v1_f1_e,v2_f1_e);

	table.submat(15,3 * v0_f1_e, 17,3 * v0_f1_e + 2) = arma::eye<arma::mat>(3,3);
	table.submat(18,3 * v1_f1_e, 20,3 * v1_f1_e + 2) = arma::eye<arma::mat>(3,3);
	table.submat(21,3 * v2_f1_e, 23,3 * v2_f1_e + 2) = arma::eye<arma::mat>(3,3);

	return table;

}



void SBGATPolyhedronGravityModelUQ::ApplyAeDeviation(arma::vec::fixed<6> delta_Ae,const int & e){
	
	int v0_index,v1_index;

	this -> model -> GetIndicesVerticesOnEdge(e,v0_index,v1_index);

	double r0[3],r1[3];
	vtkPolyData * polydata = vtkPolyData::SafeDownCast(this -> model -> GetInput());

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


	this -> model -> Modified();

	SBGATPolyhedronGravityModel::SafeDownCast(this -> model) -> Update();


}


void SBGATPolyhedronGravityModelUQ::ApplyTfDeviation(arma::vec::fixed<9> delta_Tf,const int & f){
	
	int v0_index,v1_index,v2_index;
	this -> model -> GetIndicesVerticesInFacet(f,v0_index,v1_index,v2_index);

	double r0[3],r1[3],r2[3];
	vtkPolyData * polydata = vtkPolyData::SafeDownCast(this -> model -> GetInput());

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
	this -> model -> Modified();

	SBGATPolyhedronGravityModel::SafeDownCast(this -> model) -> Update();
	

}



arma::vec SBGATPolyhedronGravityModelUQ::GetBe() const{
	
	vtkPolyData * polydata = vtkPolyData::SafeDownCast(this -> model -> GetInput());

	int N = polydata -> GetNumberOfPoints();

	double r[3];
	arma::vec points(3 * N);

	for (int i = 0; i < N; ++i){
		polydata -> GetPoint(i,r);
		points(3 * i) = r[0];
		points(3 * i + 1) = r[1];
		points(3 * i + 2) = r[2];

	}

	
	int N_edges = polydata -> GetNumberOfPoints() + polydata -> GetNumberOfCells() - 2;

	arma::vec be(24 * N_edges);
	for (int e = 0 ; e < N_edges; ++e){
		be.subvec(24 * e,24 * e + 23) = this -> PartialBePartialC(e) * points;
	}

	return be;

}



arma::vec SBGATPolyhedronGravityModelUQ::ApplyAndGetBeDeviation(const arma::vec & delta){
	
	vtkPolyData * polydata = vtkPolyData::SafeDownCast(this -> model -> GetInput());
	assert(int(delta.n_rows / 3) == polydata -> GetNumberOfPoints());

	int N = polydata -> GetNumberOfPoints();

	double r[3];

	for (int i = 0; i < N; ++i){
		polydata -> GetPoint(i,r);


		r[0] += delta(3 * i);
		r[1] += delta(3 * i + 1);
		r[2] += delta(3 * i + 2);

		polydata -> GetPoints() -> SetPoint(i,r);


	}

	polydata -> GetPoints() -> Modified();

	polydata -> Modified();
	this -> model -> Modified();

	this -> model -> Update();

	int N_edges = polydata -> GetNumberOfPoints() + polydata -> GetNumberOfCells() - 2;

	arma::vec be_deviation(24 * N_edges);
	for (int e = 0 ; e < N_edges; ++e){
		be_deviation.subvec(24 * e,24 * e + 23) = this -> PartialBePartialC(e) * delta;
	}

	return be_deviation;

}


void SBGATPolyhedronGravityModelUQ::TestPartialAtan2PartialZf(std::string filename,double tol,bool shape_in_meters){

	std::cout << "\t In TestPartialAtan2PartialZf ... ";

	int N = 100;
	int successes = 0;
	arma::arma_rng::set_seed(0);
	#pragma omp parallel for reduction(+:successes)

	for (int i =0; i < N; ++i){
		arma::vec::fixed<2> e1 = {1,0};
		arma::vec::fixed<2> e2 = {0,1};

		arma::vec::fixed<2> Zf = {0.1,-0.2};

		double at2 = 2 * std::atan(arma::dot(Zf,e2) / (arma::norm(Zf) + arma::dot(Zf,e1)));


		arma::vec::fixed<2> dZf = 1e-3 * arma::randn<arma::vec>(2);

		double at2_p = 2 * std::atan(arma::dot(Zf + dZf,e2) / (arma::norm(Zf + dZf) + arma::dot(Zf + dZf,e1)));


		double dat2 = at2_p - at2;

		double dat2_lin = arma::dot(SBGATPolyhedronGravityModelUQ::PartialAtan2PartialZf(Zf),dZf);

		if(std::abs(dat2 - dat2_lin)/std::abs(dat2_lin) < tol){
			++successes;
		}
	}
	std::cout << "\t Passed TestPartialAtan2PartialZf with " << double (successes)/N * 100<<  " \% of successes. \n";

}


void SBGATPolyhedronGravityModelUQ::TestPartialZfPartialUnitRf(std::string filename,double tol,bool shape_in_meters){

	std::cout << "\t In TestPartialZfPartialUnitRf ...";
	int N = 100;
	int successes = 0;
	arma::arma_rng::set_seed(0);
	#pragma omp parallel for reduction(+:successes)

	for (int i = 0; i < N ; ++i){
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

		arma::vec::fixed<9> dRf = 1e-2 * arma::randn<arma::vec>(9);

		arma::vec::fixed<9> Rf_p = Rf + dRf;

		arma::vec::fixed<3> r0_hat_p = Rf_p.subvec(0,2);
		arma::vec::fixed<3> r1_hat_p = Rf_p.subvec(3,5);
		arma::vec::fixed<3> r2_hat_p = Rf_p.subvec(6,8);

		double alpha_p = 1 + arma::dot(r0_hat_p,r1_hat_p) + arma::dot(r0_hat_p,r2_hat_p) + arma::dot(r1_hat_p,r2_hat_p);
		double gamma_p = arma::dot(r0_hat_p,arma::cross(r1_hat_p,r2_hat_p));

		arma::vec::fixed<2> Zf_p = {alpha_p,gamma_p};


		arma::vec::fixed<2> dZf = Zf_p - Zf;
		arma::vec::fixed<2> dZf_lin = SBGATPolyhedronGravityModelUQ::PartialZfPartialUnitRf(Rf) * dRf;

		if(arma::norm(dZf - dZf_lin)/arma::norm(dZf_lin) < tol){
			++successes;
		}
	}
	std::cout << "\t Passed TestPartialZfPartialUnitRf with " << double(successes) / N * 100 << " \% of successes. \n";

}



arma::rowvec SBGATPolyhedronGravityModelUQ::GetPartialUPartialC(const arma::vec::fixed<3> & pos) const{

	int N_C =this -> model -> GetN_vertices();

	arma::rowvec partial = arma::zeros<arma::rowvec>(3 * N_C);

	this -> AddPartialSumUePartialC(pos,partial);
	this -> AddPartialSumUfPartialC(pos,partial);

	return 0.5 * arma::datum::G * this -> model -> GetDensity() * partial;

}


void SBGATPolyhedronGravityModelUQ::AddPartialSumUePartialC(const arma::vec::fixed<3> & pos,arma::rowvec & partial) const{
	int N_f = vtkPolyData::SafeDownCast(this -> model -> GetInput()) -> GetNumberOfCells();

	int N_C =this -> model -> GetN_vertices();

	int N_e = N_C + N_f - 2;

	#pragma omp parallel for reduction(+:partial)
	for (int e = 0; e < N_e; ++e){
		partial += this -> PartialUePartialXe(pos,e) * this -> PartialXePartialBe(pos,e) * this -> PartialBePartialC(e);
	}

}

void SBGATPolyhedronGravityModelUQ::AddPartialSumUfPartialC(const arma::vec::fixed<3> & pos,arma::rowvec & partial) const{
	int N_f = vtkPolyData::SafeDownCast(this -> model -> GetInput()) -> GetNumberOfCells();

	#pragma omp parallel for reduction(+:partial)
	for (int f = 0; f < N_f; ++f){
		partial += this -> PartialUfPartialXf(pos,f) * this -> PartialXfPartialTf(pos,f) * this -> PartialTfPartialC(f);
	}


}



void SBGATPolyhedronGravityModelUQ::AddPartialSumAccePartialC(const arma::vec::fixed<3> & pos,arma::mat & partial) const{
	int N_f = vtkPolyData::SafeDownCast(this -> model -> GetInput()) -> GetNumberOfCells();

	int N_C =this -> model -> GetN_vertices();

	int N_e = N_C + N_f - 2;

	#pragma omp parallel for reduction(+:partial)
	for (int e = 0; e < N_e; ++e){
		partial += this -> PartialAccePartialXe(pos,e) * this -> PartialXePartialBe(pos,e) * this -> PartialBePartialC(e);
	}

}

void SBGATPolyhedronGravityModelUQ::AddPartialSumAccfPartialC(const arma::vec::fixed<3> & pos,arma::mat & partial) const{
	int N_f = vtkPolyData::SafeDownCast(this -> model -> GetInput()) -> GetNumberOfCells();

	#pragma omp parallel for reduction(+:partial)
	for (int f = 0; f < N_f; ++f){
		partial += this -> PartialAccfPartialXf(pos,f) * this -> PartialXfPartialTf(pos,f) * this -> PartialTfPartialC(f);
	}

}







arma::mat SBGATPolyhedronGravityModelUQ::GetPartialAPartialC(const arma::vec::fixed<3> & pos) const{


	int N_C =this -> model -> GetN_vertices();

	arma::mat partial = arma::zeros<arma::mat>(3,3 * N_C);

	this -> AddPartialSumAccePartialC(pos,partial);
	this -> AddPartialSumAccfPartialC(pos,partial);

	return arma::datum::G * this -> model -> GetDensity() * partial;

}

arma::mat::fixed<3,10> SBGATPolyhedronGravityModelUQ::PartialAccePartialXe(const arma::vec::fixed<3> & pos,const int & e) const{

	arma::mat::fixed<3,10> partial;

	SBGATPolyhedronGravityModel * pgm_model = SBGATPolyhedronGravityModel::SafeDownCast(this -> model);


	double Le = pgm_model -> GetLe(pos,e);
	arma::vec::fixed<3> r_ei_0 = pgm_model -> GetRe(pos,e); 
	arma::mat::fixed<3,6> R_ei_0 = {
		{r_ei_0[0],0,0,r_ei_0[1],r_ei_0[2],0},
		{0,r_ei_0[1],0,r_ei_0[0],0,r_ei_0[2]},
		{0,0,r_ei_0[2],0,r_ei_0[0],r_ei_0[1]}
	};

	arma::vec::fixed<6> Ee = pgm_model -> GetEeParam(e);
	
	arma::mat::fixed<3,3> Ee_mat = {
		{Ee(0),Ee(3),Ee(4)},
		{Ee(3),Ee(1),Ee(5)},
		{Ee(4),Ee(5),Ee(2)}
	};



	partial.col(0) = Ee_mat * r_ei_0;
	partial.cols(1,3) = Le * Ee_mat;
	partial.cols(4,9) = Le * R_ei_0;

	return - partial;
}

arma::mat::fixed<3,10> SBGATPolyhedronGravityModelUQ::PartialAccfPartialXf(const arma::vec::fixed<3> & pos,const int & f) const{

	arma::mat::fixed<3,10> partial;

	SBGATPolyhedronGravityModel * pgm_model = SBGATPolyhedronGravityModel::SafeDownCast(this -> model);


	double omega_f = pgm_model -> GetOmegaf(pos,f);
	arma::vec::fixed<3> r_fi_0 = pgm_model -> GetRf(pos,f); 
	arma::mat::fixed<3,6> R_fi_0 = {
		{r_fi_0[0],0,0,r_fi_0[1],r_fi_0[2],0},
		{0,r_fi_0[1],0,r_fi_0[0],0,r_fi_0[2]},
		{0,0,r_fi_0[2],0,r_fi_0[0],r_fi_0[1]}
	};

	arma::vec::fixed<6> Ff = pgm_model -> GetFfParam(f);

	arma::mat::fixed<3,3> Ff_mat = {
		{Ff(0),Ff(3),Ff(4)},
		{Ff(3),Ff(1),Ff(5)},
		{Ff(4),Ff(5),Ff(2)}
	};

	partial.col(0) = Ff_mat * r_fi_0;
	partial.cols(1,3) = omega_f * Ff_mat;
	partial.cols(4,9) = omega_f * R_fi_0;

	return partial;

}


void SBGATPolyhedronGravityModelUQ::TestAddPartialSumUePartialC(std::string filename,double tol,bool shape_in_meters){

	std::cout << "\t In TestAddPartialSumUePartialC ...";

	// MC
	int N = 100;
	int successes = 0;
	arma::arma_rng::set_seed(0);
	
	// Reading
	vtkSmartPointer<vtkOBJReader> r = vtkSmartPointer<vtkOBJReader>::New();
	r -> SetFileName(filename.c_str());
	r -> Update(); 

	// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cl =
	vtkSmartPointer<vtkCleanPolyData>::New();
	cl -> SetInputConnection (r -> GetOutputPort());
	cl -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
	cl -> Update();	


	vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
	
	mass_prop -> SetInputConnection(cl -> GetOutputPort());
	if(shape_in_meters){
		mass_prop -> SetScaleMeters();
	}
	else{
		mass_prop -> SetScaleKiloMeters();
	}

	mass_prop -> Update();
	double xmin, xmax, ymin, ymax, zmin, zmax;
	mass_prop -> GetBoundingBox( xmin, xmax, ymin, ymax, zmin, zmax);

	arma::vec::fixed<3> pos = 1.5 * arma::vec({xmax,ymax,zmax});


	#pragma omp parallel for reduction(+:successes)
	
	for (int i = 0; i < N ; ++i){
			// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);

		double SumUe = 0;

		int N_edges = pgm_filter -> GetN_edges();

		for (int e = 0; e < N_edges; ++e){
			SumUe += pgm_filter -> GetUe(pgm_filter -> GetXe(pos,e));
		}

		arma::rowvec partial = arma::zeros<arma::rowvec>(3 * pgm_filter -> GetN_vertices());
		shape_uq.AddPartialSumUePartialC(pos,partial);
		arma::vec deviation = 1e-2 * arma::randn<arma::vec>(3 * pgm_filter -> GetN_vertices()) / pgm_filter -> GetScaleFactor();

		shape_uq.ApplyDeviation(deviation);

		double SumUe_p = 0;

		for (int e = 0; e < N_edges; ++e){
			SumUe_p += pgm_filter -> GetUe(pgm_filter -> GetXe(pos,e));
		}

		double dSumUe = SumUe_p - SumUe;
		double dSumUe_lin = arma::dot(partial,pgm_filter -> GetScaleFactor() * deviation);

		if(std::abs(dSumUe - dSumUe_lin)/std::abs(dSumUe_lin) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestAddPartialSumUePartialC with " << double(successes) / N * 100 << " \% of successes.\n";

}


void SBGATPolyhedronGravityModelUQ::TestAddPartialSumAccePartialC(std::string filename,double tol,bool shape_in_meters){

	std::cout << "\t In TestAddPartialSumAccePartialC ...";

	// MC
	int N = 100;
	int successes = 0;
	arma::arma_rng::set_seed(0);
	
	// Reading
	vtkSmartPointer<vtkOBJReader> r = vtkSmartPointer<vtkOBJReader>::New();
	r -> SetFileName(filename.c_str());
	r -> Update(); 

	// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cl =
	vtkSmartPointer<vtkCleanPolyData>::New();
	cl -> SetInputConnection (r -> GetOutputPort());
	cl -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
	cl -> Update();	


	vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
	
	mass_prop -> SetInputConnection(cl -> GetOutputPort());
	if(shape_in_meters){
		mass_prop -> SetScaleMeters();
	}
	else{
		mass_prop -> SetScaleKiloMeters();
	}

	mass_prop -> Update();
	double xmin, xmax, ymin, ymax, zmin, zmax;
	mass_prop -> GetBoundingBox( xmin, xmax, ymin, ymax, zmin, zmax);

	arma::vec::fixed<3> pos = 1.5 * arma::vec({xmax,ymax,zmax});


	#pragma omp parallel for reduction(+:successes)
	
	for (int i = 0; i < N ; ++i){
			// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);

		arma::vec::fixed<3> SumAcce = arma::zeros<arma::vec>(3);

		int N_edges = pgm_filter -> GetN_edges();

		for (int e = 0; e < N_edges; ++e){
			SumAcce += pgm_filter -> GetAe(pgm_filter -> GetXe(pos,e));
		}

		arma::mat partial = arma::zeros<arma::mat>(3,3 * pgm_filter -> GetN_vertices());
		shape_uq.AddPartialSumAccePartialC(pos,partial);
		arma::vec deviation = 1e-2 * arma::randn<arma::vec>(3 * pgm_filter -> GetN_vertices()) / pgm_filter -> GetScaleFactor();

		shape_uq.ApplyDeviation(deviation);

		arma::vec::fixed<3> SumAcce_p = arma::zeros<arma::vec>(3);

		for (int e = 0; e < N_edges; ++e){
			SumAcce_p += pgm_filter -> GetAe(pgm_filter -> GetXe(pos,e));
		}

		arma::vec::fixed<3> dSumAcce = SumAcce_p - SumAcce;
		arma::vec::fixed<3> dSumAcce_lin = partial * deviation * pgm_filter -> GetScaleFactor();

		if(arma::norm(dSumAcce - dSumAcce_lin)/arma::norm(dSumAcce_lin) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestAddPartialSumAccePartialC with " << double(successes) / N * 100 << " \% of successes.\n";

}

void SBGATPolyhedronGravityModelUQ::TestAddPartialSumAccfPartialC(std::string filename,double tol,bool shape_in_meters){

	std::cout << "\t In TestAddPartialSumAccfPartialC ...";

	// MC
	int N = 100;
	int successes = 0;
	arma::arma_rng::set_seed(0);
	
	// Reading
	vtkSmartPointer<vtkOBJReader> r = vtkSmartPointer<vtkOBJReader>::New();
	r -> SetFileName(filename.c_str());
	r -> Update(); 

	// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cl =
	vtkSmartPointer<vtkCleanPolyData>::New();
	cl -> SetInputConnection (r -> GetOutputPort());
	cl -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
	cl -> Update();	


	vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
	
	mass_prop -> SetInputConnection(cl -> GetOutputPort());
	if(shape_in_meters){
		mass_prop -> SetScaleMeters();
	}
	else{
		mass_prop -> SetScaleKiloMeters();
	}

	mass_prop -> Update();
	double xmin, xmax, ymin, ymax, zmin, zmax;
	mass_prop -> GetBoundingBox( xmin, xmax, ymin, ymax, zmin, zmax);

	arma::vec::fixed<3> pos = 1.5 * arma::vec({xmax,ymax,zmax});


	#pragma omp parallel for reduction(+:successes)
	
	for (int i = 0; i < N ; ++i){
			// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);

		arma::vec::fixed<3> SumAccf = arma::zeros<arma::vec>(3);

		int N_facets = vtkPolyData::SafeDownCast(pgm_filter -> GetInput()) -> GetNumberOfCells() ;

		for (int f = 0; f < N_facets; ++f){
			SumAccf += pgm_filter -> GetAf(pgm_filter -> GetXf(pos,f));
		}

		arma::mat partial = arma::zeros<arma::mat>(3,3 * pgm_filter -> GetN_vertices());
		shape_uq.AddPartialSumAccfPartialC(pos,partial);
		arma::vec deviation = 1e-2 * arma::randn<arma::vec>(3 * pgm_filter -> GetN_vertices()) / pgm_filter -> GetScaleFactor();

		shape_uq.ApplyDeviation(deviation);

		arma::vec::fixed<3> SumAccf_p = arma::zeros<arma::vec>(3);

		for (int f = 0; f < N_facets; ++f){
			SumAccf_p += pgm_filter -> GetAf(pgm_filter -> GetXf(pos,f));
		}

		arma::vec::fixed<3> dSumAccf = SumAccf_p - SumAccf;
		arma::vec::fixed<3> dSumAccf_lin = partial *  pgm_filter -> GetScaleFactor() * deviation;

		if(arma::norm(dSumAccf - dSumAccf_lin)/arma::norm(dSumAccf_lin) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestAddPartialSumAccfPartialC with " << double(successes) / N * 100 << " \% of successes.\n";

}




void SBGATPolyhedronGravityModelUQ::TestAddPartialSumUfPartialC(std::string filename,double tol,bool shape_in_meters){

	std::cout << "\t In TestAddPartialSumUfPartialC ...";

	// MC
	int N = 100;
	int successes = 0;
	arma::arma_rng::set_seed(0);
	
	// Reading
	vtkSmartPointer<vtkOBJReader> r = vtkSmartPointer<vtkOBJReader>::New();
	r -> SetFileName(filename.c_str());
	r -> Update(); 

	// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cl =
	vtkSmartPointer<vtkCleanPolyData>::New();
	cl -> SetInputConnection (r -> GetOutputPort());
	cl -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
	cl -> Update();	


	vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
	
	mass_prop -> SetInputConnection(cl -> GetOutputPort());
	if(shape_in_meters){
		mass_prop -> SetScaleMeters();
	}
	else{
		mass_prop -> SetScaleKiloMeters();
	}

	mass_prop -> Update();
	double xmin, xmax, ymin, ymax, zmin, zmax;
	mass_prop -> GetBoundingBox( xmin, xmax, ymin, ymax, zmin, zmax);

	arma::vec::fixed<3> pos = 1.5 * arma::vec({xmax,ymax,zmax});


	#pragma omp parallel for reduction(+:successes)
	
	for (int i = 0; i < N ; ++i){
			// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);

		double SumUf = 0;
		int N_facets = pgm_filter -> GetN_facets();

		for (int f = 0; f < N_facets; ++f){
			SumUf += pgm_filter -> GetUf(pgm_filter -> GetXf(pos,f));
		}

		arma::rowvec partial = arma::zeros<arma::rowvec>(3 * pgm_filter -> GetN_vertices());
		shape_uq.AddPartialSumUfPartialC(pos,partial);
		
		arma::vec deviation = 1e-2 * arma::randn<arma::vec>(3 * pgm_filter -> GetN_vertices()) / pgm_filter -> GetScaleFactor();

		shape_uq.ApplyDeviation(deviation);

		double SumUf_p = 0;

		for (int f = 0; f < N_facets ; ++f){
			SumUf_p += pgm_filter -> GetUf(pgm_filter -> GetXf(pos,f));
		}

		double dSumUf = SumUf_p - SumUf;
		double dSumUf_lin = arma::dot(partial,pgm_filter -> GetScaleFactor() * deviation);

		if(std::abs(dSumUf - dSumUf_lin)/std::abs(dSumUf_lin) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestAddPartialSumUfPartialC with " << double(successes) / N * 100 << " \% of successes.\n";



}

void SBGATPolyhedronGravityModelUQ::TestPartialUPartialC(std::string filename,double tol,bool shape_in_meters){

	std::cout << "\t In TestPartialUPartialC ...";

	// MC
	int N = 100;
	int successes = 0;
	arma::arma_rng::set_seed(0);
	

	// Reading
	vtkSmartPointer<vtkOBJReader> r = vtkSmartPointer<vtkOBJReader>::New();
	r -> SetFileName(filename.c_str());
	r -> Update(); 

	// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cl =
	vtkSmartPointer<vtkCleanPolyData>::New();
	cl -> SetInputConnection (r -> GetOutputPort());
	cl -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
	cl -> Update();	


	vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
	
	mass_prop -> SetInputConnection(cl -> GetOutputPort());
	if(shape_in_meters){
		mass_prop -> SetScaleMeters();
	}
	else{
		mass_prop -> SetScaleKiloMeters();
	}

	mass_prop -> Update();
	double xmin, xmax, ymin, ymax, zmin, zmax;
	mass_prop -> GetBoundingBox( xmin, xmax, ymin, ymax, zmin, zmax);

	arma::vec::fixed<3> pos = 1.5 * arma::vec({xmax,ymax,zmax});

	#pragma omp parallel for reduction(+:successes)
	for (int i = 0; i < N ; ++i){

			// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);

		double U = pgm_filter -> GetPotential(pos);
		arma::rowvec dUdC = shape_uq.GetPartialUPartialC(pos);
		arma::vec deviation = 1e-2 * arma::randn<arma::vec>(3 * pgm_filter -> GetN_vertices()) / pgm_filter -> GetScaleFactor();

		shape_uq.ApplyDeviation(deviation);

		double U_p = pgm_filter -> GetPotential(pos);
		double dU = U_p - U;
		double dU_lin = arma::dot(dUdC, pgm_filter -> GetScaleFactor() * deviation);

		if(std::abs(dU - dU_lin)/std::abs(dU_lin) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestPartialUPartialC with " << double(successes) / N * 100 << " \% of successes.\n";


}


void SBGATPolyhedronGravityModelUQ::TestGetPartialAPartialC(std::string filename,double tol,bool shape_in_meters){

	std::cout << "\t In TestGetPartialAPartialC ...";

	// MC
	int N = 100;
	int successes = 0;
	arma::arma_rng::set_seed(0);

	arma::vec::fixed<3> pos = {100,200,300};
	
	#pragma omp parallel for reduction(+:successes)
	for (int i = 0; i < N ; ++i){

			// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		int N_C = pgm_filter -> GetN_vertices();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);

		arma::vec::fixed<3> A = pgm_filter -> GetAcceleration(pos);
		arma::mat dAdC = shape_uq.GetPartialAPartialC(pos);
		arma::vec deviation = 1e-1 * arma::randn<arma::vec>(N_C * 3) / pgm_filter -> GetScaleFactor();

		shape_uq.ApplyDeviation(deviation);

		arma::vec::fixed<3> A_p = pgm_filter -> GetAcceleration(pos);
		arma::vec::fixed<3> dA = A_p - A;
		arma::vec::fixed<3> dA_lin = dAdC * deviation * pgm_filter -> GetScaleFactor();

		if(arma::norm(dA - dA_lin)/arma::norm(dA_lin) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestGetPartialAPartialC with " << double(successes) / N * 100 << " \% of successes.\n";


}

void SBGATPolyhedronGravityModelUQ::TestPartialXePartialBe(std::string filename,double tol,bool shape_in_meters){

	std::cout << "\t In TestPartialXePartialBe ... ";


	
	// Reading
	vtkSmartPointer<vtkOBJReader> r = vtkSmartPointer<vtkOBJReader>::New();
	r -> SetFileName(filename.c_str());
	r -> Update(); 

	// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cl =
	vtkSmartPointer<vtkCleanPolyData>::New();
	cl -> SetInputConnection (r -> GetOutputPort());
	cl -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
	cl -> Update();	


	vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
	
	mass_prop -> SetInputConnection(cl -> GetOutputPort());
	if(shape_in_meters){
		mass_prop -> SetScaleMeters();
	}
	else{
		mass_prop -> SetScaleKiloMeters();
	}

	mass_prop -> Update();
	double xmin, xmax, ymin, ymax, zmin, zmax;
	mass_prop -> GetBoundingBox( xmin, xmax, ymin, ymax, zmin, zmax);

	arma::vec::fixed<3> pos = 1.5 * arma::vec({xmax,ymax,zmax});

	int successes = 0;
	arma::arma_rng::set_seed(0);
	int N = 100;

#pragma omp parallel for reduction(+:successes)
	for(int i = 0; i < N; ++i){
		// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner =
		vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);

		int N_edges = pgm_filter -> GetN_edges();
		arma::ivec e_vec = arma::randi<arma::ivec>(1,arma::distr_param(0,N_edges - 1));
		int e = e_vec(0);

		arma::vec deviation = 1e-2 * arma::randn<arma::vec>(3 * pgm_filter -> GetN_vertices()) / pgm_filter -> GetScaleFactor();

		auto partial = shape_uq.PartialXePartialBe(pos,e);
		arma::vec::fixed<10> Xe = pgm_filter -> GetXe(pos,e);

		// Apply global deviation and get all dBes deviation
		arma::vec dBe = shape_uq.ApplyAndGetBeDeviation(deviation);
		
		arma::vec::fixed<10> dXe_lin = partial * dBe.subvec(24 * e,24 * e + 23) *  pgm_filter -> GetScaleFactor();

		arma::vec::fixed<10> Xe_p = pgm_filter -> GetXe(pos,e);

		arma::vec::fixed<10> dXe = Xe_p - Xe;

		if(arma::norm(dXe - dXe_lin)/arma::norm(dXe_lin) < tol){
			++successes;
		}
	}
	std::cout << "\t Passed TestPartialXePartialBe with " << double(successes)/N * 100 << " \% of successes .\n";

}

void SBGATPolyhedronGravityModelUQ::TestPartialUePartialC(std::string filename,double tol,bool shape_in_meters){

	std::cout << "\t In TestPartialUePartialC ...";

	// MC

	int N = 100;
	int successes = 0;
	arma::arma_rng::set_seed(0);
	
	// Reading
	vtkSmartPointer<vtkOBJReader> r = vtkSmartPointer<vtkOBJReader>::New();
	r -> SetFileName(filename.c_str());
	r -> Update(); 

	// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cl =
	vtkSmartPointer<vtkCleanPolyData>::New();
	cl -> SetInputConnection (r -> GetOutputPort());
	cl -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
	cl -> Update();	


	vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
	
	mass_prop -> SetInputConnection(cl -> GetOutputPort());
	if(shape_in_meters){
		mass_prop -> SetScaleMeters();
	}
	else{
		mass_prop -> SetScaleKiloMeters();
	}

	mass_prop -> Update();
	double xmin, xmax, ymin, ymax, zmin, zmax;
	mass_prop -> GetBoundingBox( xmin, xmax, ymin, ymax, zmin, zmax);

	arma::vec::fixed<3> pos = 1.5 * arma::vec({xmax,ymax,zmax});


	#pragma omp parallel for reduction(+:successes)
	
	for (int i = 0; i < N ; ++i){

			// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);
		
		int N_edges = pgm_filter -> GetN_edges();
		arma::ivec e_vec = arma::randi<arma::ivec>(1,arma::distr_param(0,N_edges - 1));
		int e = e_vec(0);

		double Ue = pgm_filter -> GetUe(pgm_filter -> GetXe(pos,e));
		auto dUedC = shape_uq . PartialUePartialXe(pos,e) * shape_uq . PartialXePartialBe(pos,e) * shape_uq . PartialBePartialC(e);

		arma::vec deviation = 1e-2 * arma::randn<arma::vec>(3 * pgm_filter -> GetN_vertices()) / pgm_filter -> GetScaleFactor();

		shape_uq.ApplyDeviation(deviation);

		double Ue_p = pgm_filter -> GetUe(pgm_filter -> GetXe(pos,e));
		double dUe = Ue_p - Ue;
		double dUe_lin = arma::dot(dUedC,pgm_filter -> GetScaleFactor() * deviation);

		if(std::abs(dUe - dUe_lin)/std::abs(dUe_lin) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestPartialUePartialC with " << double(successes) / N * 100 << " \% of successes.\n";

}


void SBGATPolyhedronGravityModelUQ::TestPartialBePartialC(std::string filename,double tol,bool shape_in_meters){

	std::cout << "\t In TestPartialBePartialC ...";

	// MC

	int N = 100;
	int successes = 0;
	arma::arma_rng::set_seed(0);
	
	// Reading
	vtkSmartPointer<vtkOBJReader> r = vtkSmartPointer<vtkOBJReader>::New();
	r -> SetFileName(filename.c_str());
	r -> Update(); 

	// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cl =
	vtkSmartPointer<vtkCleanPolyData>::New();
	cl -> SetInputConnection (r -> GetOutputPort());
	cl -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
	cl -> Update();	


	vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
	
	mass_prop -> SetInputConnection(cl -> GetOutputPort());
	if(shape_in_meters){
		mass_prop -> SetScaleMeters();
	}
	else{
		mass_prop -> SetScaleKiloMeters();
	}

	mass_prop -> Update();
	double xmin, xmax, ymin, ymax, zmin, zmax;
	mass_prop -> GetBoundingBox( xmin, xmax, ymin, ymax, zmin, zmax);

	arma::vec::fixed<3> pos = 1.5 * arma::vec({xmax,ymax,zmax});


	#pragma omp parallel for reduction(+:successes)
	
	for (int i = 0; i < N ; ++i){

			// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);
		
		int N_edges = pgm_filter -> GetN_edges();
		arma::ivec e_vec = arma::randi<arma::ivec>(1,arma::distr_param(0,N_edges - 1));
		int e = e_vec(0);



		arma::vec Be = shape_uq.GetBe().subvec(24 * e, 24 * e + 23);
		arma::sp_mat dBedC = shape_uq . PartialBePartialC(e);

		arma::vec deviation = 1e-2 * arma::randn<arma::vec>(3 * pgm_filter -> GetN_vertices()) / pgm_filter -> GetScaleFactor();

		shape_uq.ApplyDeviation(deviation);

		arma::vec Be_p = shape_uq.GetBe().subvec(24 * e, 24 * e + 23);
		arma::vec dBe = Be_p - Be;
		arma::vec dBe_lin = dBedC * deviation;


		if(arma::norm(dBe - dBe_lin)/arma::norm(dBe_lin) < tol){
			++successes;
		}
		else{

			arma::mat devs = arma::zeros<arma::mat>(dBe.n_rows,2);
			devs.col(0) = dBe;
			devs.col(1) = dBe_lin;

			arma::mat Bes = arma::zeros<arma::mat>(dBe.n_rows,2);
			Bes.col(0) = Be_p;
			Bes.col(1) = Be;


			std::cout << "Edge: " << e << std::endl;
			std::cout << arma::norm(devs.col(0) - devs.col(1)) / arma::norm(devs.col(1)) * 100 << std::endl;
			std::cout << devs << std::endl;
			std::cout << Bes << std::endl << std::endl;

		}

	}
	std::string status;

	if (successes == N){
		status = "\t Passed ";
	}
	else{
		status = "\t FAILED ";
	}

	std::cout << status << "TestPartialBePartialC with " << double(successes) / N * 100 << " \% of successes.\n";

}


void SBGATPolyhedronGravityModelUQ::TestPartialUfPartialC(std::string filename,double tol,bool shape_in_meters){

	std::cout << "\t In TestPartialUfPartialC ...";

	// MC

	int N = 100;
	int successes = 0;
	arma::arma_rng::set_seed(0);
	
	// Reading
	vtkSmartPointer<vtkOBJReader> r = vtkSmartPointer<vtkOBJReader>::New();
	r -> SetFileName(filename.c_str());
	r -> Update(); 

	// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cl =
	vtkSmartPointer<vtkCleanPolyData>::New();
	cl -> SetInputConnection (r -> GetOutputPort());
	cl -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
	cl -> Update();	


	vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
	
	mass_prop -> SetInputConnection(cl -> GetOutputPort());
	if(shape_in_meters){
		mass_prop -> SetScaleMeters();
	}
	else{
		mass_prop -> SetScaleKiloMeters();
	}

	mass_prop -> Update();
	double xmin, xmax, ymin, ymax, zmin, zmax;
	mass_prop -> GetBoundingBox( xmin, xmax, ymin, ymax, zmin, zmax);

	arma::vec::fixed<3> pos = 1.5 * arma::vec({xmax,ymax,zmax});

	#pragma omp parallel for reduction(+:successes)
	
	
	for (int i = 0; i < N ; ++i){

			// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);
		
		int N_f = pgm_filter -> GetN_facets();
		arma::ivec f_vec = arma::randi(1,arma::distr_param(0,N_f - 1));
		int f = f_vec(0);

		double Uf = pgm_filter -> GetUf(pgm_filter -> GetXf(pos,f));
		auto dUfdC = shape_uq . PartialUfPartialXf(pos,f) * shape_uq . PartialXfPartialTf(pos,f) * shape_uq . PartialTfPartialC(f);

		arma::vec deviation = 1e-2 * arma::randn<arma::vec>(3 * pgm_filter -> GetN_vertices()) / pgm_filter -> GetScaleFactor();

		shape_uq.ApplyDeviation(deviation);

		double Uf_p = pgm_filter -> GetUf(pgm_filter -> GetXf(pos,f));
		double dUf = Uf_p - Uf;
		double dUf_lin = arma::dot(dUfdC,pgm_filter -> GetScaleFactor() * deviation);

		if(std::abs(dUf - dUf_lin)/std::abs(dUf_lin) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestPartialUfPartialC with " << double(successes) / N * 100 << " \% of successes.\n";

}



void SBGATPolyhedronGravityModelUQ::RunMCUQ(std::string path_to_shape,
	const double & density,
	const bool & shape_in_meters,
	const arma::mat & P_CC,
	const unsigned int & N_samples,
	const arma::vec::fixed<3> & position,
	std::vector<arma::vec> & deviations,
	std::vector<arma::vec::fixed<3> > & accelerations,
	std::vector<double> & potentials,
	std::vector<vtkSmartPointer<vtkPolyData> > & saved_shapes){

	std::vector< std::vector < arma::vec::fixed<3> > > all_accelerations(N_samples);
	std::vector< std::vector< double> > all_potentials(N_samples);
	std::vector< arma::vec::fixed<3> > all_positions = {position};
	
	accelerations = std::vector<arma::vec::fixed<3> >(N_samples);
	potentials = std::vector<double>(N_samples);
	deviations = std::vector<arma::vec>(N_samples);


	SBGATPolyhedronGravityModelUQ::RunMCUQ(path_to_shape,
		density,
		shape_in_meters,
		P_CC,
		N_samples,
		all_positions,
		deviations,
		all_accelerations,
		all_potentials,
		saved_shapes);

	for (int i = 0; i < N_samples; ++i){
		accelerations[i] = all_accelerations[i][0];
		potentials[i] = all_potentials[i][0];
	}

}


void SBGATPolyhedronGravityModelUQ::RunMCUQ(std::string path_to_shape,
	const double & density,
	const bool & shape_in_meters,
	const arma::mat & P_CC,
	const unsigned int & N_samples,
	const std::vector<arma::vec::fixed<3> > & all_positions,
	std::vector<arma::vec> & deviations,
	std::vector<std::vector<arma::vec::fixed<3> >> & all_accelerations,
	std::vector<std::vector<double > > & all_potentials,
	std::vector<vtkSmartPointer<vtkPolyData> > & saved_shapes){

	if (all_accelerations.size() == 0)
		all_accelerations = std::vector< std::vector < arma::vec::fixed<3> > >(N_samples);
	if (all_potentials.size() == 0)
		all_potentials = std::vector< std::vector < double > > (N_samples);
	if (deviations.size() == 0)
		deviations = std::vector<arma::vec>(N_samples);

	// Reading
	vtkSmartPointer<vtkOBJReader> reader_mc = vtkSmartPointer<vtkOBJReader>::New();
	reader_mc -> SetFileName(path_to_shape.c_str());
	reader_mc -> Update(); 

		// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cleaner_mc = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner_mc -> SetInputConnection (reader_mc -> GetOutputPort());
	cleaner_mc -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );

	cleaner_mc -> Update();


	#pragma omp parallel for
	for (int i = 0; i < N_samples ; ++i){

		vtkSmartPointer<vtkPolyData> shape_copy = vtkSmartPointer<vtkPolyData>::New();
		shape_copy -> DeepCopy(cleaner_mc -> GetOutput());

		// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter_mc = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter_mc -> SetInputData(shape_copy);
		pgm_filter_mc -> SetDensity(density); 
		if (shape_in_meters){
			pgm_filter_mc -> SetScaleMeters();
		}
		else{
			pgm_filter_mc -> SetScaleKiloMeters();

		}
		pgm_filter_mc -> Update();
		SBGATPolyhedronGravityModelUQ shape_uq_mc;
		
		shape_uq_mc.SetModel(pgm_filter_mc);
		shape_uq_mc.SetCovariance(P_CC);

		arma::mat C_CC = shape_uq_mc.GetCovarianceSquareRoot(false);

		deviations[i] = C_CC * arma::randn<arma::vec>(3 * pgm_filter_mc -> GetN_vertices());
		shape_uq_mc.ApplyDeviation(deviations[i]);

		for (auto pos : all_positions){

			arma::vec::fixed<3> acc;
			double pot;

			pgm_filter_mc -> GetPotentialAcceleration(pos,pot,acc);
			all_accelerations[i].push_back(acc);
			all_potentials[i].push_back(pot);

		}


		if (i < saved_shapes.size()){
			vtkSmartPointer<vtkPolyData> shape = vtkSmartPointer<vtkPolyData>::New();
			shape -> DeepCopy(pgm_filter_mc -> GetInput());
			saved_shapes[i] = shape;
		}
	}


}




void SBGATPolyhedronGravityModelUQ::TestGetPartialSlopePartialwPartialC(std::string filename,double tol,bool shape_in_meters){

	std::cout << "\t In TestGetPartialSlopePartialwPartialC ...";

	// MC
	int N = 100;
	int successes = 0;
	arma::arma_rng::set_seed(0);

	#pragma omp parallel for reduction(+:successes)
	for (int i = 0; i < N ; ++i){

			// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);


		double w = 2 * arma::datum::pi / (12 * 3600);
		arma::vec::fixed<3> rotation_axis_principal_frame = arma::normalise(arma::randn<arma::vec>(3));
		arma::vec::fixed<3> Omega = w * pgm_filter -> GetPrincipalAxes().t() * rotation_axis_principal_frame;

		int N_facets = pgm_filter -> GetN_facets();
		int N_C = pgm_filter -> GetN_vertices();;

		arma::ivec f_vec = arma::randi<arma::ivec>(1,arma::distr_param(0,N_facets - 1));
		int f = f_vec(0);

		double slope = pgm_filter -> GetSlope(f,Omega);

		arma::rowvec dSlopedwC = shape_uq.GetPartialSlopePartialwPartialC(f,Omega);
		arma::vec deviation = 1e-1 * arma::randn<arma::vec>(N_C * 3) / pgm_filter -> GetScaleFactor();
		arma::vec dw_vector = arma::randn<arma::vec>(1) * arma::norm(Omega) / 100;
		double dw = dw_vector(0);

		shape_uq.ApplyDeviation(deviation);

		arma::vec all_deviations(1 + 3 * pgm_filter -> GetN_vertices());

		all_deviations(0) = dw;
		all_deviations.subvec(1,all_deviations.n_rows - 1) = pgm_filter -> GetScaleFactor() * deviation;

		arma::vec::fixed<3> Omega_p = (w + dw) * pgm_filter -> GetPrincipalAxes().t() * rotation_axis_principal_frame;
		
		double slope_p = pgm_filter -> GetSlope(f,Omega_p);

		double dslope = slope_p - slope;
		double dslope_lin = arma::dot(dSlopedwC, all_deviations);

		if(std::abs(dslope - dslope_lin)/std::abs(dslope_lin) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestGetPartialSlopePartialwPartialC with " << double(successes) / N * 100 << " \% of successes.\n";


}



void SBGATPolyhedronGravityModelUQ::TestPartialBodyFixedAccelerationfPartialC(std::string filename , double tol, bool shape_in_meters){

	std::cout << "\t In TestPartialBodyFixedAccelerationfPartialC ...";

	// MC
	int N = 100;
	int successes = 0;
	arma::arma_rng::set_seed(0);


	#pragma omp parallel for reduction(+:successes)
	for (int i = 0; i < N ; ++i){

			// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);

		double w = 2 * arma::datum::pi / (12 * 3600);
		arma::vec::fixed<3> rotation_axis_principal_frame = arma::normalise(arma::randn<arma::vec>(3));
		arma::vec::fixed<3> Omega = w * pgm_filter -> GetPrincipalAxes().t() * rotation_axis_principal_frame;


		int N_facets = pgm_filter -> GetN_facets();
		int N_C = pgm_filter -> GetN_vertices();
		arma::ivec f_vec = arma::randi<arma::ivec>(1,arma::distr_param(0,N_facets - 1));
		int f = f_vec(0);

		arma::vec::fixed<3> body_fixed_acc = pgm_filter -> GetBodyFixedAccelerationf(f,Omega);

		arma::mat dAdC =  shape_uq.PartialBodyFixedAccelerationfPartialC(f,Omega);

		arma::vec deviation = 1e-1 * arma::randn<arma::vec>(N_C * 3) / pgm_filter -> GetScaleFactor();

		shape_uq.ApplyDeviation(deviation);

		arma::vec::fixed<3> body_fixed_acc_p = pgm_filter -> GetBodyFixedAccelerationf(f,Omega);
		
		arma::vec::fixed<3> dbody_fixed_acc = body_fixed_acc_p - body_fixed_acc;
		arma::vec::fixed<3> dbody_fixed_acc_lin = dAdC * deviation * pgm_filter -> GetScaleFactor();

		

		if(arma::norm(dbody_fixed_acc - dbody_fixed_acc_lin)/arma::norm(dbody_fixed_acc_lin) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestPartialBodyFixedAccelerationfPartialC with " << double(successes) / N * 100 << " \% of successes.\n";



}



void SBGATPolyhedronGravityModelUQ::TestPartialBodyFixedAccelerationfPartialwC(std::string filename , double tol, bool shape_in_meters){

	std::cout << "\t In TestPartialBodyFixedAccelerationfPartialwC ...";

	// MC
	int N = 100;
	int successes = 0;
	arma::arma_rng::set_seed(0);


	// Reading
	vtkSmartPointer<vtkOBJReader> r = vtkSmartPointer<vtkOBJReader>::New();
	r -> SetFileName(filename.c_str());
	r -> Update(); 

	// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cl =
	vtkSmartPointer<vtkCleanPolyData>::New();
	cl -> SetInputConnection (r -> GetOutputPort());
	cl -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
	cl -> Update();	


	vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
	
	mass_prop -> SetInputConnection(cl -> GetOutputPort());
	if(shape_in_meters){
		mass_prop -> SetScaleMeters();
	}
	else{
		mass_prop -> SetScaleKiloMeters();
	}

	mass_prop -> Update();
	double xmin, xmax, ymin, ymax, zmin, zmax;
	mass_prop -> GetBoundingBox( xmin, xmax, ymin, ymax, zmin, zmax);

	arma::vec::fixed<3> pos = 1.5 * arma::vec({xmax,ymax,zmax});

	#pragma omp parallel for reduction(+:successes)
	for (int i = 0; i < N ; ++i){

			// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);

		double w = 2 * arma::datum::pi / (12 * 3600);
		arma::vec::fixed<3> rotation_axis_principal_frame = arma::normalise(arma::randn<arma::vec>(3));
		arma::vec::fixed<3> Omega = w * pgm_filter -> GetPrincipalAxes().t() * rotation_axis_principal_frame;

		int N_facets = pgm_filter -> GetN_facets();
		int N_C = pgm_filter -> GetN_vertices();
		arma::ivec f_vec = arma::randi<arma::ivec>(1,arma::distr_param(0,N_facets - 1));
		int f = f_vec(0);

		arma::vec::fixed<3> body_fixed_acc = pgm_filter -> GetBodyFixedAccelerationf(f,Omega);

		arma::mat partial = shape_uq.PartialBodyFixedAccelerationfPartialOmegaC(f,Omega);
		arma::vec deviation = 1e-1 * arma::randn<arma::vec>(N_C * 3) / pgm_filter -> GetScaleFactor();

		arma::vec::fixed<3> Omega_p = Omega + arma::normalise(arma::randn<arma::vec>(3)) * arma::norm(Omega) / 100;


		shape_uq.ApplyDeviation(deviation);

		arma::vec all_deviations(3 + 3 * N_C);

		all_deviations.subvec(0,2) = Omega_p - Omega;
		all_deviations.subvec(3,all_deviations.n_rows - 1) = pgm_filter -> GetScaleFactor() * deviation;

		arma::vec::fixed<3> body_fixed_acc_p = pgm_filter -> GetBodyFixedAccelerationf(f,Omega_p);

		arma::vec::fixed<3> dbody_fixed_acc = body_fixed_acc_p - body_fixed_acc;
		arma::vec::fixed<3> dbody_fixed_acc_lin = partial * all_deviations ;

		if(arma::norm(dbody_fixed_acc - dbody_fixed_acc_lin)/arma::norm(dbody_fixed_acc_lin) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestPartialBodyFixedAccelerationfPartialwC with " << double(successes) / N * 100 << " \% of successes.\n";

}


void SBGATPolyhedronGravityModelUQ::TestPartialOmegaPartialwC(std::string input , double tol, bool shape_in_meters){



	std::cout << "\t In TestPartialOmegaPartialwC ...";

	// MC
	int N = 100;
	int successes = 0;
	arma::arma_rng::set_seed(0);

	#pragma omp parallel for reduction(+:successes)
	for (int i = 0; i < N ; ++i){

			// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(input.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);

		arma::vec::fixed<3> rotation_axis_principal_frame = arma::normalise(arma::randn<arma::vec>(3));
		double w = 2 * arma::datum::pi / (12 * 3600);

		arma::mat::fixed<3,3> PB = pgm_filter -> GetPrincipalAxes();
		arma::vec::fixed<3> Omega = w * PB.t() * rotation_axis_principal_frame;

		int N_C = pgm_filter -> GetN_vertices();

		arma::mat partial = shape_uq.PartialOmegaPartialwC(Omega);
		arma::vec deviation = 1e-1 * arma::randn<arma::vec>(N_C * 3) / pgm_filter -> GetScaleFactor();

		arma::vec dw_vector = arma::randn<arma::vec>(1) * arma::norm(Omega) / 100;
		double dw = dw_vector(0);

		shape_uq.ApplyDeviation(deviation);

		arma::mat::fixed<3,3> PB_p = pgm_filter -> GetPrincipalAxes();
		arma::vec::fixed<3> Omega_p = (w + dw) * PB_p.t() * rotation_axis_principal_frame;

		arma::vec all_deviations(1 + 3 * N_C);

		all_deviations(0) = dw;
		all_deviations.subvec(1,all_deviations.n_rows - 1) = pgm_filter -> GetScaleFactor() * deviation;

		arma::vec::fixed<3> dOmega = Omega_p - Omega;
		arma::vec::fixed<3> dOmega_lin = partial * all_deviations;

		
		if(arma::norm(dOmega - dOmega_lin)/arma::norm(dOmega_lin) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestPartialOmegaPartialwC with " << double(successes) / N * 100 << " \% of successes.\n";

}





void SBGATPolyhedronGravityModelUQ::TestPartialSlopeArgumentPartialOmegaC(std::string input , double tol, bool shape_in_meters){

	std::cout << "\t In TestPartialSlopeArgumentPartialOmegaC ...";

	// MC
	int N = 100;
	int successes = 0;
	arma::arma_rng::set_seed(0);

	#pragma omp parallel for reduction(+:successes)
	for (int i = 0; i < N ; ++i){

			// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(input.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		if(shape_in_meters){
			pgm_filter -> SetScaleMeters();
		}
		else{
			pgm_filter -> SetScaleKiloMeters();
		}
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetModel(pgm_filter);

		arma::vec::fixed<3> rotation_axis_principal_frame = arma::normalise(arma::randn<arma::vec>(3));
		double w = arma::datum::pi / (12 * 3600);
		arma::vec::fixed<3> Omega = w * pgm_filter -> GetPrincipalAxes().t() * rotation_axis_principal_frame;

		int N_facets = pgm_filter -> GetN_facets();
		arma::ivec f_vec = arma::randi<arma::ivec>(1,arma::distr_param(0,N_facets - 1));
		int f = f_vec(0);

		arma::vec::fixed<3> body_fixed_acc = pgm_filter -> GetBodyFixedAccelerationf(f,Omega);
		double slope_argument = arma::dot(arma::normalise(body_fixed_acc),
			arma::normalise(pgm_filter -> GetNonNormalizedFacetNormal(f)));

		arma::rowvec dSlope_argumentdOmegaC = shape_uq.PartialSlopeArgumentPartialOmegaC(f,Omega,body_fixed_acc);
		arma::vec deviation = 1e-2 * arma::randn<arma::vec>(3 * pgm_filter -> GetN_vertices()) / pgm_filter -> GetScaleFactor();

		shape_uq.ApplyDeviation(deviation);

		arma::vec::fixed<3> Omega_p = Omega + arma::randn<arma::vec>(3) * arma::norm(Omega)/100;

		double slope_argument_p = arma::dot(arma::normalise(pgm_filter -> GetBodyFixedAccelerationf(f,Omega_p)),
			arma::normalise(pgm_filter -> GetNonNormalizedFacetNormal(f)));

		arma::vec all_deviations(3 + 3 * pgm_filter -> GetN_vertices());

		all_deviations.subvec(0,2) = Omega_p - Omega;
		all_deviations.subvec(3,all_deviations.n_rows - 1) = pgm_filter -> GetScaleFactor() * deviation;
		
		double dslope_argument = slope_argument_p - slope_argument;
		double dslope_argument_lin = arma::dot(dSlope_argumentdOmegaC, all_deviations);

		if(std::abs(dslope_argument - dslope_argument_lin)/std::abs(dslope_argument_lin) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestPartialSlopeArgumentPartialOmegaC with " << double(successes) / N * 100 << " \% of successes.\n";

}


double SBGATPolyhedronGravityModelUQ::  GetVarianceSlope(const int & f ,const arma::vec::fixed<3> & Omega) const{

	arma::mat augmented_PCC = arma::zeros<arma::mat>(this -> P_CC.n_rows + 1, this -> P_CC.n_rows + 1);


	augmented_PCC(0,0) = std::pow( arma::dot(Omega,Omega) / (2 * arma::datum::pi),2) * std::pow(this -> period_standard_deviation,2);
	augmented_PCC.submat(1,1,augmented_PCC.n_rows -1,augmented_PCC.n_rows -1) = this -> P_CC;

	arma::rowvec partial = this -> GetPartialSlopePartialwPartialC(f,Omega);

	return arma::dot(partial, augmented_PCC * partial.t()); 


}



