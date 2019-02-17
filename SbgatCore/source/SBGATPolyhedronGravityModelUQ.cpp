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

	arma::vec::fixed<3> Nf = arma::cross(R1 - R0,R2 - R0);


	return SBGATPolyhedronGravityModelUQ::PartialFfPartialnf(arma::normalise(Nf)) * SBGATPolyhedronGravityModelUQ::PartialNormalizedVPartialNonNormalizedV(Nf) * this -> PartialNfPartialTf(f);

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

	arma::vec::fixed<3> N_e_f0 = this -> pgm_model -> GetNonNormalizedFacetNormal(f0);
	arma::vec::fixed<3> N_e_f1 = this -> pgm_model -> GetNonNormalizedFacetNormal(f1);

	arma::vec::fixed<3> n_e_f0 = arma::normalise(N_e_f0);
	arma::vec::fixed<3> n_e_f1 = arma::normalise(N_e_f1);

	arma::vec::fixed<3> Vr = this -> pgm_model -> GetEdgeFlipping(e) * RBK::tilde(r1_arma - r0_arma) * e_r;

	double le = arma::norm(r0_arma - r1_arma);

	double E_qr = arma::dot(e_q,(n_e_f1 * n_e_f1.t() - n_e_f0 * n_e_f0.t()) * Vr) / le;


	arma::mat::fixed<3,6> M = arma::zeros<arma::mat>(3,6);
	M.cols(0,2) = - arma::eye<arma::mat>(3,3);
	M.cols(3,5) = arma::eye<arma::mat>(3,3);

	arma::vec::fixed<13> first_partial;

	first_partial(0) = - E_qr;
	first_partial.subvec(1,6) = this -> pgm_model -> GetEdgeFlipping(e) * M.t() * RBK::tilde(e_r) * (arma::dot(n_e_f1,e_q) * n_e_f1 - arma::dot(n_e_f0,e_q) * n_e_f0);
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

void SBGATPolyhedronGravityModelUQ::TestPartials(double tol){
	SBGATPolyhedronGravityModelUQ::TestPartialNfPartialTf(tol);
	SBGATPolyhedronGravityModelUQ::TestPartialNormalizedVPartialNonNormalizedV(tol);	
	SBGATPolyhedronGravityModelUQ::TestPartialAtan2PartialZf(tol);
	SBGATPolyhedronGravityModelUQ::TestPartialZfPartialUnitRf(tol);
	SBGATPolyhedronGravityModelUQ::TestPartialOmegafPartialTf(tol);
	SBGATPolyhedronGravityModelUQ::TestPartialFfPartialnf(tol);
	SBGATPolyhedronGravityModelUQ::TestPartialFfPartialNonNormalizedNf(tol);
	SBGATPolyhedronGravityModelUQ::TestPartialFfPartialTf(tol);
	SBGATPolyhedronGravityModelUQ::TestPartialXfPartialTf(tol);
	SBGATPolyhedronGravityModelUQ::TestPartialUfPartialXf(tol);
	SBGATPolyhedronGravityModelUQ::TestPartialUfPartialTf(tol);

	SBGATPolyhedronGravityModelUQ::TestPartialLePartialAe(tol);

	SBGATPolyhedronGravityModelUQ::TestPartialEdgeLengthPartialAe(tol);
	SBGATPolyhedronGravityModelUQ::TestPartialEPartialBe(tol);
	SBGATPolyhedronGravityModelUQ::TestPartialUePartialXe(tol);
	SBGATPolyhedronGravityModelUQ::TestPartialUePartialBe(tol);



}

void SBGATPolyhedronGravityModelUQ::TestPartialUfPartialTf(double tol){



	std::cout << "\t In TestPartialUfPartialTf ... ";
	int successes = 0;

	for (int i = 0; i < 100 ; ++i){

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
		pgm_filter -> SetDensity(1970); 
		pgm_filter -> SetScaleMeters();
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetPGM(pgm_filter);

		int N_facets = vtkPolyData::SafeDownCast(pgm_filter -> GetInput()) -> GetNumberOfCells();
		int f = arma::randi(arma::distr_param(0,N_facets - 1));

		arma::vec pos = {2,3,4};

	// Nominal Uf
		double Uf = shape_uq.GetPGMModel() -> GetUf(shape_uq.GetPGMModel() -> GetXf(pos,f));

	// Deviation
		arma::vec::fixed<9> delta_Tf = 1e-3 * arma::randn<arma::vec>(9);

	// Linear dUf
		double dUf_lin = arma::dot(shape_uq.PartialUfPartialXf(pos,f) * shape_uq.PartialXfPartialTf(pos,f), delta_Tf);

	// Apply Tf deviation
		shape_uq. ApplyTfDeviation(delta_Tf,f);

	// Perturbed Uf
		double Uf_p = shape_uq.GetPGMModel() -> GetUf(shape_uq.GetPGMModel() -> GetXf(pos,f));
		

	// Non-linear dUf
		double dUf = Uf_p - Uf;

		if(std::abs(dUf - dUf_lin) / std::abs(dUf_lin) < tol){
			++successes;
		}



	}

	std::cout << "\t Passed TestPartialUfPartialTf with " << successes << " \% of sucesses. \n";


}
void SBGATPolyhedronGravityModelUQ::TestPartialUePartialBe(double tol){



	std::cout << "\t In TestPartialUePartialBe ... ";
	std::string filename  = "../input/cube.obj";

	int successes = 0;
	int N = 1000;
	arma::vec::fixed<3> pos = {1,2,3};
	for(int i = 0; i < N; ++i){
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
		pgm_filter -> SetDensity(1970); 
		pgm_filter -> SetScaleMeters();
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetPGM(pgm_filter);

		int N_edges = vtkPolyData::SafeDownCast(pgm_filter -> GetInput()) -> GetNumberOfCells() + vtkPolyData::SafeDownCast(pgm_filter -> GetInput()) -> GetNumberOfPoints()  -2;
		int e = arma::randi(arma::distr_param(0,N_edges - 1));

		arma::vec::fixed<6> Ee_param = shape_uq.GetPGMModel() -> GetEeParam(e);

		arma::vec deviation = 1e-3 * arma::randn<arma::vec>(3 * vtkPolyData::SafeDownCast(pgm_filter -> GetInput()) -> GetNumberOfPoints());

		arma::rowvec::fixed<24> partial = shape_uq.PartialUePartialXe(pos,e) * shape_uq.PartialXePartialBe(pos,e);


		double Ue = shape_uq.GetPGMModel() -> GetUe( shape_uq.GetPGMModel() -> GetXe(pos,e));

		// Apply global deviation and get all dBes deviation
		arma::vec dBe = shape_uq.ApplyAndGetBeDeviation(deviation);

		double dUe = shape_uq.GetPGMModel() -> GetUe( shape_uq.GetPGMModel() -> GetXe(pos,e)) - Ue;
		double dUe_lin = arma::dot(partial, dBe.subvec(24 * e, 24 * e + 23));
		
		if(std::abs(dUe - dUe_lin)/std::abs(dUe_lin) < tol){
			++successes;
		}


	}
	std::cout << "\t Passed TestPartialUePartialBe with " << double(successes)/N * 100 << " \% of successes .\n";






}


void SBGATPolyhedronGravityModelUQ::TestPartialUePartialXe(double tol){

	std::cout << "\t In TestPartialUePartialXe ... ";

	int successes = 0;
	int N = 1000;
	for (int i = 0; i < N; ++i){
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
		pgm_filter -> SetDensity(1970); 
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
		arma::vec::fixed<10> dXe = 1e-3 * arma::randn<arma::vec>(10);

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

		if(std::abs(dUe - dUe_lin) / std::abs(dUe_lin) < tol){
			++successes;
		}
	}

	std::cout << "\t Passed TestPartialUePartialXe with " << double(successes) / N * 100 <<  " \% of successes .\n";


}
void SBGATPolyhedronGravityModelUQ::TestPartialUfPartialXf(double tol){


	std::cout << "\t In TestPartialUfPartialXf ... ";

	int successes = 0;

	for (int i =0; i < 100; ++i){
		std::string filename  = "../input/cube.obj";

	// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(filename.c_str());
		reader -> Update(); 

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(reader -> GetOutputPort());
		pgm_filter -> SetDensity(1970); 
		pgm_filter -> SetScaleMeters();
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetPGM(pgm_filter);

	// Test point 
		arma::vec::fixed<3> pos = {3,2,3};

  	// Pick f
		int N_facets = vtkPolyData::SafeDownCast(pgm_filter -> GetInput()) -> GetNumberOfCells();

		int f = arma::randi(arma::distr_param(0,N_facets - 1));

  	// Edge point
		arma::vec::fixed<3> rfi_0 = shape_uq.GetPGMModel() -> GetRf(pos,f);

  	// Xe before
		arma::vec::fixed<10> Xf = shape_uq.GetPGMModel() -> GetXf(pos,f);

  	// Ue before
		double Uf = shape_uq.GetPGMModel() -> GetUf(Xf);

  	// Perturbation to Xe 
		arma::vec::fixed<10> dXf = 1e-3 * arma::randn<arma::vec>(10);

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

		if(std::abs(dUf - dUf_lin) / std::abs(dUf_lin) < tol){
			++successes;
		}
	}


	std::cout << "\t Passed TestPartialUfPartialXf with " << successes<< "\% of successes. \n";


}

void SBGATPolyhedronGravityModelUQ::TestPartialXfPartialTf(double tol){


	std::cout << "\t In TestPartialXfPartialTf ... ";
	int successes = 0;

	for (int i = 0; i < 100 ; ++i){

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
		pgm_filter -> SetDensity(1970); 
		pgm_filter -> SetScaleMeters();
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetPGM(pgm_filter);
		int N_facets = vtkPolyData::SafeDownCast(pgm_filter -> GetInput()) -> GetNumberOfCells();

		int f = arma::randi(arma::distr_param(0,N_facets - 1));
		arma::vec pos = {2,3,4};

	// Nominal Xf
		arma::vec::fixed<10> Xf;
		Xf(0) = pgm_filter -> GetOmegaf(pos,f);
		Xf.subvec(1,3) = pgm_filter -> GetRf(pos,f);
		Xf.subvec(4,9) = pgm_filter -> GetFfParam(f);

	// Deviation
		arma::vec::fixed<9> delta_Tf = 1e-3 * arma::randn<arma::vec>(9);

	// Linear dXf
		arma::vec::fixed<10> dXf_lin = shape_uq.PartialXfPartialTf(pos,f) * delta_Tf;

	// Apply Tf deviation
		shape_uq. ApplyTfDeviation(delta_Tf,f);

	// Perturbed Xf
		arma::vec::fixed<10> Xf_p;
		Xf_p(0) = pgm_filter -> GetOmegaf(pos,f);
		Xf_p.subvec(1,3) = pgm_filter -> GetRf(pos,f);
		Xf_p.subvec(4,9) = pgm_filter -> GetFfParam(f);

	// Non-linear dXf
		arma::vec::fixed<10> dXf = Xf_p - Xf;

		if(arma::norm(dXf - dXf_lin) / arma::norm(dXf_lin) < tol){
			++successes;
		}



	}

	std::cout << "\t Passed TestPartialXfPartialTf with " << successes << " \% of sucesses. \n";

}
void SBGATPolyhedronGravityModelUQ::TestPartialOmegafPartialTf(double tol){


	std::cout << "\t In TestPartialOmegafPartialTf ... ";

	int successes = 0;

	for (int i = 0; i < 100; ++i){
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

		vtkSmartPointer<vtkPolyData> polydata = cleaner -> GetOutput();
		
	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputData(polydata);
		pgm_filter -> SetDensity(1970); 
		pgm_filter -> SetScaleMeters();
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetPGM(pgm_filter);
		int N_facets = vtkPolyData::SafeDownCast(pgm_filter -> GetInput()) -> GetNumberOfCells();

		int f = arma::randi(arma::distr_param(0,N_facets - 1));
		arma::vec pos = {2,3,4};

	// Nominal omega_f
		double omega_f = pgm_filter -> GetOmegaf(pos,f);

	// Deviation
		arma::vec::fixed<9> delta_Tf = 1e-3 * arma::randn<arma::vec>(9);

	// Linear dXf
		double domega_f_lin = arma::dot(shape_uq.PartialOmegafPartialTf(pos,f), delta_Tf);

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
void SBGATPolyhedronGravityModelUQ::TestPartialFfPartialTf(double tol){


	std::cout << "\t In TestPartialFfPartialTf ... ";
	int successes = 0;
	for (int i = 0; i < 100; ++i){
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
		pgm_filter -> SetDensity(1970); 
		pgm_filter -> SetScaleMeters();
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetPGM(pgm_filter);
		int N_facets = vtkPolyData::SafeDownCast(pgm_filter -> GetInput()) -> GetNumberOfCells();

		int f = arma::randi(arma::distr_param(0,N_facets - 1));

	// Nominal dyad
		arma::vec::fixed<6> Ff = pgm_filter -> GetFfParam(f);

	// Deviation
		arma::vec::fixed<9> delta_Tf = 1e-3 * arma::randn<arma::vec>(9);

	// Linear difference
		arma::vec::fixed<6> dFf_lin = shape_uq.PartialFfPartialTf(f) * delta_Tf;

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


void SBGATPolyhedronGravityModelUQ::TestPartialNormalizedVPartialNonNormalizedV(double tol){

	std::cout << "\t In TestPartialNormalizedVPartialNonNormalizedV ... ";

	int N = 1000;
	int successes = 0;
	for (int i = 0; i < N ;++i){
  	// non-normalized vector
		arma::vec::fixed<3> X = arma::randn<arma::vec>(3);

  	// normalized vector
		arma::vec::fixed<3> x = arma::normalise(X);


  	// disturbance 
		arma::vec::fixed<3> dX = 1e-3 * arma::randn<arma::vec>(3);

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
void SBGATPolyhedronGravityModelUQ::TestPartialNfPartialTf(double tol){

	std::cout << "\t In TestPartialNfPartialTf ... ";

	std::string filename  = "../input/cube.obj";

	int successes = 0;
	for (int i = 0; i < 100; ++i){
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
		pgm_filter -> SetDensity(1970); 
		pgm_filter -> SetScaleMeters();
		pgm_filter -> Update();


		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetPGM(pgm_filter);

		int N_facets = vtkPolyData::SafeDownCast(pgm_filter -> GetInput()) -> GetNumberOfCells();

		int f = arma::randi(arma::distr_param(0,N_facets - 1));
		arma::vec::fixed<9> delta_Tf = 1e-3 * arma::randn<arma::vec>(9);

		arma::vec::fixed<3> N = shape_uq.GetPGMModel() -> GetNonNormalizedFacetNormal(f);

	// Linear difference
		arma::vec::fixed<3> dN_lin = shape_uq . PartialNfPartialTf(f) * delta_Tf;

	// Apply Tf deviation
		shape_uq. ApplyTfDeviation(delta_Tf,f);

	// Perturbed normal
		arma::vec::fixed<3> N_p = shape_uq.GetPGMModel() -> GetNonNormalizedFacetNormal(f);
		
	// Non-linear difference
		arma::vec::fixed<3> dN = N_p - N;

		if(arma::norm(dN - dN_lin)/arma::norm(dN_lin) < tol){
			++successes;
		}

	}

	std::cout << "\t Passed TestPartialNfPartialTf with " << successes <<" \% of successes. \n";



}
void SBGATPolyhedronGravityModelUQ::TestPartialFfPartialnf(double tol){

	std::cout << "\t In TestPartialFfPartialnf ... ";
	int N = 1000;
	int successes = 0;
	for (int i = 0; i < N; ++i){

		arma::vec::fixed<3> nf = {1,2,-3};
		nf = arma::normalise(nf);


	// Dyad
		arma::mat::fixed<3,3> Ff = nf * nf.t();

	// Dyad parametrization
		arma::vec::fixed<6> Ff_vec = {Ff(0,0),Ff(1,1),Ff(2,2),Ff(0,1),Ff(0,2),Ff(1,2)};

	// Perturbation

		arma::vec::fixed<3> dnf = 1e-3 * arma::randn<arma::vec>(3);

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

void SBGATPolyhedronGravityModelUQ::TestPartialFfPartialNonNormalizedNf(double tol){

	std::cout << "\t In TestPartialFfPartialNonNormalizedNf ... ";

	int N = 1000;
	int successes = 0;
	for (int i =0; i < N; ++i){
		arma::vec::fixed<3> Nf = {1,2,3};
		arma::vec::fixed<3> nf = arma::normalise(Nf);


	// Dyad
		arma::mat::fixed<3,3> Ff = nf * nf.t();

	// Dyad parametrization
		arma::vec::fixed<6> Ff_vec = {Ff(0,0),Ff(1,1),Ff(2,2),Ff(0,1),Ff(0,2),Ff(1,2)};

	// Perturbation

		arma::vec::fixed<3> dNf = 1e-3 * arma::randn<arma::vec>(3);

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







void SBGATPolyhedronGravityModelUQ::TestPartialLePartialAe(double tol){


	std::cout << "\t In TestPartialLePartialAe ... ";
	std::string filename  = "../input/cube.obj";

	int successes = 0;
	arma::vec pos = {1,2,3};

	for (int i = 0; i < 100; ++i){
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
		pgm_filter -> SetDensity(1970); 
		pgm_filter -> SetScaleMeters();
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetPGM(pgm_filter);

		int e = 1;
		arma::vec::fixed<6> delta_Ae = 1e-3 * arma::randn<arma::vec>(6);
	
		double Le = shape_uq .GetPGMModel() -> GetLe(pos,e);
		

	// Apply Ae deviation
		shape_uq.ApplyAeDeviation(delta_Ae,e);

	// Perturbed length
		double Le_p = shape_uq .GetPGMModel() -> GetLe(pos,e);


	// Non-linear difference
		double dLe = Le_p - Le;

	// Linear difference
		double dLe_lin = arma::dot(shape_uq.PartialLePartialAe(pos,e), delta_Ae);


		if(std::abs(dLe - dLe_lin)/std::abs(dLe_lin) < tol){
			++successes;
		}
	}
	std::cout << "\t Passed TestPartialLePartialAe with " << successes<< " \% of successes. \n";


}


void SBGATPolyhedronGravityModelUQ::TestPartialEdgeLengthPartialAe(double tol){

	std::cout << "\t In TestPartialEdgeLengthPartialAe ... ";
	std::string filename  = "../input/cube.obj";

	int successes = 0;

	for (int i = 0; i < 100; ++i){
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
		pgm_filter -> SetDensity(1970); 
		pgm_filter -> SetScaleMeters();
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetPGM(pgm_filter);

		int e = 1;
		arma::vec::fixed<6> delta_Ae = 1e-3 * arma::randn<arma::vec>(6);

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


		if(std::abs(dle - dle_lin)/std::abs(dle_lin) < tol){
			++successes;
		}
	}
	std::cout << "\t Passed TestPartialEdgeLengthPartialAe with " << successes<< " \% of successes. \n";

}
void SBGATPolyhedronGravityModelUQ::TestPartialEPartialBe(double tol){


	std::cout << "\t In TestPartialEPartialBe ... ";
	std::string filename  = "../input/cube.obj";

	int successes = 0;
	int N = 1000;

	for(int i = 0; i < N; ++i){
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
		pgm_filter -> SetDensity(1970); 
		pgm_filter -> SetScaleMeters();
		pgm_filter -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq;
		shape_uq.SetPGM(pgm_filter);

		int N_edges = vtkPolyData::SafeDownCast(pgm_filter -> GetInput()) -> GetNumberOfCells() + vtkPolyData::SafeDownCast(pgm_filter -> GetInput()) -> GetNumberOfPoints()  -2;
		int e = arma::randi(arma::distr_param(0,N_edges - 1));


		arma::vec::fixed<6> Ee_param = shape_uq.GetPGMModel() -> GetEeParam(e);

		arma::vec deviation = 1e-3 * arma::randn<arma::vec>(3 * vtkPolyData::SafeDownCast(pgm_filter -> GetInput()) -> GetNumberOfPoints());

		auto partial = shape_uq.PartialEPartialBe(e);


		arma::vec::fixed<24> Be_before = arma::zeros<arma::vec>(24);
		int f0_e,f1_e;

		double r0_e[3],r1_e[3];
		double r0_f0_e[3],r1_f0_e[3], r2_f0_e[3];
		double r0_f1_e[3],r1_f1_e[3], r2_f1_e[3];

		shape_uq.GetPGMModel() -> GetVerticesOnEdge(e,r0_e,r1_e);
		shape_uq.GetPGMModel() -> GetIndicesOfAdjacentFacets(e,f0_e,f1_e);
		shape_uq.GetPGMModel() -> GetVerticesInFacet(f0_e,r0_f0_e,r1_f0_e,r2_f0_e);
		shape_uq.GetPGMModel() -> GetVerticesInFacet(f1_e,r0_f1_e,r1_f1_e,r2_f1_e);

		Be_before.subvec(0,2) = arma::vec::fixed<3>({r0_e[0],r0_e[1],r0_e[2]});
		Be_before.subvec(3,5) = arma::vec::fixed<3>({r1_e[0],r1_e[1],r1_e[2]});

		Be_before.subvec(6,8) = arma::vec::fixed<3>({r0_f0_e[0],r0_f0_e[1],r0_f0_e[2]});
		Be_before.subvec(9,11) = arma::vec::fixed<3>({r1_f0_e[0],r1_f0_e[1],r1_f0_e[2]});
		Be_before.subvec(12,14) = arma::vec::fixed<3>({r2_f0_e[0],r2_f0_e[1],r2_f0_e[2]});

		Be_before.subvec(15,17) = arma::vec::fixed<3>({r0_f1_e[0],r0_f1_e[1],r0_f1_e[2]});
		Be_before.subvec(18,20) = arma::vec::fixed<3>({r1_f1_e[0],r1_f1_e[1],r1_f1_e[2]});
		Be_before.subvec(21,23) = arma::vec::fixed<3>({r2_f1_e[0],r2_f1_e[1],r2_f1_e[2]});

		// Apply global deviation and get all dBes deviation
		arma::vec dBe = shape_uq.ApplyAndGetBeDeviation(deviation);

		arma::vec::fixed<6> Ee_param_p = shape_uq.GetPGMModel() -> GetEeParam(e);

		// Difference
		arma::vec::fixed<24> Be_after = arma::zeros<arma::vec>(24);

		double r0_e_p[3],r1_e_p[3];
		double r0_f0_e_p[3],r1_f0_e_p[3], r2_f0_e_p[3];
		double r0_f1_e_p[3],r1_f1_e_p[3], r2_f1_e_p[3];
		
		shape_uq.GetPGMModel() -> GetVerticesOnEdge(e,r0_e_p,r1_e_p);
		shape_uq.GetPGMModel() -> GetVerticesInFacet(f0_e,r0_f0_e_p,r1_f0_e_p,r2_f0_e_p);
		shape_uq.GetPGMModel() -> GetVerticesInFacet(f1_e,r0_f1_e_p,r1_f1_e_p,r2_f1_e_p);

		Be_after.subvec(0,2) = arma::vec::fixed<3>({r0_e_p[0],r0_e_p[1],r0_e_p[2]});
		Be_after.subvec(3,5) = arma::vec::fixed<3>({r1_e_p[0],r1_e_p[1],r1_e_p[2]});

		Be_after.subvec(6,8) = arma::vec::fixed<3>({r0_f0_e_p[0],r0_f0_e_p[1],r0_f0_e_p[2]});
		Be_after.subvec(9,11) = arma::vec::fixed<3>({r1_f0_e_p[0],r1_f0_e_p[1],r1_f0_e_p[2]});
		Be_after.subvec(12,14) = arma::vec::fixed<3>({r2_f0_e_p[0],r2_f0_e_p[1],r2_f0_e_p[2]});

		Be_after.subvec(15,17) = arma::vec::fixed<3>({r0_f1_e_p[0],r0_f1_e_p[1],r0_f1_e_p[2]});
		Be_after.subvec(18,20) = arma::vec::fixed<3>({r1_f1_e_p[0],r1_f1_e_p[1],r1_f1_e_p[2]});
		Be_after.subvec(21,23) = arma::vec::fixed<3>({r2_f1_e_p[0],r2_f1_e_p[1],r2_f1_e_p[2]});



		arma::vec::fixed<6> dEe_param = Ee_param_p - Ee_param;
		arma::vec::fixed<6> dEe_param_lin = partial * dBe.subvec(24 * e,24 * e + 23);


		assert( arma::norm(Be_after - Be_before - dBe.subvec(24 * e,24 * e + 23)) /arma::norm(dBe.subvec(24 * e,24 * e + 23)) < tol);

		if(arma::norm(dEe_param - dEe_param_lin)/arma::norm(dEe_param_lin) < tol){
			++successes;
		}


	}
	std::cout << "\t Passed TestPartialEPartialBe with " << double(successes)/N * 100 << " \% of successes .\n";



}


arma::sp_mat  SBGATPolyhedronGravityModelUQ::PartialBePartialC(int e) const{

	arma::sp_mat table(24, 3 * vtkPolyData::SafeDownCast(this -> pgm_model -> GetInput()) -> GetNumberOfPoints());

	// Ae
	int v0_e,v1_e;
	this -> pgm_model -> GetIndicesVerticesOnEdge(e,v0_e,v1_e);
	table.submat(0,3 * v0_e, 2,3 * v0_e + 2) = arma::eye<arma::mat>(3,3);
	table.submat(3,3 * v1_e, 5,3 * v1_e + 2) = arma::eye<arma::mat>(3,3);

	// Ts
	int f0_e,f1_e;
	this -> pgm_model -> GetIndicesOfAdjacentFacets(e,f0_e,f1_e);

	// T0
	int v0_f0_e,v1_f0_e,v2_f0_e;
	this -> pgm_model -> GetIndicesVerticesInFacet(f0_e, v0_f0_e,v1_f0_e,v2_f0_e);

	table.submat(6,3 * v0_f0_e, 8,3 * v0_f0_e + 2) = arma::eye<arma::mat>(3,3);
	table.submat(9,3 * v1_f0_e, 11,3 * v1_f0_e + 2) = arma::eye<arma::mat>(3,3);
	table.submat(12,3 * v2_f0_e, 14,3 * v2_f0_e + 2) = arma::eye<arma::mat>(3,3);

	// T1
	int v0_f1_e,v1_f1_e,v2_f1_e;
	this -> pgm_model -> GetIndicesVerticesInFacet(f1_e, v0_f1_e,v1_f1_e,v2_f1_e);

	table.submat(15,3 * v0_f1_e, 17,3 * v0_f1_e + 2) = arma::eye<arma::mat>(3,3);
	table.submat(18,3 * v1_f1_e, 20,3 * v1_f1_e + 2) = arma::eye<arma::mat>(3,3);
	table.submat(21,3 * v2_f1_e, 23,3 * v2_f1_e + 2) = arma::eye<arma::mat>(3,3);

	return table;

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



arma::vec SBGATPolyhedronGravityModelUQ::ApplyAndGetBeDeviation(const arma::vec & delta){
	


	vtkPolyData * polydata = vtkPolyData::SafeDownCast(this -> pgm_model -> GetInput());
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
	this -> pgm_model -> Modified();

	this -> pgm_model -> Update();

	int N_edges = polydata -> GetNumberOfPoints() + polydata -> GetNumberOfCells() - 2;

	arma::vec be_deviation(24 * N_edges);
	for (int e = 0 ; e < N_edges; ++e){
		be_deviation.subvec(24 * e,24 * e + 23) = this -> PartialBePartialC(e) * delta;
	}

	return be_deviation;

}


void SBGATPolyhedronGravityModelUQ::TestPartialAtan2PartialZf(double tol){

	std::cout << "\t In TestPartialAtan2PartialZf ... ";

	int N = 1000;
	int successes = 0;
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


void SBGATPolyhedronGravityModelUQ::TestPartialZfPartialUnitRf(double tol){

	std::cout << "\t In TestPartialZfPartialUnitRf ...";
	int N = 1000;
	int successes = 0;

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


		arma::vec::fixed<9> dRf = 1e-3 * arma::randn<arma::vec>(9);

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


