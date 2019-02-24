#include <SBGATMassPropertiesUQ.hpp>
#include <RigidBodyKinematics.hpp>
#include <vtkCleanPolyData.h>
#include <vtkOBJReader.h>

void SBGATMassPropertiesUQ::SetMassProperties(vtkSmartPointer<SBGATMassProperties> mass_prop){
	this -> mass_prop = mass_prop;
}

arma::rowvec::fixed<9> SBGATMassPropertiesUQ::PartialDeltaVfPartialTf(const int & f) const{

	double r0[3];
	double r1[3];
	double r2[3];

	arma::rowvec::fixed<9> partial;
	this -> mass_prop -> GetVerticesInFacet(f,r0,r1,r2);
	arma::vec::fixed<3> C0 = {r0[0],r0[1],r0[2]};
	arma::vec::fixed<3> C1 = {r1[0],r1[1],r1[2]};
	arma::vec::fixed<3> C2 = {r2[0],r2[1],r2[2]};

	partial.subvec(0,2) = arma::cross(C1,C2).t();
	partial.subvec(3,5) = - C0.t() * RBK::tilde(C2);
	partial.subvec(6,8) = C0.t() * RBK::tilde(C1);

	return 1./6 * partial;

	
}

arma::mat::fixed<3,9> SBGATMassPropertiesUQ::PartialDeltaCOMfPartialTf(const int & f) const{

	double volume = this -> mass_prop -> GetVolume();
	const arma::vec::fixed<3> & cm = this -> mass_prop -> GetCenterOfMass();

	arma::mat::fixed<3,9> triple_I;
	triple_I.cols(0,2) = arma::eye<arma::mat>(3,3);
	triple_I.cols(3,5) = arma::eye<arma::mat>(3,3);
	triple_I.cols(6,8) = arma::eye<arma::mat>(3,3);

	return ((- cm * this -> PartialDeltaVfPartialTf(f) + 1./4 * triple_I)/volume);

}

void SBGATMassPropertiesUQ::TestPartials(std::string input,double tol){

	SBGATMassPropertiesUQ::TestPartialDeltaVfPartialTf(input,tol);
	SBGATMassPropertiesUQ::TestPartialDeltaCOMfPartialTf(input,tol);

}

void SBGATMassPropertiesUQ::TestPartialDeltaVfPartialTf(std::string input,double tol){

	std::cout << "\t In TestPartialDeltaVfPartialTf ... ";
	int successes = 0;
	arma::arma_rng::set_seed(0);
	int N = 1000;
	#pragma omp parallel for reduction(+:successes)
	
	for (int i = 0; i < N ; ++i){

		// Reading
		vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
		reader -> SetFileName(input.c_str());
		reader -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner =
		vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner -> SetInputConnection (reader -> GetOutputPort());
		cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
		cleaner -> Update();


	// Creating the PGM dyads
		vtkSmartPointer<SBGATMassProperties> mass_prop = vtkSmartPointer<SBGATMassProperties>::New();
		mass_prop -> SetInputConnection(cleaner -> GetOutputPort());
		mass_prop -> SetDensity(1970); 
		mass_prop -> SetScaleMeters();
		mass_prop -> Update();

		SBGATMassPropertiesUQ shape_uq;
		shape_uq.SetMassProperties(mass_prop);

		int N_facets = vtkPolyData::SafeDownCast(mass_prop -> GetInput()) -> GetNumberOfCells();
		int f = arma::randi(arma::distr_param(0,N_facets - 1));

	// Nominal Volume
		double volume = shape_uq.GetMassProperties() -> GetDeltaV(f);

	// Deviation
		arma::vec::fixed<9> delta_Tf = 1e-3 * arma::randn<arma::vec>(9);

	// Linear dUf
		double dV_lin = arma::dot(shape_uq.PartialDeltaVfPartialTf(f) , delta_Tf);

	// Apply Tf deviation
		shape_uq. ApplyTfDeviation(delta_Tf,f);

	// Perturbed Uf
		double volume_p = shape_uq.GetMassProperties() -> GetDeltaV(f);
		
	// Non-linear dUf
		double dV = volume_p - volume;

		if(std::abs(dV - dV_lin) / std::abs(dV_lin) < tol){
			++successes;
		}



	}

	std::cout << "\t Passed TestPartialDeltaVfPartialTf with " << double(successes)/N * 100 << " \% of successes. \n";




}

void SBGATMassPropertiesUQ::TestPartialDeltaCOMfPartialTf(std::string input,double tol){

}



void SBGATMassPropertiesUQ::ApplyTfDeviation(arma::vec::fixed<9> delta_Tf,const int & f){
	
	int v0_index,v1_index,v2_index;
	this -> mass_prop -> GetIndicesVerticesInFacet(f,v0_index,v1_index,v2_index);

	double r0[3],r1[3],r2[3];
	vtkPolyData * polydata = vtkPolyData::SafeDownCast(this -> mass_prop -> GetInput());

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
	this -> mass_prop -> Modified();

	this -> mass_prop -> Update();

}







