
#include <SBGATFilterUQ.hpp>
#include <json.hpp>
#include <RigidBodyKinematics.hpp>
#include <vtkPolyDataNormals.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>

arma::sp_mat  SBGATFilterUQ::PartialTfPartialC(const int & f) const{

	arma::sp_mat table(9, 3 * this -> model -> GetN_vertices());

	
	int v0_f,v1_f,v2_f;
	this -> model -> GetIndicesVerticesInFacet(f, v0_f,v1_f,v2_f);

	table.submat(0,3 * v0_f, 2,3 * v0_f + 2) = arma::eye<arma::mat>(3,3);
	table.submat(3,3 * v1_f, 5,3 * v1_f + 2) = arma::eye<arma::mat>(3,3);
	table.submat(6,3 * v2_f, 8,3 * v2_f + 2) = arma::eye<arma::mat>(3,3);

	return table;

}


arma::mat SBGATFilterUQ::GetCovarianceSquareRoot(bool use_cholesky) const{

	if (use_cholesky){

		try{
			return arma::chol(this -> P_CC,"lower") / this -> model ->  GetScaleFactor() ;
		}
		catch(std::runtime_error & e){
			return arma::eye<arma::mat>(this -> P_CC.n_rows,this -> P_CC.n_rows);
		}
	}
	else{
		try{
			arma::vec eigenvalues;
			arma::mat eigenvector;

			arma::eig_sym(eigenvalues,eigenvector,this -> P_CC);

			for (int e =0; e < eigenvalues.size(); ++e){
				if(eigenvalues(e) < 0){
					eigenvalues(e) = 0;
				}
			}


			return eigenvector * arma::diagmat(arma::sqrt(eigenvalues)) /  this -> model ->  GetScaleFactor() * eigenvector.t();
		}
		catch(std::runtime_error & e){
			return arma::eye<arma::mat>(this -> P_CC.n_rows,this -> P_CC.n_rows);
		}

	}


}

void SBGATFilterUQ::SetCovarianceComponent(const arma::mat::fixed<3,3> & P,const int & v0, const int & v1){

	this -> P_CC.submat(3 * v0,3 * v1,3 * v0 + 2,3 * v1 + 2) = P * std::pow(this -> model ->  GetScaleFactor(),2);

	// this -> P_CC_sparse.submat(3 * v0,3 * v1,3 * v0 + 2,3 * v1 + 2) = this -> P_CC.submat(3 * v0,3 * v1,3 * v0 + 2,3 * v1 + 2);

}


void SBGATFilterUQ::SaveVerticesCovariance(std::string path) const{
	this -> P_CC.save(path,arma::raw_ascii);
}

arma::mat SBGATFilterUQ::GetVerticesCovariance() const{
	return this -> P_CC / std::pow(this -> model -> GetScaleFactor(),2);
}

void SBGATFilterUQ::SetCovariance(const arma::mat & P_CC){
	this -> P_CC = P_CC * std::pow(this -> model -> GetScaleFactor(),2);
}

int SBGATFilterUQ::LoadVerticesCovarianceFromJson(std::string path){


	// The JSON container is created
	nlohmann::json vertices_covariances_json;

  // The file is loaded into the container
	std::ifstream i(path);  
	i >> vertices_covariances_json;

	int N_C = this -> model -> GetN_vertices(); 


	if (vertices_covariances_json["N_C"] != N_C){
		return 0;
	}
	else{
		this -> P_CC.fill(0);

		nlohmann::json covariance_partitions_json = vertices_covariances_json.at("COVARIANCE_PARTITIONS");
		for (auto covariance_partition : covariance_partitions_json){

			int i = covariance_partition["i"];
			int j = covariance_partition["j"];

			std::vector<double> P_CiCj_vector = {
				covariance_partition["value"][0],
				covariance_partition["value"][1],
				covariance_partition["value"][2],
				covariance_partition["value"][3],
				covariance_partition["value"][4],
				covariance_partition["value"][5],
				covariance_partition["value"][6],
				covariance_partition["value"][7],
				covariance_partition["value"][8]
			};


			this -> P_CC.submat(3 * i, 3 * j, 3 * i + 2, 3 * j + 2) = arma::mat::fixed<3,3>({
				{P_CiCj_vector[0],P_CiCj_vector[3],P_CiCj_vector[6]},
				{P_CiCj_vector[1],P_CiCj_vector[4],P_CiCj_vector[7]},
				{P_CiCj_vector[2],P_CiCj_vector[5],P_CiCj_vector[8]}

			});
			if (i != j){
				this -> P_CC.submat(3 * j, 3 * i, 3 * j + 2, 3 * i + 2) = this -> P_CC.submat(3 * i, 3 * j, 3 * i + 2, 3 * j + 2).t();
			}

		}
	}

	return 1;

}





void SBGATFilterUQ::SaveNonZeroVerticesCovariance(std::string path) const {


	nlohmann::json vertices_covariances_json;
	
	int N_C = this -> model -> GetN_vertices(); 
	vertices_covariances_json["N_C"] = N_C;

	nlohmann::json covariance_partitions;
	

	for (int i = 0; i < N_C; ++i){

		const arma::mat::fixed<3,3> & P_CiCi = this -> P_CC.submat(3 * i, 3 * i,3 * i + 2, 3 * i + 2 );

		if(arma::abs(P_CiCi).max() > 1e-10){
			arma::vec::fixed<9> P_CiCi_vectorized = arma::vectorise( P_CiCi ) ;
			std::vector<double> P_CiCi_vector = {
				P_CiCi_vectorized(0),
				P_CiCi_vectorized(1),
				P_CiCi_vectorized(2),
				P_CiCi_vectorized(3),
				P_CiCi_vectorized(4),
				P_CiCi_vectorized(5),
				P_CiCi_vectorized(6),
				P_CiCi_vectorized(7),
				P_CiCi_vectorized(8)
			};
			nlohmann::json cov_CiCi = { 
				{"i", i}, 
				{"j", i}, 
				{"value", P_CiCi_vector},
			};
			covariance_partitions.push_back(cov_CiCi);
		}

		for (int j = i + 1; j < N_C; ++j){

			const arma::mat::fixed<3,3> & P_CiCj = this -> P_CC.submat(3 * i, 3 * j,3 * i + 2, 3 * j + 2 );

			if(arma::abs(P_CiCj).max() > 1e-10){
				arma::vec::fixed<9> P_CiCj_vectorized = arma::vectorise( P_CiCj ) ;
				std::vector<double> P_CiCj_vector = {
					P_CiCj_vectorized(0),
					P_CiCj_vectorized(1),
					P_CiCj_vectorized(2),
					P_CiCj_vectorized(3),
					P_CiCj_vectorized(4),
					P_CiCj_vectorized(5),
					P_CiCj_vectorized(6),
					P_CiCj_vectorized(7),
					P_CiCj_vectorized(8)
				};
				nlohmann::json cov_CiCj = { 
					{"i", i}, 
					{"j", j}, 
					{"value", P_CiCj_vector},
				};
				covariance_partitions.push_back(cov_CiCj);
			}
			

		}
		

	}

	vertices_covariances_json["COVARIANCE_PARTITIONS"] = covariance_partitions;
	std::ofstream o(path);
	o << std::setw(4) << vertices_covariances_json << std::endl;


}



arma::mat::fixed<3,9> SBGATFilterUQ::PartialNfPartialTf(const int & f) const{

	arma::mat::fixed<3,9> partial;

	double r0[3], r1[3], r2[3];

	this -> model -> GetVerticesInFacet(f,r0,r1,r2);

	arma::vec::fixed<3> r0_arma = {r0[0],r0[1],r0[2]};
	arma::vec::fixed<3> r1_arma = {r1[0],r1[1],r1[2]};
	arma::vec::fixed<3> r2_arma = {r2[0],r2[1],r2[2]};

	partial.cols(0,2) = RBK::tilde(r2_arma - r1_arma);
	partial.cols(3,5) = RBK::tilde(r0_arma - r2_arma);
	partial.cols(6,8) = RBK::tilde(r1_arma - r0_arma);

	return partial * this -> model -> GetScaleFactor();

}


void SBGATFilterUQ::ComputeVerticesCovarianceGlobal(const double & standard_dev,const double & correl_distance){

	double epsilon = 1e-4;

	vtkPolyData * input = vtkPolyData::SafeDownCast(this -> model -> GetInput());
	int N_C = input -> GetNumberOfPoints();
	
	vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();

	normalGenerator -> SetInputData(input);
	normalGenerator -> ComputePointNormalsOn();
	normalGenerator -> ComputeCellNormalsOff();
	normalGenerator -> SplittingOff ();
	normalGenerator -> Update();

	vtkPolyData * input_with_normals = normalGenerator -> GetOutput();

	vtkFloatArray * normals =  vtkFloatArray::SafeDownCast(input_with_normals->GetPointData()->GetArray("Normals"));

	assert(input_with_normals -> GetNumberOfPoints() == N_C);
	assert(normals -> GetNumberOfTuples() == N_C);

	this -> P_CC.fill(0);
	
	for (unsigned int i = 0; i < N_C; ++i){

		double ni_[3];
		double Pi_[3];
		input -> GetPoint(i,Pi_);
		normals -> GetTuple(i,ni_);
		arma::vec::fixed<3> ni = {ni_[0],ni_[1],ni_[2]};
		arma::vec::fixed<3> u_2 = arma::normalise(arma::cross(ni,arma::randn<arma::vec>(3)));
		arma::vec u_1 = arma::cross(u_2,ni);
		arma::mat::fixed<3,3> P = std::pow(standard_dev,2) * (ni * ni.t() + epsilon * (u_1 * u_1.t() + u_2 * u_2.t()));

		this -> P_CC.submat(3 * i, 3 * i, 3 * i + 2, 3 * i + 2) = P;

		for (unsigned int j = i + 1; j < N_C; ++j){

			double nj_[3];
			double Pj_[3];
			input -> GetPoint(j,Pj_);
			normals -> GetTuple(j,nj_);
			arma::vec::fixed<3> nj = {nj_[0],nj_[1],nj_[2]};

			double distance = this -> model -> GetScaleFactor() * std::sqrt(vtkMath::Distance2BetweenPoints(Pj_,Pi_));

			if ( distance < 3 * correl_distance){
				double decay = std::exp(- std::pow(distance / correl_distance,2)) ;

				arma::mat::fixed<3,3> P_correlation = std::pow(standard_dev,2) * decay * ni * nj.t();

				this -> P_CC.submat(3 * i, 3 * j, 3 * i + 2, 3 * j + 2) = P_correlation;
				this -> P_CC.submat(3 * j, 3 * i, 3 * j + 2, 3 * i + 2) = P_correlation.t();

			}

		}

	}

}

