
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


arma::mat SBGATFilterUQ::GetCovarianceSquareRoot(std::string method) const{

	if (method == "chol"){

		try{
			return arma::chol(this -> P_CC,"lower") / this -> model ->  GetScaleFactor() ;
		}
		catch(std::runtime_error & e){
			return this -> GetCovarianceSquareRoot("eigen");
		}
	}
	else if (method == "eigen"){
		try{
			arma::vec eigenvalues;
			arma::mat eigenvector;
			arma::eig_sym(eigenvalues,eigenvector,this -> P_CC,"std");
			for (int e = 0; e < eigenvalues.size(); ++e){
				if(eigenvalues(e) < 0){
					eigenvalues(e) = 0;
				}
			}

			return eigenvector * arma::diagmat(arma::sqrt(eigenvalues)) /  this -> model ->  GetScaleFactor() * eigenvector.t();

			
			
		}
		catch(std::runtime_error & e){
			return arma::zeros<arma::mat>(this -> P_CC.n_rows,this -> P_CC.n_rows);
		}

	}
	else if(method == "svd"){

		try{
			
			arma::mat U,V;
			arma::vec s;

			arma::svd(U,s,V, this -> P_CC ) ;


			return V * (arma::diagmat(1./s)/  this -> model ->  GetScaleFactor()) * V.t();
			
		}
		catch(std::runtime_error & e){
			return arma::zeros<arma::mat>(this -> P_CC.n_rows,this -> P_CC.n_rows);
		}

	}


}

void SBGATFilterUQ::SetCovarianceComponent(const arma::mat::fixed<3,3> & P,const int & v0, const int & v1){

	if (this -> P_CC.n_rows != 3 * this -> model -> GetN_vertices()){
		this -> P_CC.clear();
		this -> P_CC = arma::zeros<arma::mat>(3 * this -> model -> GetN_vertices(),3 * this -> model -> GetN_vertices());
	}


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

	return partial ;

}


void SBGATFilterUQ::ComputeVerticesCovarianceGlobal(const double & standard_dev,const double & correl_distance){


	if (this -> P_CC.n_rows != 3 * this -> model -> GetN_vertices()){
		this -> P_CC.clear();
		this -> P_CC = arma::zeros<arma::mat>(3 * this -> model -> GetN_vertices(),3 * this -> model -> GetN_vertices());
	}
	

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


void SBGATFilterUQ::TakeAndSaveSlice(int axis,std::string path, const double & c) const {
	std::vector<std::vector<arma::vec> > lines;
	this -> TakeSlice(axis,lines,c);
	this -> SaveSlice(axis,path,lines);
}


void SBGATFilterUQ::SaveSlice(int axis, 
	std::string path, const std::vector<std::vector<arma::vec> > & lines) const{

	int a_1,a_2;

	if (axis == 0){
		a_1 = 1;
		a_2 = 2;
	}
	else if (axis == 1){
		a_1 = 0;
		a_2 = 2;
	}
	else if (axis == 2 ){
		a_1 = 0;
		a_2 = 1;	
	}
	else{
		a_1 = a_2 = 0;
		throw(std::runtime_error("Specified incorrect axis: " + std::to_string(axis)));
	}



	arma::mat lines_arma;

	if (lines.size() > 0){
		lines_arma = arma::mat(lines.size(),4);

		for (int i = 0; i < lines.size() ; ++i){

			arma::rowvec rowvec = {lines[i][0](a_1),lines[i][0](a_2),lines[i][1](a_1),lines[i][1](a_2)};
			lines_arma.row(i) = rowvec;
		}

		lines_arma.save(path, arma::raw_ascii);

	}

}



void SBGATFilterUQ::TakeSlice(int axis,
	std::vector<std::vector<arma::vec> > & lines,
	const double & c) const{

	arma::vec n_plane = {0,0,0};
	n_plane(axis) = 1;


	arma::mat::fixed<3,2> T = {
		{1,0},
		{0,1},
		{-1,-1}
	};

	arma::vec::fixed<3> e3 = {0,0,1};
	arma::mat::fixed<3,3> C;

	// Each surface element is "sliced"
	for (int el = 0; el < this -> model -> GetN_facets(); ++el){


		double rf0[3];
		double rf1[3];
		double rf2[3];

		int p0,p1,p2;
		this -> model -> GetIndicesVerticesInFacet(el,p0,p1,p2);
		this -> model -> GetVerticesInFacet(el,rf0,rf1,rf2);

		C.col(0) = arma::vec({rf0[0],rf0[1],rf0[2]});
		C.col(1) = arma::vec({rf1[0],rf1[1],rf1[2]});
		C.col(2) = arma::vec({rf2[0],rf2[1],rf2[2]});

		arma::rowvec M = n_plane.t() * C * T;
		double e = c - arma::dot(n_plane,C * e3);
		arma::vec intersect;

		std::vector <arma::vec> intersects;

		// Looking for an intersect along the u = 0 edge
		if (std::abs(M(1)) > 1e-16){
			double v_intersect = e / M(1);
			if (v_intersect >= 0 && v_intersect <= 1 ){
				arma::vec Y = {0,v_intersect};

				intersect = C * (T * Y + e3);
				intersects.push_back(intersect);
			}
		}

		// Looking for an intersect along the v = 0 edge
		if (std::abs(M(0)) > 1e-16){
			double u_intersect = e / M(0);
			if (u_intersect >= 0 && u_intersect <= 1 ){
				arma::vec Y = {u_intersect,0};

				intersect = C * (T * Y + e3);
				intersects.push_back(intersect);
			}
		}

		// Looking for an intersect along the w = 0 edge
		// using u as the parameter

		if (std::abs(M(0) - M(1)) > 1e-16 && intersects.size() < 2){
			double u_intersect = (e - M(1)) / (M(0) - M(1));
			if (u_intersect >= 0 && u_intersect <= 1 ){
				arma::vec Y = {u_intersect,1 - u_intersect};

				intersect = C * (T * Y + e3);
				intersects.push_back(intersect);
			}
		}

		if (intersects.size() == 2){
			lines.push_back(intersects);
		}

		

	}

}




void SBGATFilterUQ::AddUncertaintyRegionToCovariance(int region_center_index,const double & standard_dev,const double & correl_distance){


	if (this -> P_CC.n_rows != 3 * this -> model -> GetN_vertices()){
		this -> P_CC.clear();
		this -> P_CC = arma::zeros<arma::mat>(3 * this -> model -> GetN_vertices(),3 * this -> model -> GetN_vertices());
	}
	
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

	double center[3];
	input -> GetPoint(region_center_index,center);

	for (unsigned int i = 0; i < N_C; ++i){

		double ni_[3];
		double Pi_[3];
		input -> GetPoint(i,Pi_);
		normals -> GetTuple(i,ni_);
		arma::vec::fixed<3> ni = {ni_[0],ni_[1],ni_[2]};
		arma::vec::fixed<3> u_2 = arma::normalise(arma::cross(ni,arma::randn<arma::vec>(3)));
		arma::vec u_1 = arma::cross(u_2,ni);
		
		double distance_from_center = this -> model -> GetScaleFactor() * std::sqrt(vtkMath::Distance2BetweenPoints(center,Pi_));
		
		// If the following is false, skip $i, it is outside of the uncertainty region
		if(distance_from_center < 3 * correl_distance){
			double decay_from_center = std::exp(- std::pow(distance_from_center / correl_distance,2)) ;

			arma::mat::fixed<3,3> P = decay_from_center * std::pow(standard_dev,2) * (ni * ni.t() + epsilon * (u_1 * u_1.t() + u_2 * u_2.t()));

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

					arma::mat::fixed<3,3> P_correlation = decay_from_center * std::pow(standard_dev,2) * decay * ni * nj.t();

					this -> P_CC.submat(3 * i, 3 * j, 3 * i + 2, 3 * j + 2) = P_correlation;
					this -> P_CC.submat(3 * j, 3 * i, 3 * j + 2, 3 * i + 2) = P_correlation.t();

				}

			}
		}

	}

}



