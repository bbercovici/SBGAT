#include "SBGATRefFrame.hpp"




SBGATRefFrame::SBGATRefFrame(std::string name) {

	this -> name = name;

	arma::vec mrp = {0, 0, 0};
	arma::vec origin = {0, 0, 0};
	arma::mat dcm = arma::eye<arma::mat>(3, 3);


	this -> mrp_from_parent = std::make_shared<arma::vec> (mrp);
	this -> origin_from_parent = std::make_shared<arma::vec> (origin);
	this -> dcm_from_parent = std::make_shared<arma::mat>(dcm);
}


SBGATRefFrame& SBGATRefFrame::operator=(const SBGATRefFrame & other) {
	this -> name = other.name;
	*this -> mrp_from_parent = *(other.mrp_from_parent);
	*this -> origin_from_parent = *(other.origin_from_parent);
	*this -> dcm_from_parent = *(other . dcm_from_parent);


	return *this;
}

SBGATRefFrame::SBGATRefFrame( const SBGATRefFrame &ref_frame) {
	this -> name = ref_frame.get_name();

	arma::vec mrp = {0, 0, 0};
	arma::vec origin = {0, 0, 0};
	arma::mat dcm_from_parent = arma::eye<arma::mat>(3, 3);


	this -> mrp_from_parent = std::make_shared<arma::vec> (mrp);
	this -> origin_from_parent = std::make_shared<arma::vec> (origin);
	this -> dcm_from_parent = std::make_shared<arma::mat>(dcm_from_parent);


	*this -> mrp_from_parent = *(ref_frame . mrp_from_parent);
	*this -> origin_from_parent = *(ref_frame . origin_from_parent);
	*this -> dcm_from_parent = *(ref_frame . dcm_from_parent);
}


arma::vec * SBGATRefFrame::get_mrp_from_parent() {
	return this -> mrp_from_parent.get();
}


arma::mat * SBGATRefFrame::get_dcm_from_parent() {
	return this -> dcm_from_parent.get();
}



arma::vec * SBGATRefFrame::get_origin_from_parent() {
	return this -> origin_from_parent.get();
}


void SBGATRefFrame::set_mrp_from_parent(arma::vec & mrp) {
	*this -> mrp_from_parent = mrp;
	*this -> dcm_from_parent = RBK::mrp_to_dcm(mrp);
}

void SBGATRefFrame::set_origin_from_parent(arma::vec & origin) {
	*this -> origin_from_parent = origin;
}


std::string SBGATRefFrame::get_name() const {
	return this -> name;
};




SBGATRefFrame::~SBGATRefFrame() {

}
