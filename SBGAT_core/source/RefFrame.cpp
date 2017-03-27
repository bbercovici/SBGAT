#include "RefFrame.hpp"

RefFrame::RefFrame(std::string name) {

	this -> name = name;

	arma::vec mrp = {0, 0, 0};
	arma::vec origin = {0, 0, 0};
	arma::mat dcm = arma::eye<arma::mat>(3, 3);


	this -> mrp_from_parent = std::make_shared<arma::vec> (mrp);
	this -> origin_from_parent = std::make_shared<arma::vec> (origin);
	this -> dcm_from_parent = std::make_shared<arma::mat>(dcm);
}


RefFrame& RefFrame::operator=(const RefFrame & other) {
	this -> name = other.name;
	*this -> mrp_from_parent = *(other.mrp_from_parent);
	*this -> origin_from_parent = *(other.origin_from_parent);
	*this -> dcm_from_parent = *(other . dcm_from_parent);


	return *this;
}

RefFrame::RefFrame( const RefFrame &ref_frame) {
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


arma::vec * RefFrame::get_mrp_from_parent() {
	return this -> mrp_from_parent.get();
}


arma::mat * RefFrame::get_dcm_from_parent() {
	return this -> dcm_from_parent.get();
}



arma::vec * RefFrame::get_origin_from_parent() {
	return this -> origin_from_parent.get();
}


void RefFrame::set_mrp_from_parent(arma::vec & mrp) {
	*this -> mrp_from_parent = mrp;
	*this -> dcm_from_parent = mrp_to_dcm(mrp);
}

void RefFrame::set_origin_from_parent(arma::vec & origin) {
	*this -> origin_from_parent = origin;
}


std::string RefFrame::get_name() const {
	return this -> name;
};




RefFrame::~RefFrame() {

}
