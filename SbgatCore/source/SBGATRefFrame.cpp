/** MIT License

Copyright (c) 2018 Benjamin Bercovici and Jay McMahon

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

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
