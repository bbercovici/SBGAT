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


#ifndef HEADER_REFFRAME
#define HEADER_REFFRAME

#include <memory>
#include <armadillo>
#include <RigidBodyKinematics.hpp>






class SBGATRefFrame {

public:

	/**
	Constructor
	@param name Reference frame name
	*/
	SBGATRefFrame(std::string name);


	/**
	Copy constructor (RO3)
	@param ref_frame Reference frame
	*/
	SBGATRefFrame( const SBGATRefFrame &ref_frame);

	/**
	Assignment operator (RO3)
	*/
	SBGATRefFrame& operator=(const SBGATRefFrame & other);

	/**
	Destructor (RO3)
	*/
	~SBGATRefFrame();

	/**
	Returns the reference frame name
	@return name of reference frame
	*/
	std::string get_name() const;

	/**
	Sets the mrp of $this relative to its parent
	to the prescribed value
	@param mrp MRP
	*/
	void set_mrp_from_parent(arma::vec & mrp) ;


	/**
	Sets the origin of $this relative to its parent
	to the prescribed value
	@param origin Origin of $this with respect to its
	parent expressed in the parent frame
	*/
	void set_origin_from_parent(arma::vec & origin) ;

	/**
	Returns the mrp of $this with respect to its parent
	@return mrp of $this with respect to its parent
	*/
	arma::vec * get_mrp_from_parent();


	/**
	Returns the dcm of $this with respect to its parent
	@return dcm of $this with respect to its parent
	*/
	arma::mat * get_dcm_from_parent();

	/**
	Returns the origin of $this with respect to its parent
	expressed in the parent reference frame
	@return origin of $this with respect to its parent
	*/
	arma::vec * get_origin_from_parent();


protected:

	std::string name;

	std::shared_ptr<arma::vec> mrp_from_parent;
	std::shared_ptr<arma::vec> origin_from_parent;
	std::shared_ptr<arma::mat> dcm_from_parent;



};




#endif
