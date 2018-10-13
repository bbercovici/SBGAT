
/**
@file SBGATRefFrame.hpp
@class  SBGATRefFrame
@author Benjamin Bercovici
@author Jay McMahon
@date October 2018

@brief  Defines the SBGATRefFrame class
@details Defines the SBGATRefFrame class, a right-handed set of basis vector positioned
at some origin in space. The position of the origin and the orientation of the basis vectors
are specified with respect to the parent frame of the present frame.
@copyright MIT License, Benjamin Bercovici and Jay McMahon
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
	const arma::vec::fixed<3> & get_mrp_from_parent() const;


	/**
	Returns the dcm of $this with respect to its parent
	@return dcm of $this with respect to its parent
	*/
	const arma::mat::fixed<3,3> & get_dcm_from_parent() const;

	/**
	Returns the origin of $this with respect to its parent
	expressed in the parent reference frame
	@return origin of $this with respect to its parent
	*/
	const arma::vec::fixed<3> & get_origin_from_parent() const;


protected:

	std::string name;

	arma::vec::fixed<3> mrp_from_parent;
	arma::vec::fixed<3> origin_from_parent;
	arma::mat::fixed<3,3> dcm_from_parent;



};




#endif
