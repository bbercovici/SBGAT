
#ifndef HEADER_REFFRAME
#define HEADER_REFFRAME

#include <memory>
#include <armadillo>
#include "RigidBodyKinematics.hpp"
class RefFrame {

public:

	/**
	Constructor
	@param name Reference frame name
	*/
	RefFrame(std::string name);


	/**
	Copy constructor (RO3)
	@param ref_frame Reference frame
	*/
	RefFrame( const RefFrame &ref_frame);

	/**
	Assignment operator (RO3)
	*/
	RefFrame& operator=(const RefFrame & other);

	/**
	Destructor (RO3)
	*/
	~RefFrame();

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
