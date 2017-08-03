#ifndef FRAMEGRAPH_HEADER
#define FRAMEGRAPH_HEADER

#include "RefFrame.hpp"
#include <memory>
#include "Adjacency_List.hpp"

namespace SBGAT_CORE {

class FrameGraph {

public:

	/**
	Constructor. Creates a undirected reference frame graph
	*/
	FrameGraph();

	/**
	Creates a frame and adds it to the graph
	@param name of grame to be added. A warning is issued if that frame name was already used (no frame is added then)
	*/
	void add_frame(std::string frame_name);


	/**
	Converts the coordinates of the provided vector from frame
	$from to frame $to. For this to work, $from and $to must be
	in the FrameGraph and a path must be connecting them
	@param input vector to convert
	@param from name of reference frame to convert from
	@param to name of reference frame to convert to
	@param conserve_norm defines whether the provided coordinates are that of a vector whose norm should be conserved or not

	True if the expected conversion is of the form x_B = [BN]x_N where
	- x_N are the provided coordinates in the departure frame N, 
	- x_N are the expected coordinates in the departure frame B
	- N the departure frame, 
	- B the arrival frame
	- [BN] the direction cosine matrix orienting the two frames.

	If False, then 
	a translational part T is added as in x_B = [BN]x_N + T, where 
	-  x_B are the coordinates of point x expressed in the B frame
	-  x_N are the coordinates of point x expressed in the N frame
	- [BN] being the direction cosine matrix orienting the two frames. 
	- T the coordinates of the displacement vector from the origin of frame B to that of frame N, expressed 
	in the B frame
	@return converted coordinates
	*/
	arma::vec convert(arma::vec & input, std::string from,
	                  std::string to, bool conserve_norm = false);


	/**
	Sets the mrp of the transform to the one provided as argument
	@param parent_name Name of parent frame
	@param child_name Name of chuld frame
	@param mrp MRP set
	@param check True if consistency check should
	*/
	void set_transform_mrp(std::string parent_name,
	                       std::string child_name,
	                       arma::vec & mrp);


	/**
	Sets the origin of the transform to the one provided as argument
	@param parent_name Name of parent frame
	@param child_name Name of chuld frame
	@param origin Origin of child frame expressed in parent frame
	@param check True if consistency check should
	*/
	void set_transform_origin(std::string parent_name,
	                          std::string child_name,
	                          arma::vec & origin);


	/**
	Returns a pointer to the reference frame whose name is passed as argument
	WARNINGS:
		- an exception will be thrown if the frame is not present in the graph
		- this method should only be used to read from the reference frame
		and not to modify if
	@param frame_name Name of reference frame
	@param Pointer to reference frame
	*/
	RefFrame * get_frame(std::string frame_name);


	/**
	Creates a reference frame transform between a
	parent frame $parent_name and and child frame $child_name.

	If either of the input names are not present in the graph,
	a warning is issued and nothing happens.

	If the provided transform or its opposite was already present,
	a warning is issued and nothing happens

	The transform is stored internaly  in a pair formed with the following members:
	- first: parent frame
	- second: child frame

	@param child_name Name of the children frame
	@param parent_name Name of the parent frame
	*/
	void add_transform(std::string parent_name, std::string child_name) ;


protected:
	Adjacency_List<std::shared_ptr <RefFrame> , std::pair< std::string, std::string > > adjacency_list;
	std::map< std::string , std::shared_ptr <RefFrame> > ref_names_to_ref_ptrs;


	void convert_to_parent_of_provided_child_frame(arma::vec & coords, RefFrame * ref_frame,
	        bool conserve_norm) const;
	void convert_to_child_of_provided_parent_frame(arma::vec & coords, RefFrame * ref_frame,
	        bool conserve_norm) const;


};

}

#endif