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

#ifndef FRAMEGRAPH_HEADER
#define FRAMEGRAPH_HEADER

#include "SBGATRefFrame.hpp"
#include <memory>
#include "Adjacency_List.hpp"

class SBGATFrameGraph {

public:

	/**
	Constructor. Creates a undirected reference frame graph
	*/
	SBGATFrameGraph();

	/**
	Creates a frame and adds it to the graph
	@param name of grame to be added. A warning is issued if that frame name was already used (no frame is added then)
	*/
	void add_frame(std::string frame_name);


	/**
	Converts the coordinates of the provided vector from frame
	$from to frame $to. For this to work, $from and $to must be
	in the SBGATFrameGraph and a path must be connecting them
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
	SBGATRefFrame * get_frame(std::string frame_name);


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
	Adjacency_List<std::shared_ptr <SBGATRefFrame> , std::pair< std::string, std::string > > adjacency_list;
	std::map< std::string , std::shared_ptr <SBGATRefFrame> > ref_names_to_ref_ptrs;


	void convert_to_parent_of_provided_child_frame(arma::vec & coords, SBGATRefFrame * ref_frame,
	        bool conserve_norm) const;
	void convert_to_child_of_provided_parent_frame(arma::vec & coords, SBGATRefFrame * ref_frame,
	        bool conserve_norm) const;


};



#endif