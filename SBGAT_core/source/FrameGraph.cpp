#include "FrameGraph.hpp"

FrameGraph::FrameGraph() {

}



arma::vec FrameGraph::convert(arma::vec & input, std::string from, std::string to,
                              bool is_unit_vector) {

	std::deque<std::shared_ptr<RefFrame > > path = this -> adjacency_list.dfs(
	            ref_names_to_ref_ptrs[from], ref_names_to_ref_ptrs[to]);

	arma::vec coords = input;
	if (from == to) {
		return coords;
	}

	for (auto it_current_frame = path.begin();
	        it_current_frame != --path.end();
	        ++it_current_frame) {

		auto it_next_frame = std::next(it_current_frame);

		std::pair<std::string, std::string> transform = this -> adjacency_list.getedge(*it_current_frame, *it_next_frame);

		// current frame is parent frame
		if (transform.second == (*it_next_frame) -> get_name()) {
			this -> convert_to_child_of_provided_parent_frame(coords, (*it_next_frame).get(),
			        is_unit_vector);
		}

		// current frame is child frame
		else if (transform.first == (*it_next_frame) -> get_name()) {

			this -> convert_to_parent_of_provided_child_frame(coords, (*it_current_frame).get(),
			        is_unit_vector);
		}
		else {
			throw (std::runtime_error("Illegal frame conversion"));
		}


	}

	return coords;
}


void FrameGraph::convert_to_parent_of_provided_child_frame(arma::vec & coords,
        RefFrame * ref_frame, bool is_unit_vector) const {

	if (is_unit_vector == false) {
		coords = *ref_frame -> get_origin_from_parent() +  (*ref_frame -> get_dcm_from_parent()).t() * coords;
	}
	else {
		coords = (*ref_frame -> get_dcm_from_parent()).t() * coords;

	}
}

void FrameGraph::convert_to_child_of_provided_parent_frame(arma::vec & coords,
        RefFrame * ref_frame, bool is_unit_vector) const {

	if (is_unit_vector == false) {
		coords = (*ref_frame -> get_dcm_from_parent()) * ( coords - (*ref_frame -> get_origin_from_parent()) );
	}
	else {
		coords = (*ref_frame -> get_dcm_from_parent()) * ( coords );
	}
}



void FrameGraph::add_frame(std::string frame_name) {

	std::shared_ptr<RefFrame> frame = std::make_shared<RefFrame>(RefFrame(frame_name));

	std::set<std::shared_ptr<RefFrame> > frames = this -> adjacency_list.get_vertices();

	for (auto const & existing_frame : frames) {
		if (existing_frame -> get_name() == frame_name) {
			std::cerr << "The reference frame name ' " << frame_name << " 'is already in use" << std::endl;
			return;
		}
	}

	this -> adjacency_list.addvertex(frame);
	this -> ref_names_to_ref_ptrs[frame_name] = frame;
}


void FrameGraph::set_transform_mrp(std::string parent_name,
                                   std::string child_name,
                                   arma::vec & mrp) {


	std::set<std::pair<std::string, std::string > > transforms = this -> adjacency_list.get_edges();

	for (auto const & transform : transforms) {

		if (transform.first == parent_name) {

			if (transform.second == child_name) {
				this -> ref_names_to_ref_ptrs[child_name]-> set_mrp_from_parent(mrp) ;
				return;
			}
		}

	}

	std::cerr << "The transform relating parent '" << parent_name << "' to child '" <<  child_name << "' was not found in the graph" << std::endl;
	return;
}


void FrameGraph::set_transform_origin(std::string parent_name,
                                      std::string child_name,
                                      arma::vec & origin) {

	// Consistency check: if $child_name == "N", something is wrong
	if (child_name == "N") {
		std::cerr << " The inertial frame of reference N cannot have its origin edited" << std::endl;
		return;
	}

	std::set<std::pair<std::string, std::string > > transforms = this -> adjacency_list.get_edges();

	for (auto const & transform : transforms) {

		if (transform.first == parent_name) {

			if (transform.second == child_name) {
				this -> ref_names_to_ref_ptrs[child_name]-> set_origin_from_parent(origin) ;
				return;
			}
		}

	}

	std::cerr << "The transform relating parent '" << parent_name << "' to child '" <<  child_name << "' was not found in the graph" << std::endl;
	return;
}

RefFrame * FrameGraph::get_frame(std::string frame_name) {
	return this -> ref_names_to_ref_ptrs[frame_name].get();
}

void FrameGraph::add_transform(std::string parent_name, std::string child_name) {

	std::set<std::shared_ptr<RefFrame> > frames = this -> adjacency_list.get_vertices();

	std::shared_ptr<RefFrame> child_frame;
	std::shared_ptr<RefFrame> parent_frame;

	//########################################################################
	//####################### Consistency checks #############################
	//########################################################################

	// Consistency check: are the two provided frames identically called?
	if (parent_name == child_name) {
		std::cerr << "The parent name has to be different from the child name. Both were called " << child_name << std::endl;
		return;
	}


	// Consistency check: are both frames present?
	bool found_parent = false;
	for (auto const & existing_frame : frames) {
		if (existing_frame -> get_name() == parent_name) {
			found_parent = true;
			parent_frame = existing_frame;
		}
	}

	if (found_parent == false) {
		std::cerr << "The parent reference frame name '" << parent_name << "' was not found in the graph" << std::endl;
		return;
	}

	bool found_child = false;
	for (auto const & existing_frame : frames) {
		if (existing_frame -> get_name() == child_name) {
			found_child = true;
			child_frame = existing_frame;

		}
	}


	if (found_child == false) {
		std::cerr << "The child reference frame name '" << child_name << "' was not found in the graph" << std::endl;
		return;
	}

	// Consistency check: is this transform already present?
	std::set<std::pair<std::string, std::string > > transforms = this -> adjacency_list.get_edges();

	for (auto const & transform : transforms) {
		if (transform.first == parent_name || transform.first == child_name) {
			if (transform.second == parent_name || transform.second == child_name) {
				std::cerr << "A transform relating '" << parent_name << "' and '" <<  child_name << "' was found in the graph" << std::endl;
				return;
			}
		}
	}


	//########################################################################
	//########################################################################
	//########################################################################


	// The edge is then created by storing the name of the corresponding tansform
	// The transform is stored in a pair formed with the following members:
	// - first: parent frame
	// - second: child frame
	std::pair<std::string, std::string> transform_name = std::make_pair(parent_name, child_name);
	this -> adjacency_list.addedge(parent_frame, child_frame, transform_name);


}