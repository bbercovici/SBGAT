#ifndef HEADER_INHERITEDPICKER
#define HEADER_INHERITEDPICKER


/**
Definition of the InheritedPicker class. This class directly inherits from vtkAreaPicker.
Used as an interface that allows the user to get access to the dimensions of the picking area.
*/
class InheritedPicker : public vtkAreaPicker {
public:
	InheritedPicker() {};
	int get_X0() {
		return this -> X0;
	}
	int get_X1() {
		return this -> X1;
	}
	int get_Y0() {
		return this -> Y0;
	}
	int get_Y1() {
		return this -> Y1;
	}

	std::vector<unsigned int> get_dimensions_uint() {

		// WARNING: check order of dimensions
		std::vector<unsigned int> dimensions;
		dimensions.push_back(this -> get_X0());
		dimensions.push_back(this -> get_Y0());
		dimensions.push_back(this -> get_X1());
		dimensions.push_back(this -> get_Y1());
		return dimensions;
	}

	std::vector<int> get_dimensions_int() {
		
		// WARNING: check order of dimensions
		std::vector<int> dimensions;
		dimensions.push_back(this -> get_X0());
		dimensions.push_back(this -> get_X1());
		dimensions.push_back(this -> get_Y0());
		dimensions.push_back(this -> get_Y1());
		return dimensions;
	}

};
#endif