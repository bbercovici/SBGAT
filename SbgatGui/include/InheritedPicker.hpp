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