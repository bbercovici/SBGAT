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

#include "ModelDataWrapper.hpp"

using namespace SBGAT_GUI;


ModelDataWrapper::ModelDataWrapper() {

}


vtkSmartPointer<vtkPolyData> ModelDataWrapper::get_polydata() const {
	return this -> polydata;
}


vtkSmartPointer<vtkPolyDataMapper> ModelDataWrapper::get_mapper() const {
	return this -> mapper;
}


vtkSmartPointer<vtkActor> ModelDataWrapper::get_actor() const {
	return this -> actor;
}


vtkSmartPointer<vtkPoints> ModelDataWrapper::get_points() const {
	return this -> points;
}

void ModelDataWrapper::set_light(vtkSmartPointer<vtkLight> light){
	this -> light = light;
}

vtkSmartPointer<vtkLight> ModelDataWrapper::get_light() const{
	return this -> light;
}



void ModelDataWrapper::set_polydata(vtkSmartPointer<vtkPolyData> polydata) {
	this -> polydata = polydata;
}


void ModelDataWrapper::set_points(vtkSmartPointer<vtkPoints> points){
	this -> points = points;
}


void ModelDataWrapper::set_mapper(vtkSmartPointer<vtkPolyDataMapper> mapper) {
	this -> mapper = mapper;
}


void ModelDataWrapper::set_actor(vtkSmartPointer<vtkActor> actor) {
	this -> actor = actor;
}



bool ModelDataWrapper::get_global_pgm_acc() const {
	return this -> consistent_global_acceleration;
}


bool ModelDataWrapper::get_global_pgm_pot() const {
	return this -> consistent_global_potential;
}


bool ModelDataWrapper::get_grav_slopes() const {
	return this -> consistent_grav_slope;
}

void ModelDataWrapper::set_global_pgm_acc(bool global_pgm_acc) {
	this -> consistent_global_acceleration = global_pgm_acc;
}


void ModelDataWrapper::set_global_pgm_pot(bool global_pgm_pot) {
	this -> consistent_global_potential = global_pgm_pot;
}


void ModelDataWrapper::set_grav_slopes(bool grav_slopes) {
	this -> consistent_grav_slope = grav_slopes;
}

void ModelDataWrapper::shape_was_modified(){
	this -> consistent_shape_model = false;
}

bool ModelDataWrapper::get_consistent_shape_model() const {
	return this -> consistent_shape_model;
}

void ModelDataWrapper::set_light_actor(vtkSmartPointer<vtkLightActor> light_actor){
	this -> light_actor = light_actor;
}

vtkSmartPointer<vtkLightActor> ModelDataWrapper::get_light_actor() const{
	return this -> light_actor;
}


void ModelDataWrapper::set_slopes(std::vector<double> slopes){

	// Create cell data
	if (this -> slopes == nullptr){

		this -> slopes = vtkSmartPointer<vtkFloatArray>::New();
		
		for (int i = 0; i < this -> polydata -> GetNumberOfCells(); i++){
			this -> slopes -> InsertNextValue(slopes[i]);
		}

	}

	else{
		for (int i = 0; i < this -> polydata -> GetNumberOfCells(); i++){
			this -> slopes -> SetValue(i,slopes[i]);
			
		}
	}

	this -> slopes -> Modified();

	
}


void ModelDataWrapper::set_potentials(std::vector<double> potentials){

	
	// Create cell data
	if (this -> potentials == nullptr){
		this -> potentials = vtkSmartPointer<vtkFloatArray>::New();
		
		for (int i = 0; i < this -> polydata -> GetNumberOfCells(); i++){
			this -> potentials -> InsertNextValue(potentials[i]);
		}

	}

	else{

		for (int i = 0; i < this -> polydata -> GetNumberOfCells(); i++){
			this -> potentials -> SetValue(i,potentials[i]);
		}

	}
	this -> potentials -> Modified();


	
}


void ModelDataWrapper::set_acc_magnitudes(std::vector<double> acc_magnitudes){

	// Create cell data
	if (this -> acc_magnitudes == nullptr){
		this -> acc_magnitudes = vtkSmartPointer<vtkFloatArray>::New();
		
		for (int i = 0; i < this -> polydata -> GetNumberOfCells(); i++){
			this -> acc_magnitudes -> InsertNextValue(acc_magnitudes[i]);
		}

	}

	else{

		for (int i = 0; i < this -> polydata -> GetNumberOfCells(); i++){
			this -> acc_magnitudes -> SetValue(i,acc_magnitudes[i]);
		}

	}
	this -> acc_magnitudes -> Modified();

}


void ModelDataWrapper::set_acc_body_fixed_magnitudes(std::vector<double> acc_body_fixed_magnitudes){

	// Create cell data
	if (this -> acc_body_fixed_magnitudes == nullptr){
		this -> acc_body_fixed_magnitudes = vtkSmartPointer<vtkFloatArray>::New();
		
		for (int i = 0; i < this -> polydata -> GetNumberOfCells(); i++){
			this -> acc_body_fixed_magnitudes -> InsertNextValue(acc_body_fixed_magnitudes[i]);
		}

	}

	else{

		for (int i = 0; i < this -> polydata -> GetNumberOfCells(); i++){
			this -> acc_body_fixed_magnitudes -> SetValue(i,acc_body_fixed_magnitudes[i]);
		}

	}
	this -> acc_body_fixed_magnitudes -> Modified();


}

vtkSmartPointer<vtkFloatArray> ModelDataWrapper::get_slopes(){
	return this -> slopes;
}

vtkSmartPointer<vtkFloatArray> ModelDataWrapper::get_potentials(){
	return this -> potentials;
}

vtkSmartPointer<vtkFloatArray> ModelDataWrapper::get_acc_magnitudes(){
	return this -> acc_magnitudes;
}

vtkSmartPointer<vtkFloatArray> ModelDataWrapper::get_acc_body_fixed_magnitudes(){
	return this -> acc_body_fixed_magnitudes;
}


vtkSmartPointer<vtkScalarBarActor> ModelDataWrapper::get_colorbar_actor(){
	return this -> colorbar_actor;
}

void ModelDataWrapper::set_colorbar_actor(vtkSmartPointer<vtkScalarBarActor> actor){
	this -> colorbar_actor = actor;
}

void ModelDataWrapper::set_scale_factor(double scale_factor){
	this -> scale_factor =scale_factor;
}

double ModelDataWrapper::get_scale_factor() const{
	return this -> scale_factor;
}


