#include "ModelDataWrapper.hpp"

using namespace SBGAT_GUI;


ModelDataWrapper::ModelDataWrapper() {

}

std::shared_ptr<SBGAT_CORE::ShapeModel> ModelDataWrapper::get_shape_model() const {
	return this -> shape_model;
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


void ModelDataWrapper::set_shape_model(std::shared_ptr<SBGAT_CORE::ShapeModel> shape_model) {
	this -> shape_model = shape_model;
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

