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


void ModelDataWrapper::set_shape_model(std::shared_ptr<SBGAT_CORE::ShapeModel> shape_model) {
	this -> shape_model = shape_model;
}


void ModelDataWrapper::set_polydata(vtkSmartPointer<vtkPolyData> polydata) {
	this -> polydata = polydata;
}


void ModelDataWrapper::set_mapper(vtkSmartPointer<vtkPolyDataMapper> mapper) {
	this -> mapper = mapper;
}


void ModelDataWrapper::set_actor(vtkSmartPointer<vtkActor> actor) {
	this -> actor = actor;
}



