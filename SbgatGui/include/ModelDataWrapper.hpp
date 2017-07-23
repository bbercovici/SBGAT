/**
ModelDataWrapper.hpp
\author Benjamin Bercovici
\date July 22, 2017
\brief ModelDataWrapper class. Convenient wrapper around SbgatCore shape models, VTK entities and 
other shape-related that must remain together

\details This class stores pointers to various elements belonging to the same "shape model" and
that should not be separated. This includes the SbgatCore::ShapeModel object itself, its associated vtkPolydata,
vtkActor, vtkMapper, and also several boolean variables that are used to keep the GUI consistent at all time.
Using this wrapping structure thus helps condensing the code by providing unified access to all these
elements through a unique interface.
*/



#ifndef HEADER_MODELDATAWRAPPER
#define HEADER_MODELDATAWRAPPER

#include <ShapeModel.hpp>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>



namespace SBGAT_GUI {


class ModelDataWrapper {

public:
	ModelDataWrapper();

	/**
	Accessor to shape model.
	@return pointer to shape model.
	*/
	std::shared_ptr<SBGAT_CORE::ShapeModel> get_shape_model() const;

	/**
	Accessor to polydata.
	@return pointer to polydata.
	*/
	vtkSmartPointer<vtkPolyData> get_polydata() const;

	/**
	Accessor to mapper.
	@return pointer to mapper.
	*/
	vtkSmartPointer<vtkPolyDataMapper> get_mapper() const;


	/**
	Accessor to actor.
	@return pointer to actor.
	*/
	vtkSmartPointer<vtkActor> get_actor() const;


	/**
	Setter of shape model.
	@param shape_model pointer to shape model to assign.
	*/
	void set_shape_model(std::shared_ptr<SBGAT_CORE::ShapeModel> shape_model);


	/**
	Setter of polydata.
	@param polydata pointer to polydata to assign.
	*/
	void set_polydata(vtkSmartPointer<vtkPolyData> polydata);


	/**
	Setter of mapper.
	@param mapper pointer to mapper to assign.
	*/
	void set_mapper(vtkSmartPointer<vtkPolyDataMapper> mapper);


	/**
	Setter of actor.
	@param actor pointer to actor to assign.
	*/
	void set_actor(vtkSmartPointer<vtkActor> actor);


protected:

	std::shared_ptr<SBGAT_CORE::ShapeModel>  shape_model;
	vtkSmartPointer<vtkPolyData>  polydata;
	vtkSmartPointer<vtkPolyDataMapper>  mapper;
	vtkSmartPointer<vtkActor>  actor;

	bool consistent_global_acceleration;
	bool consistent_grav_slope;


};




}



#endif