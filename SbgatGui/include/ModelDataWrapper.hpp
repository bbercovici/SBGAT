


#ifndef HEADER_MODELDATAWRAPPER
#define HEADER_MODELDATAWRAPPER

#include <ShapeModel.hpp>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>



namespace SBGAT_GUI {

/*!
@class ModelDataWrapper
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
	Accessor to points.
	@return pointer to points.
	*/
	vtkSmartPointer<vtkPoints> get_points() const;


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
	Setter of points.
	@param points pointer to vtkPoints to assign.
	*/
	void set_points(vtkSmartPointer<vtkPoints> points);


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


	/**
	Accessor to consistent_global_acceleration.
	@return true if global pgm accelerations have been computed
	*/
	bool get_global_pgm_acc() const;


	/**
	Accessor to consistent_global_potential
	@return true if global pgm potentials have been computed
	*/
	bool get_global_pgm_pot() const;

	/**
	Accessor to consistent_grav_slope.
	@return true gravity slopes have been computed
	*/
	bool get_grav_slopes() const;

	/**
	Setter to consistent_global_acceleration
	@param global_pgm true if consistent global accelerations have been
	computed
	*/

	void set_global_pgm_acc(bool global_pgm_acc);


	/**
	Setter to consistent_global_potential
	@param global_pgm_pot true if consistent global potentials have been
	computed
	*/
	void set_global_pgm_pot(bool global_pgm_pot);


	/**
	Setter to consistent_global_acceleration
	@param grav_slopes true if consistent global accelerations have been
	computed
	*/

	void set_grav_slopes(bool grav_slopes);

	/**
	Sets the internal flag consistent_shape_model to false, 
	indicating that the VTK shape and the SbgatCore shape are different
	*/
	void shape_was_modified();

	/**
	Returns the internal flag consistent_shape_model
	@return true if the VTK shape and the SbgatCore are identical, false otherwise
	*/
	bool get_consistent_shape_model() const;



protected:

	std::shared_ptr<SBGAT_CORE::ShapeModel>  shape_model;
	vtkSmartPointer<vtkPolyData>  polydata;
	vtkSmartPointer<vtkPolyDataMapper>  mapper;
	vtkSmartPointer<vtkActor>  actor;
	vtkSmartPointer<vtkPoints>  points;


	bool consistent_shape_model = false;
	bool consistent_global_acceleration = false;
	bool consistent_global_potential = false;
	bool consistent_grav_slope = false;



};




}



#endif