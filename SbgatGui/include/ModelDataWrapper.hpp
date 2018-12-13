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




#ifndef HEADER_MODELDATAWRAPPER
#define HEADER_MODELDATAWRAPPER

#include <vtkLight.h>
#include <vtkLightActor.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkFloatArray.h>
#include <vtkScalarBarActor.h>


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
	Setter of polydata.
	@param polydata pointer to polydata to assign.
	*/
		void set_polydata(vtkSmartPointer<vtkPolyData> polydata);

		/**
	Setter of scale factor 
	@param scale_factor scale factor applied to the shape model coordinates to have them scaled to meters
		*/	

		void set_scale_factor(double scale_factor);

		/**
	Getter to scale factor 
	@return  scale factor applied to the shape model coordinates to have them scaled to meters
		*/	

		double get_scale_factor() const;



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


	/**
	Setter to light
	*/
		void set_light(vtkSmartPointer<vtkLight> light);

	/**
	Getter to light
	*/
		vtkSmartPointer<vtkLight> get_light() const;


	/**
	Getter to light actor
	*/
		vtkSmartPointer<vtkLightActor> get_light_actor() const;

	/**
	Set light actor
	*/
		void set_light_actor(vtkSmartPointer<vtkLightActor> light_actor);


	/**
	Set surface slopes
	*/
		void set_slopes(std::vector<double> slopes);


	/**
	Set inertial surface gravity potentials
	*/
		void set_inertial_potentials(std::vector<double> potentials);


	/**
	Set body fixed surface gravity potentials
	*/
		void set_body_fixed_potentials(std::vector<double> potentials);


	/**
	Set surface inertial acceleration magnitudes
	*/
		void set_inertial_acc_magnitudes(std::vector<double> acc_magnitudes);


	/**
	Set surface body-fixed acceleration magnitudes
	*/
		void set_body_fixed_acc_magnitudes(std::vector<double> acc_body_fixed_magnitudes);





	/**
	Set surface slopes
	*/
		vtkSmartPointer<vtkFloatArray> get_slopes();


	/**
	Set inertial surface gravity potentials
	*/
		vtkSmartPointer<vtkFloatArray> get_inertial_potentials();


	/**
	Set body-fixed surface gravity potentials
	*/
		vtkSmartPointer<vtkFloatArray> get_body_fixed_potentials();


	/**
	Set surface inertial acceleration magnitudes
	*/
		vtkSmartPointer<vtkFloatArray> get_inertial_acc_magnitudes();


	/**
	Set surface body-fixed acceleration magnitudes
	*/
		vtkSmartPointer<vtkFloatArray> get_body_fixed_acc_magnitudes();



	/**
	Get colorbar actor
	*/
		vtkSmartPointer<vtkScalarBarActor> get_colorbar_actor();


	/**
	Set colorbar actor 
	*/
		void set_colorbar_actor(vtkSmartPointer<vtkScalarBarActor> actor);


	protected:

		vtkSmartPointer<vtkPolyData>  polydata;
		vtkSmartPointer<vtkPolyDataMapper>  mapper;
		vtkSmartPointer<vtkActor>  actor;
		vtkSmartPointer<vtkPoints>  points;
		vtkSmartPointer<vtkLight> light;
		vtkSmartPointer<vtkLightActor> light_actor;
		vtkSmartPointer<vtkScalarBarActor> colorbar_actor;


		vtkSmartPointer<vtkFloatArray> slopes = nullptr;
		vtkSmartPointer<vtkFloatArray> inertial_potentials = nullptr;
		vtkSmartPointer<vtkFloatArray> body_fixed_potentials = nullptr;

		vtkSmartPointer<vtkFloatArray> inertial_acc_magnitudes = nullptr;
		vtkSmartPointer<vtkFloatArray> body_fixed_acc_magnitudes = nullptr;


		bool consistent_shape_model = false;
		bool consistent_global_acceleration = false;
		bool consistent_global_potential = false;
		bool consistent_grav_slope = false;

		double scale_factor;



	};




}



#endif