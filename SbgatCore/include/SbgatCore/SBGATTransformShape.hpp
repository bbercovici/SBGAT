#ifndef SBGATTRANSFORMSHAPE_HEADER
#define SBGATTRANSFORMSHAPE_HEADER

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <armadillo>

class SBGATTransformShape{


public:


	/**
	Shifts origin of provided shape's coordinates to shape's barycenter
	@param shape pointer to shape to align
	*/
	static void ShiftShapeToBarycenter(vtkSmartPointer<vtkPolyData> shape);

	/**
	Shifts origin of provided shape's coordinates to shape's barycenter
	and rotates coordinates post-translation so that the shape's coordinates are expressed 
	in the shape's barycentered principal frame
	@param shape pointer to shape to align
	*/
	static void ShiftShapeToBarycenterAndRotateToPrincipalAxes(vtkSmartPointer<vtkPolyData> shape);

	/**
	Applies translation to provided shape
	@param x translation vector
	@param shape pointer to shape to align
	*/
	static void Translate(const arma::vec::fixed<3> x,vtkSmartPointer<vtkPolyData> shape);

	/**
	Applies rotation to provided shape
	@param dcm rotation matrix. If B designates the frame in which the coordinates are expressed BEFORE
	the rotation and B' the frame in in which the coordinates are expressed AFTER
	the rotation, then dcm == [B'B]
	@param shape pointer to shape to align
	*/
	static void Rotate(const arma::mat::fixed<3,3> dcm,vtkSmartPointer<vtkPolyData> shape);









};



#endif
