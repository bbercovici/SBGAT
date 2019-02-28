#include <SBGATTransformShape.hpp>
#include <SBGATMassProperties.hpp>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <RigidBodyKinematics.hpp>


void SBGATTransformShape::ShiftShapeToBarycenter(vtkSmartPointer<vtkPolyData> & shape){

	vtkSmartPointer<SBGATMassProperties> center_of_mass_filter =
	vtkSmartPointer<SBGATMassProperties>::New();

	center_of_mass_filter -> SetInputData(shape);
	center_of_mass_filter -> Update();

	arma::vec com_translation = - center_of_mass_filter -> GetCenterOfMass();
	SBGATTransformShape::Translate(com_translation,shape);

}

void SBGATTransformShape::Translate(const arma::vec::fixed<3> x,vtkSmartPointer<vtkPolyData> & shape){


	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();

	transform -> Translate( x(0),  x(1),  x(2));

	vtkSmartPointer<vtkTransformPolyDataFilter> filter =
	vtkSmartPointer<vtkTransformPolyDataFilter>::New();

	filter -> SetInputData(shape);
	filter -> SetTransform(transform);
	filter -> Update();

	shape = filter -> GetOutput();


}


void SBGATTransformShape::Rotate(const arma::mat::fixed<3,3> dcm,vtkSmartPointer<vtkPolyData> & shape){

	auto prv = RBK::dcm_to_prv(dcm);
	double angle = 180 / arma::datum::pi * arma::norm(prv);
	arma::vec axis = arma::normalise(prv);

	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	transform -> RotateWXYZ(angle,axis.colptr(0));

	vtkSmartPointer<vtkTransformPolyDataFilter> filter =
	vtkSmartPointer<vtkTransformPolyDataFilter>::New();

	filter -> SetInputData(shape);
	filter -> SetTransform(transform);
	filter -> Update();

	shape = filter -> GetOutput();



}
