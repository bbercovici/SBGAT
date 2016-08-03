
#ifndef HEADER_INTERACTORSTYLE
#define HEADER_INTERACTORSTYLE

#include <vtkSmartPointer.h>
#include <vtkInteractorStyleRubberBandPick.h>
#include <vtkPolyData.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkIdTypeArray.h>
#include <vtkRendererCollection.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkAreaPicker.h>
#include <vtkExtractGeometry.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkIdFilter.h>
#include <vtkRenderWindow.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkExtractGeometry.h>
#include <vtkPlanes.h>
#include <vtkHardwareSelector.h>
#include <vtkSelection.h>
#include <vtkExtractSelection.h>
#include <vtkSelectVisiblePoints.h>
#include <vtkPointSource.h>

#include "vtkObjectFactory.h"
#include "mainwindow.h"

// forward declaration of the MainWindow class
class MainWindow;

/**
Declaration of the InteractorStyle. This class enables the user to access (read AND write) the underlying data
displayed in the main window's QVTKWidget
*/

class InteractorStyle : public vtkInteractorStyleRubberBandPick {
public:
	static InteractorStyle * New();
	vtkTypeMacro(InteractorStyle, vtkInteractorStyleRubberBandPick);
	InteractorStyle();
	virtual ~InteractorStyle();
	virtual void OnLeftButtonUp();

	/**
	Enables the interactor to get access to the vtkPolyData storing the shape's vertices
	@param points_polydata Pointer to the vtkPolyData storing the vertices of the shape model
	*/
	void SetPoints(vtkSmartPointer<vtkPolyData> points_polydata);

	/**
	Enables the interactor to get access to the GUI's mainwindow
	@param points_polydata Pointer to the GUI's mainwindow
	*/
	void SetMainWindow(MainWindow * mainwindow);

	/**
	Sets the interactor style
	@param mode Int representing the mode (0: orient, 1:select)
	*/
	void SetCurrentMode(const int mode);


	void transform_points();


private:
	vtkSmartPointer<vtkPolyData> points_polydata;
	vtkSmartPointer<vtkPolyData> selected_points_polydata;
	vtkSmartPointer<vtkActor> SelectedActor;
	vtkSmartPointer<vtkDataSetMapper> SelectedMapper;
	MainWindow * mainwindow ;
	bool * widget_is_open;


};

#endif