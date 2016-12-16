
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
#include "Mainwindow.hpp"
#include "InheritedPicker.hpp"

// forward declaration of the Mainwindow class
class Mainwindow;

/**
Declaration of the InteractorStyle. This class enables the user to access (read AND write) the underlying data
displayed in the main window's QVTKWidget
*/

class InteractorStyle : public vtkInteractorStyleRubberBandPick {
public:
	static InteractorStyle * New();
	vtkTypeMacro(InteractorStyle, vtkInteractorStyleRubberBandPick);
	
	InteractorStyle();
	virtual void OnLeftButtonUp();

	/**
	Enables the interactor to get access to the GUI's mainwindow
	@param points_polydata Pointer to the GUI's mainwindow
	*/
	void set_mainwindow(Mainwindow * mainwindow);

	/**
	Returns a pointer to the GUI's mainwindow (i.e the highest-level QT Widget)
	*/
	Mainwindow * get_mainwindow();

	/**
	Sets the interactor style
	@param mode Int representing the mode (0: orient, 1:select)
	*/
	void set_current_mode(const int mode);

	vtkSmartPointer<vtkPolyData> get_selected_points_polydata();


private:
	vtkSmartPointer<vtkPolyData> selected_points_polydata;
	Mainwindow * mainwindow ;


};

#endif