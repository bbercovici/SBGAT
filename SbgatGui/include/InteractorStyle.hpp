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
#include <vtkHardwareSelector.h>
#include <vtkExtractSelectedPolyDataIds.h>
#include <vtkSelectionNode.h>
#include <vtkInformation.h>
#include <vtkGeometryFilter.h>
#include <vtkFeatureEdges.h>

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
	Selects facets by grabbing the connected vertices
	inside a rectangular selection
	*/
	void grab_area() ;


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