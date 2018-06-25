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


#ifndef HEADER_LC_VIZUALIZER
#define HEADER_LC_VIZUALIZER


#include "LCWindow.hpp"
#include <QVTKOpenGLWidget.h>
#include <vtkRenderer.h>
#include <vtkImageData.h>
#include <QPushButton>

#include <vtkContextView.h>
#include <QVTKWidget.h>


namespace SBGAT_GUI {

	class RadarWindow;

/*!
@class LCVisualizer
\author Benjamin Bercovici
\date March, 2018
\brief LCVisualizer class defining a window where a user can visualize previously computed lightcurves
\details TODO
*/

	class LCVisualizer : public QDialog {
		Q_OBJECT

	public:

	/**
	Creates the settings window
	@param parent pointer to parent window.
	@param measurements reference to vector of arrays of (time,luminosity)
	*/
		LCVisualizer(LCWindow * parent,const std::vector<std::array<double, 2> > & measurements) ;

		/**
		Initializes the visualizer window
		*/	

		void init(const std::vector<std::array<double, 2> > & measurements);

		private slots:



	protected:


		QVTKOpenGLWidget * qvtkWidget;

		LCWindow * parent;
		vtkSmartPointer<vtkContextView>  view;

		QDialogButtonBox * button_box;

		
	};
}
#endif