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


#ifndef HEADER_VertexEditionWindow
#define HEADER_VertexEditionWindow

#include <QMainWindow>
#include <QComboBox>
#include <QVTKOpenGLWidget.h>
#include <QDialog.h>
#include <QSlider>
#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <QSpinBox.h>

#include "Mainwindow.hpp"

namespace SBGAT_GUI {

	class Mainwindow;

/*!
@class VertexEditionWindow
\author Benjamin Bercovici
\date February 2019
\brief VertexEditionWindow
\details A window where the user can move the selected vertex or vertex neighborhood
*/

	class VertexEditionWindow : public QDialog {
		Q_OBJECT

	public:

	/**
	Creates the window
	@param parent pointer to parent window.
	*/
		VertexEditionWindow(Mainwindow * parent) ;

		private slots:
		void move_vertex();
		void update_direction(int index);
		void close();
		void accept();
		void update_neighborhood_size(int neighborhood_size);
		void update_neighbors_positions();


	protected:

		vtkSmartPointer<vtkActor> line_actor;
		vtkSmartPointer<vtkDataSetMapper> line_mapper;
		std::vector< vtkSmartPointer<vtkActor> > neighbors_actors_vector;
		std::vector< double > neighbors_distances_vector;


		std::vector< vtkIdType > neighbors_indices_vector;

		QSlider * position_slider;
		QComboBox * direction_combo_box;
		QSpinBox * neighborhood_spin_box;
		Mainwindow * parent;
		double dir[3];
		double queried_point[3];
		void init();



	};
}
#endif