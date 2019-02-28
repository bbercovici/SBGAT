/*=========================================================================
 *
 *  Copyright David Doria 2011 daviddoria@gmail.com
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

/** SBGATObjWriter - Write a Wavefront OBJ (ASCII) file. The core of this class
  * was taken from vtkOBJExporter. The goal is to make it follow the standard
  * "writer" syntax that is present for the rest of the VTK writers (vtkPLYWriter, etc).
  * It only saves the geometry (vertices, lines, polylines, strips) with normals and tcoords (if any).
  * It does not write triangle strips properly.
  */

#ifndef __SBGATObjWriter_h
#define __SBGATObjWriter_h

#include "vtkPolyDataAlgorithm.h" //superclass

class SBGATObjWriter : public vtkPolyDataAlgorithm
{
public:
  static SBGATObjWriter *New();
  vtkTypeMacro(SBGATObjWriter,vtkPolyDataAlgorithm)
  void PrintSelf(ostream& os, vtkIndent indent) override;
  void PrintHeader(std::ostream& os, vtkIndent indent) override;
  void PrintTrailer(std::ostream& os, vtkIndent indent) override;


  // Description:
  // Specify the name of the file to write out.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

protected:
  SBGATObjWriter();
  ~SBGATObjWriter();

private:
  SBGATObjWriter(const SBGATObjWriter&);  // Not implemented.
  void operator=(const SBGATObjWriter&);  // Not implemented.

  char *FileName;
};

#endif
