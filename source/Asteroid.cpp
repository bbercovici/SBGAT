#include <iostream>
#include <cassert>
#include <cmath>
#include <set>
#include <map>

#include <armadillo>


#include "Asteroid.hpp"
#include "Matrix.hpp"
#include "Vect.hpp"

Asteroid::Asteroid(std::ifstream& GravityFile)
{

  assert(GravityFile.is_open());

  double Gs;

  GravityFile >> Gs;

  mGs = Gs;

  int NOV, NOF, NOE;

  GravityFile >> NOV >> NOF >> NOE;

  mNOV = NOV;
  mNOF = NOF;
  mNOE = NOE;

  assert(mNOF = 2 * mNOV - 4);
  assert(mNOE = 3 * (mNOV - 2));

  mX = new double [mNOV];
  mY = new double [mNOV];
  mZ = new double [mNOV];

  for (int i = 0; i < mNOV; i++)
  {
    GravityFile >> mX[i] >> mY[i] >> mZ[i];
  }


  mListTri = new int* [mNOF];
  mListN = new double* [mNOF];
  mF = new double* [mNOF];

  for (int i = 0; i < mNOF; i++)
  {
    mListTri[i] = new int [3];
    mListN[i]   = new double [3];

    mF[i] = new double [9];

    GravityFile >> mListTri[i][0] >> mListTri[i][1] >> mListTri[i][2];
    GravityFile >> mListN[i][0] >> mListN[i][1] >> mListN[i][2];
    GravityFile >> mF[i][0] >> mF[i][1] >> mF[i][2] >> mF[i][3] >> mF[i][4] >> mF[i][5] >> mF[i][6] >> mF[i][7] >> mF[i][8];
  }

  mListE = new int* [mNOE];

  mE = new double* [mNOE];

  for (int i = 0; i < mNOE; i++)
  {
    mListE[i] = new int [2];

    mE[i] = new double [9];

    GravityFile >> mListE[i][0] >> mListE[i][1];
    GravityFile >> mE[i][0] >> mE[i][1] >> mE[i][2] >> mE[i][3] >> mE[i][4] >> mE[i][5] >> mE[i][6] >> mE[i][7] >> mE[i][8];
  }
}

Asteroid::Asteroid(vtkSmartPointer<vtkPolyData> input_polydata, double Gs) {

  unsigned int mNOV = input_polydata -> GetNumberOfPoints();
  unsigned int mNOF = input_polydata -> GetNumberOfPolys();
  std::map<std::set<unsigned int>, std::pair<int, int> > edge_info;
  vtkSmartPointer<vtkIdList> facet_vertex_indices = vtkSmartPointer<vtkIdList>::New();

  assert(mNOF = 2 * mNOV - 4);

  mX = new double [mNOV];
  mY = new double [mNOV];
  mZ = new double [mNOV];

  // Vertex coordinates are fetched in
  for (vtkIdType i = 0; i < mNOV; ++i) {
    double p[3];
    input_polydata -> GetPoint(i, p);

    mX[i] = p[0];
    mY[i] = p[1];
    mZ[i] = p[2];
  }

  // Facet indices are fetched in
  mListTri = new int* [mNOF];
  mListN = new double* [mNOF];
  mF = new double* [mNOF];

  for (vtkIdType i = 0; i < mNOF; ++i) {

    mListTri[i] = new int [3];
    mListN[i]   = new double [3];

    mF[i] = new double [9];

    input_polydata -> GetCellPoints(i, facet_vertex_indices);


    // Vertex indices in facet
    unsigned int P1_index = facet_vertex_indices -> GetId(0);
    unsigned int P2_index = facet_vertex_indices -> GetId(1);
    unsigned int P3_index = facet_vertex_indices -> GetId(2);

    mListTri[i][0] = P1_index;
    mListTri[i][1] = P2_index;
    mListTri[i][2] = P3_index;


    // Facet normal coordinates
    arma::vec P1 = {
      input_polydata -> GetPoints() -> GetData() -> GetComponent(P1_index, 0),
      input_polydata -> GetPoints() -> GetData() -> GetComponent(P1_index, 1),
      input_polydata -> GetPoints() -> GetData() -> GetComponent(P1_index, 2)
    };
    arma::vec P2 = {
      input_polydata -> GetPoints() -> GetData() -> GetComponent(P2_index, 0),
      input_polydata -> GetPoints() -> GetData() -> GetComponent(P2_index, 1),
      input_polydata -> GetPoints() -> GetData() -> GetComponent(P2_index, 2)
    };
    arma::vec P3 = {
      input_polydata -> GetPoints() -> GetData() -> GetComponent(P3_index, 0),
      input_polydata -> GetPoints() -> GetData() -> GetComponent(P3_index, 1),
      input_polydata -> GetPoints() -> GetData() -> GetComponent(P3_index, 2)

    };

    arma::vec n = arma::cross(P3 - P2, P1 - P2) / (arma::norm(P3 - P2) *  arma::norm(P1 - P2));
    mListN[i][0] = n(0);
    mListN[i][1] = n(1);
    mListN[i][2] = n(2);

    // Facet Dyad
    arma::mat F = n * n.t();

    mF[i][0] = F.row(0)(0);
    mF[i][1] = F.row(0)(1);
    mF[i][2] = F.row(0)(2);
    mF[i][3] = F.row(1)(0);
    mF[i][4] = F.row(1)(1);
    mF[i][5] = F.row(1)(2);
    mF[i][6] = F.row(2)(0);
    mF[i][7] = F.row(2)(1);
    mF[i][8] = F.row(2)(2);

    // For each edge, info is computed
    std::set<unsigned int> edge_1_indices;
    std::set<unsigned int> edge_2_indices;
    std::set<unsigned int> edge_3_indices;

    edge_1_indices.insert(P1_index);
    edge_1_indices.insert(P2_index);

    edge_2_indices.insert(P2_index);
    edge_2_indices.insert(P3_index);

    edge_3_indices.insert(P3_index);
    edge_3_indices.insert(P1_index);

    // First insertion
    if (edge_info.find(edge_1_indices) !=  edge_info.end()) {
      edge_info.insert(std::make_pair(edge_1_indices, std::make_pair(i, -1)));

    }
    // Pair is already there: index of other facet is set
    else {
      edge_info[edge_1_indices].second = i;
    }

    // First insertion
    if (edge_info.find(edge_2_indices) !=  edge_info.end()) {
      edge_info.insert(std::make_pair(edge_2_indices, std::make_pair(i, -1)));
    }
    // Pair is already there: index of other facet is set
    else {
      edge_info[edge_2_indices].second = i;
    }

    // First insertion
    if (edge_info.find(edge_3_indices) !=  edge_info.end()) {
      edge_info.insert(std::make_pair(edge_3_indices, std::make_pair(i, -1)));
    }
    // Pair is already there: index of other facet is set
    else {
      edge_info[edge_3_indices].second = i;
    }

  }

  unsigned int mNOE = edge_info.size();
  assert(mNOE = 3 * (mNOV - 2));

  mListE = new int* [mNOE];

  mE = new double* [mNOE];

  // edge counter is reset
  unsigned int edge = 0;

  for (std::map<std::set<unsigned int>, std::pair<int, int> >::const_iterator edge_iterator = edge_info.begin();
       edge_iterator != edge_info.end();
       ++edge_iterator) {

    mListE[edge] = new int [2];
    mE[edge] = new double [9];

    unsigned int edge_end_index = 0;

    for (const int & edge_end : edge_iterator -> first) {
      mListE[edge][edge_end_index] = edge_end;
      ++edge_end_index;
    }

    // Normal of the two facets owning this edge
    arma::vec n_A = {
      mListN[edge_iterator -> second.first][0],
      mListN[edge_iterator -> second.first][1],
      mListN[edge_iterator -> second.first][2]
    };
    arma::vec n_B = {
      mListN[edge_iterator -> second.second][0],
      mListN[edge_iterator -> second.second][1],
      mListN[edge_iterator -> second.second][2]
    };

    // Unit vector directing the edge
    arma::vec e = {
      mX[mListE[edge][0]] - mX[mListE[edge][1]],
      mY[mListE[edge][0]] - mY[mListE[edge][1]],
      mZ[mListE[edge][0]] - mZ[mListE[edge][1]]
    };
    e = e / arma::norm(e);

    // In-facets edge normals
    arma::vec n_AA = arma::cross(e, n_A);
    n_AA = n_AA / arma::norm(n_AA);

    arma::vec n_BB = arma::cross(e, n_B);
    n_BB = n_BB / arma::norm(n_BB);

    // Edge dyad
    arma::mat E = n_A * n_AA.t() + n_B * n_BB.t();

    mE[edge][0] = E.row(0)(0);
    mE[edge][1] = E.row(0)(1);
    mE[edge][2] = E.row(0)(2);
    mE[edge][3] = E.row(1)(0);
    mE[edge][4] = E.row(1)(1);
    mE[edge][5] = E.row(1)(2);
    mE[edge][6] = E.row(2)(0);
    mE[edge][7] = E.row(2)(1);
    mE[edge][8] = E.row(2)(2);

    ++edge;
  }



}

Asteroid::~Asteroid()
{

  delete[] mX;
  delete[] mY;
  delete[] mZ;

  for (int i = 0; i < mNOF; i++)
  {
    delete[] mListTri[i];
    delete[] mListN[i];
    delete[] mF[i];
  }

  delete[] mListTri;
  delete[] mListN;
  delete[] mF;

  for (int i = 0; i < mNOE; i++)
  {
    delete[] mListE[i];
    delete[] mE[i];
  }

  delete[] mListE;
  delete[] mE;

}

// Get G*sigma value
double Asteroid::GetGs() const
{
  return mGs;
}

// Get No. of Vertices
int Asteroid::GetNOV() const
{
  return mNOV;
}

// Get No. of Facets
int Asteroid::GetNOF() const
{
  return mNOF;
}

// Get No. of Edges
int Asteroid::GetNOE() const
{
  return mNOE;
}

// Get Vertex X coordinates
Vect Asteroid::GetX() const
{
  Vect x(mNOV);

  for (int i = 0; i < mNOV; i++)
  {
    x[i] = mX[i];
  }

  return x;
}

// Get Vertex Y coordinates
Vect Asteroid::GetY() const
{
  Vect y(mNOV);

  for (int i = 0; i < mNOV; i++)
  {
    y[i] = mY[i];
  }

  return y;
}

// Get Vertex Z coordinates
Vect Asteroid::GetZ() const
{
  Vect z(mNOV);

  for (int i = 0; i < mNOV; i++)
  {
    z[i] = mZ[i];
  }

  return z;
}


Vect Asteroid::GetListTri() const
{

  Vect ListTri(mNOF * 3);

  for (int i = 0; i < mNOF; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      ListTri[3 * i + j] = (double) mListTri[i][j];
    }
  }

  return ListTri;
}


Vect Asteroid::GetListN() const
{

  Vect ListN(mNOF * 3);

  for (int i = 0; i < mNOF; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      ListN[3 * i + j] = mListN[i][j];
    }
  }

  return ListN;
}


Vect Asteroid::GetF() const
{

  Vect F(mNOF * 9);

  for (int i = 0; i < mNOF; i++)
  {
    for (int j = 0; j < 9; j++)
    {
      F[9 * i + j] = mF[i][j];
    }
  }

  return F;
}


Vect Asteroid::GetListE() const
{

  Vect ListE(mNOE * 2);

  for (int i = 0; i < mNOE; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      ListE[2 * i + j] = mListE[i][j];
    }
  }

  return ListE;
}


Vect Asteroid::GetE() const
{

  Vect E(mNOE * 9);

  for (int i = 0; i < mNOE; i++)
  {
    for (int j = 0; j < 9; j++)
    {
      E[9 * i + j] = mE[i][j];
    }
  }

  return E;
}


Vect PolyGrav(Vect& Xsc, Asteroid& Body)
{

  Vect a_grav(3);

  int v1, v2, v3;
  Vect r1(3), r2(3), r3(3);
  double R1, R2, R3;
  Vect  nf(3);
  Matrix F(3, 3);
  int l;
  double wf;

  //Sum over Polyhedron Facets
  for (int i = 0; i < Body.mNOF; i++)
  {
    // Vertex no. 1
    v1 = Body.mListTri[i][0];
    r1[0] = Body.mX[v1 - 1];
    r1[1] = Body.mY[v1 - 1];
    r1[2] = Body.mZ[v1 - 1];

    // Vertex no. 2
    v2 = Body.mListTri[i][1];
    r2[0] = Body.mX[v2 - 1];
    r2[1] = Body.mY[v2 - 1];
    r2[2] = Body.mZ[v2 - 1];

    // Vertex no. 3
    v3 = Body.mListTri[i][2];
    r3[0] = Body.mX[v3 - 1];
    r3[1] = Body.mY[v3 - 1];
    r3[2] = Body.mZ[v3 - 1];

    // Normal Unit Vector
    nf[0] = Body.mListN[i][0];
    nf[1] = Body.mListN[i][1];
    nf[2] = Body.mListN[i][2];

    // Face Dyad
    l = 0;
    for (int j = 0; j < 3; j++)
    {
      for (int k = 0; k < 3; k++)
      {
        F(j, k) = Body.mF[i][l];
        l++;
      }
    }


    //Pos vector wrt S/C
    r1 = r1 - Xsc;
    r2 = r2 - Xsc;
    r3 = r3 - Xsc;

    R1 = norm(r1);
    R2 = norm(r2);
    R3 = norm(r3);

    wf = 2 * atan2(dot(r1, cross(r2, r3)), R1 * R2 * R3 + R1 * dot(r2, r3) + R2 * dot(r3, r1) + R3 * dot(r1, r2));

    a_grav = a_grav + F * r1 * wf;
  }


  // Sum over the polyhedron edges
  int e1, e2;
  double Re, Le;
  Matrix E(3, 3);

  for (int i = 0; i < Body.mNOE; i++)
  {

    e1 = Body.mListE[i][0];
    e2 = Body.mListE[i][1];

    //Vertex 1 pos vector from spacecraft
    r1[0] = Body.mX[e1 - 1];
    r1[1] = Body.mY[e1 - 1];
    r1[2] = Body.mZ[e1 - 1];
    r1 = r1 - Xsc;

    R1 = norm(r1);

    //Vertex 2 pos vector from s/c
    r2[0] = Body.mX[e2 - 1];
    r2[1] = Body.mY[e2 - 1];
    r2[2] = Body.mZ[e2 - 1];
    r2 = r2 - Xsc;

    R2 = norm(r2);

    Re = norm(r2 - r1);

    Le = log((R1 + R2 + Re) / (R1 + R2 - Re));

    // Edge dyad
    l = 0;
    for (int j = 0; j < 3; j++)
    {
      for (int k = 0; k < 3; k++)
      {
        E(j, k) = Body.mE[i][l];
        l++;
      }
    }

    // Gravitational Acceleration
    a_grav = a_grav - E * r1 * Le;
  }

  return Body.mGs * a_grav;

}
