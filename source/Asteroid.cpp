// Adapted from the original work
// of Nicola Baresi, CSML, University of Colorado Boulder (2014)

#include "Asteroid.hpp"

Asteroid::Asteroid(std::ifstream& GravityFile) {

  assert(GravityFile.is_open());

  double Gs;

  GravityFile >> Gs;

  this -> mGs = Gs;

  int NOV, NOF, NOE;

  GravityFile >> NOV >> NOF >> NOE;

  this -> mNOV = NOV;
  this -> mNOF = NOF;
  this -> mNOE = NOE;

  assert(mNOF == 2 * mNOV - 4);
  assert(mNOE == 3 * (mNOV - 2));


  this -> mX = new double [mNOV];
  this -> mY = new double [mNOV];
  this -> mZ = new double [mNOV];

  for (int i = 0; i < this -> mNOV; i++) {
    GravityFile >> this -> mX[i] >> this -> mY[i] >> this -> mZ[i];
  }

  this -> mListTri = new int* [this -> mNOF];
  this -> mListN = new double* [this -> mNOF];
  this -> mF = new double* [this -> mNOF];
  this -> surface_grav = new double * [this -> mNOF];


  for (int i = 0; i < this -> mNOF; i++) {

    this -> mListTri[i] = new int [3];
    this -> mListN[i]   = new double [3];
    this -> mF[i] = new double [9];
    this -> surface_grav[i]  = new double [3];

    GravityFile >> this -> mListTri[i][0] >> this -> mListTri[i][1] >> this -> mListTri[i][2];
    GravityFile >> this -> mListN[i][0] >> this -> mListN[i][1] >> this -> mListN[i][2];
    GravityFile >> this -> mF[i][0] >> this -> mF[i][1] >> this -> mF[i][2]
                >> this -> mF[i][3] >> this -> mF[i][4] >> this -> mF[i][5]
                >> this -> mF[i][6] >> this -> mF[i][7] >> this -> mF[i][8];


  }


  this -> mListE = new int* [this -> mNOE];

  this -> mE = new double* [this -> mNOE];

  for (int i = 0; i < this -> mNOE; i++) {
    this -> mListE[i] = new int [2];

    this -> mE[i] = new double [9];

    GravityFile >> this -> mListE[i][0] >> this -> mListE[i][1];
    GravityFile >> this -> mE[i][0] >> this -> mE[i][1] >> this -> mE[i][2]
                >> this -> mE[i][3] >> this -> mE[i][4] >> this -> mE[i][5]
                >> this -> mE[i][6] >> this -> mE[i][7] >> this -> mE[i][8];


  }

}

Asteroid::Asteroid(vtkSmartPointer<vtkPolyData> input_polydata, double Gs) {

  this -> mGs = Gs;
  this -> mNOV = input_polydata -> GetNumberOfPoints();
  this -> mNOF = input_polydata -> GetNumberOfPolys();

  assert(mNOF == 2 * mNOV - 4);

  std::map<std::set<unsigned int>, std::pair<int, int> > edge_info;
  vtkSmartPointer<vtkIdList> facet_vertex_indices = vtkSmartPointer<vtkIdList>::New();

  this -> mX = new double [this -> mNOV];
  this -> mY = new double [this -> mNOV];
  this -> mZ = new double [this -> mNOV];

  // Vertex coordinates are fetched in
  for (vtkIdType vertex = 0; vertex < this -> mNOV; ++vertex) {
    double p[3];
    input_polydata -> GetPoint(vertex, p);
    this -> mX[vertex] = p[0];
    this -> mY[vertex] = p[1];
    this -> mZ[vertex] = p[2];
  }

  // Facet indices are fetched in
  this -> mListTri = new int * [this -> mNOF];
  this -> mListN = new double * [this -> mNOF];
  this -> mF = new double * [this -> mNOF];
  this -> surface_grav = new double * [this -> mNOF];

  for (unsigned int facet = 0; facet < this -> mNOF; ++facet) {

    this -> mListTri[facet] = new int [3];
    this -> mListN[facet]   = new double [3];
    this -> surface_grav[facet]  = new double [3];

    this -> mF[facet] = new double [9];

    input_polydata -> GetCellPoints(facet, facet_vertex_indices);

    // Vertex indices in facet
    unsigned int P1_index = facet_vertex_indices -> GetId(0);
    unsigned int P2_index = facet_vertex_indices -> GetId(1);
    unsigned int P3_index = facet_vertex_indices -> GetId(2);

    this -> mListTri[facet][0] = P1_index;
    this -> mListTri[facet][1] = P2_index;
    this -> mListTri[facet][2] = P3_index;

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

    arma::vec n = arma::cross(P3 - P2, P1 - P2);
    n = n / arma::norm(n);

    this -> mListN[facet][0] = n(0);
    this -> mListN[facet][1] = n(1);
    this -> mListN[facet][2] = n(2);

    // Facet Dyad
    arma::mat F = n * n.t();

    this -> mF[facet][0] = F.row(0)(0);
    this -> mF[facet][1] = F.row(0)(1);
    this -> mF[facet][2] = F.row(0)(2);
    this -> mF[facet][3] = F.row(1)(0);
    this -> mF[facet][4] = F.row(1)(1);
    this -> mF[facet][5] = F.row(1)(2);
    this -> mF[facet][6] = F.row(2)(0);
    this -> mF[facet][7] = F.row(2)(1);
    this -> mF[facet][8] = F.row(2)(2);

    // For each edge, the indices of the two
    // end vertices are extracted
    // along with the vertices of the two facets the
    // edge belongs to
    std::set<unsigned int> edge_1_indices;
    std::set<unsigned int> edge_2_indices;
    std::set<unsigned int> edge_3_indices;

    edge_1_indices.insert(P1_index);
    edge_1_indices.insert(P2_index);

    edge_2_indices.insert(P2_index);
    edge_2_indices.insert(P3_index);

    edge_3_indices.insert(P3_index);
    edge_3_indices.insert(P1_index);

    // First insertion: if the edge has not been inserted yet
    if (edge_info.find(edge_1_indices) ==  edge_info.end()) {
      edge_info.insert(std::make_pair(edge_1_indices, std::make_pair(facet , -1)));
    }

    // Edge has already been inserted: index of other facet is set
    else {
      edge_info[edge_1_indices].second = facet;
    }

    // First insertion: if the edge has not been inserted yet
    if (edge_info.find(edge_2_indices) ==  edge_info.end()) {
      edge_info.insert(std::make_pair(edge_2_indices, std::make_pair(facet , -1)));
    }
    // Edge has already been inserted: index of other facet is set
    else {
      edge_info[edge_2_indices].second = facet ;
    }

    // First insertion: if the edge has not been inserted yet
    if (edge_info.find(edge_3_indices) ==  edge_info.end()) {
      edge_info.insert(std::make_pair(edge_3_indices, std::make_pair(facet , -1)));
    }
    // Edge has already been inserted: index of other facet is set
    else {
      edge_info[edge_3_indices].second = facet ;
    }

  }


  this -> mNOE = edge_info.size();
  assert(mNOE == 3 * (mNOV - 2));

  this -> mListE = new int * [this -> mNOE];

  this -> mE = new double * [this -> mNOE];

  // edge counter is reset
  unsigned int edge = 0;

  // Iterating over edges
  for (std::map<std::set<unsigned int>, std::pair<int, int> >::const_iterator edge_iterator = edge_info.begin();
       edge_iterator != edge_info.end();
       ++edge_iterator) {

    this -> mListE[edge] = new int [2];
    this -> mE[edge] = new double [9];

    unsigned int edge_end_index = 0;

    for (const int & edge_end : edge_iterator -> first) {
      this -> mListE[edge][edge_end_index] = edge_end;
      ++edge_end_index;
    }


    // Vertices on edge
    arma::vec P1 = {
      this -> mX[mListE[edge][0]],
      this -> mY[mListE[edge][0]],
      this -> mZ[mListE[edge][0]],
    };

    arma::vec P2 = {
      this -> mX[mListE[edge][1]],
      this -> mY[mListE[edge][1]],
      this -> mZ[mListE[edge][1]],
    };

    // Vertices not on edge
    arma::vec P3, P4;
    unsigned int P3_index, P4_index;

    // Indices of facets owning the edge
    unsigned int facet_A = edge_iterator -> second.first;
    unsigned int facet_B = edge_iterator -> second.second;

    if (this -> mListTri[facet_A][0] != mListE[edge][0]
        && this -> mListTri[facet_A][0] != mListE[edge][1]) {
      P3_index = this -> mListTri[facet_A][0] ;
    }

    else if (this -> mListTri[facet_A][1] != mListE[edge][0]
             && this -> mListTri[facet_A][1] != mListE[edge][1]) {
      P3_index = this -> mListTri[facet_A][1] ;
    }

    else {
      P3_index = this -> mListTri[facet_A][2] ;
    }

    if (this -> mListTri[facet_B][0] != mListE[edge][0]
        && this -> mListTri[facet_B][0] != mListE[edge][1]) {
      P4_index = this -> mListTri[facet_B][0] ;
    }

    else if (this -> mListTri[facet_B][1] != mListE[edge][0]
             && this -> mListTri[facet_B][1] != mListE[edge][1]) {
      P4_index = this -> mListTri[facet_B][1] ;
    }

    else {
      P4_index = this -> mListTri[facet_B][2] ;
    }

    P3 = {
      this -> mX[P3_index],
      this -> mY[P3_index],
      this -> mZ[P3_index],
    };

    P4 = {
      this -> mX[P4_index],
      this -> mY[P4_index],
      this -> mZ[P4_index],
    };


    // Normal of the two facets owning this edge
    arma::vec n_A = {
      this -> mListN[facet_A][0],
      this -> mListN[facet_A][1],
      this -> mListN[facet_A][2]
    };

    arma::vec n_B = {
      this -> mListN[facet_B][0],
      this -> mListN[facet_B][1],
      this -> mListN[facet_B][2]
    };

    // Unit vector directing the edge
    arma::vec e = P2 - P1;
    e = e / arma::norm(e);

    // Facet-tangent edge normals
    arma::vec n_AA = arma::cross(e, n_A);
    n_AA = n_AA / arma::norm(n_AA);

    arma::vec n_BB = arma::cross(e, n_B);
    n_BB = n_BB / arma::norm(n_BB);

    // Check if facet-tangent edge normals are correctly
    // oriented
    if (arma::dot (P1 - P3, n_AA) < 0) {
      n_AA = - n_AA;
    }

    if (arma::dot (P1 - P4, n_BB) < 0) {
      n_BB = - n_BB;
    }

    // Edge dyad
    arma::mat E = n_A * n_AA.t() + n_B * n_BB.t();

    this -> mE[edge][0] = E.row(0)(0);
    this -> mE[edge][1] = E.row(0)(1);
    this -> mE[edge][2] = E.row(0)(2);
    this -> mE[edge][3] = E.row(1)(0);
    this -> mE[edge][4] = E.row(1)(1);
    this -> mE[edge][5] = E.row(1)(2);
    this -> mE[edge][6] = E.row(2)(0);
    this -> mE[edge][7] = E.row(2)(1);
    this -> mE[edge][8] = E.row(2)(2);
    ++edge;
  }


}

Asteroid::~Asteroid() {

  delete[] mX;
  delete[] mY;
  delete[] mZ;

  for (int i = 0; i < mNOF; i++) {
    delete[] mListTri[i];
    delete[] mListN[i];
    delete[] mF[i];
    delete[] surface_grav[i];

  }

  delete[] mListTri;
  delete[] mListN;
  delete[] mF;
  delete[] surface_grav;


  for (int i = 0; i < mNOE; i++) {
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

void Asteroid::setmGs(double mGs) {
  assert(mGs > 0);
  this -> mGs = mGs;
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

double ** Asteroid::get_ListN() {
  return this -> mListN;
}

int ** Asteroid::get_ListTri() {
  return this -> mListTri;
}

double * Asteroid::get_X() {
  return this -> mX;
}

double * Asteroid::get_Y() {
  return this -> mY;
}

double * Asteroid::get_Z() {
  return this -> mZ;
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


Vect PolyGrav(Vect& Xsc, Asteroid& Body, bool zero_indexed) {

  int start_index;
  if (zero_indexed == true) {
    start_index = 0;

  }
  else {
    start_index = 1;
  }
  Vect a_grav(3);


  int v1, v2, v3;
  Vect r1(3), r2(3), r3(3);
  double R1, R2, R3;
  Vect  nf(3);
  Matrix F(3, 3);
  int l;
  double wf;

  //Sum over Polyhedron Facets
  for (int i = 0; i < Body.mNOF; i++) {

    // Vertex no. 1
    v1 = Body.mListTri[i][0];
    r1[0] = Body.mX[v1 - start_index];
    r1[1] = Body.mY[v1 - start_index];
    r1[2] = Body.mZ[v1 - start_index];

    // Vertex no. 2
    v2 = Body.mListTri[i][1];
    r2[0] = Body.mX[v2 - start_index];
    r2[1] = Body.mY[v2 - start_index];
    r2[2] = Body.mZ[v2 - start_index];

    // Vertex no. 3
    v3 = Body.mListTri[i][2];
    r3[0] = Body.mX[v3 - start_index];
    r3[1] = Body.mY[v3 - start_index];
    r3[2] = Body.mZ[v3 - start_index];

    // Normal Unit Vector
    nf[0] = Body.mListN[i][0];
    nf[1] = Body.mListN[i][1];
    nf[2] = Body.mListN[i][2];

    // Face Dyad
    l = 0;
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
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
    r1[0] = Body.mX[e1 - start_index];
    r1[1] = Body.mY[e1 - start_index];
    r1[2] = Body.mZ[e1 - start_index];
    r1 = r1 - Xsc;

    R1 = norm(r1);

    //Vertex 2 pos vector from s/c
    r2[0] = Body.mX[e2 - start_index];
    r2[1] = Body.mY[e2 - start_index];
    r2[2] = Body.mZ[e2 - start_index];
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


Vect Asteroid::PolyGrav(Vect & Xsc, bool zero_indexed) {

  int start_index;
  if (zero_indexed == true) {
    start_index = 0;

  }
  else {
    start_index = 1;
  }
  Vect a_grav(3);


  int v1, v2, v3;
  Vect r1(3), r2(3), r3(3);
  double R1, R2, R3;
  Vect  nf(3);
  Matrix F(3, 3);
  int l;
  double wf;

  //Sum over Polyhedron Facets
  for (int i = 0; i < this -> mNOF; i++) {

    // Vertex no. 1
    v1 = this -> mListTri[i][0];
    r1[0] = this -> mX[v1 - start_index];
    r1[1] = this -> mY[v1 - start_index];
    r1[2] = this -> mZ[v1 - start_index];

    // Vertex no. 2
    v2 = this -> mListTri[i][1];
    r2[0] = this -> mX[v2 - start_index];
    r2[1] = this -> mY[v2 - start_index];
    r2[2] = this -> mZ[v2 - start_index];

    // Vertex no. 3
    v3 = this -> mListTri[i][2];
    r3[0] = this -> mX[v3 - start_index];
    r3[1] = this -> mY[v3 - start_index];
    r3[2] = this -> mZ[v3 - start_index];

    // Normal Unit Vector
    nf[0] = this -> mListN[i][0];
    nf[1] = this -> mListN[i][1];
    nf[2] = this -> mListN[i][2];

    // Face Dyad
    l = 0;
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        F(j, k) = this -> mF[i][l];
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

  for (int i = 0; i < this -> mNOE; i++) {

    e1 = this -> mListE[i][0];
    e2 = this -> mListE[i][1];

    //Vertex 1 pos vector from spacecraft
    r1[0] = this -> mX[e1 - start_index];
    r1[1] = this -> mY[e1 - start_index];
    r1[2] = this -> mZ[e1 - start_index];
    r1 = r1 - Xsc;

    R1 = norm(r1);

    //Vertex 2 pos vector from s/c
    r2[0] = this -> mX[e2 - start_index];
    r2[1] = this -> mY[e2 - start_index];
    r2[2] = this -> mZ[e2 - start_index];
    r2 = r2 - Xsc;

    R2 = norm(r2);

    Re = norm(r2 - r1);

    Le = log((R1 + R2 + Re) / (R1 + R2 - Re));

    // Edge dyad
    l = 0;
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        E(j, k) = this -> mE[i][l];
        l++;
      }
    }

    // Gravitational Acceleration
    a_grav = a_grav - E * r1 * Le;
  }

  return this -> mGs * a_grav;

}

void Asteroid::compute_global_pgm() {



  for (unsigned int facet = 0; facet < this -> mNOF; ++facet) {
    unsigned int P1_index = this -> mListTri[facet][0];
    unsigned int P2_index = this -> mListTri[facet][1];
    unsigned int P3_index = this -> mListTri[facet][2];

    double Px = 1. / 3. * (this -> mX[P1_index] + this -> mX[P2_index]
                           + this -> mX[P3_index]);
    double Py = 1. / 3. * (this -> mY[P1_index] + this  -> mY[P2_index]
                           + this -> mY[P3_index]);
    double Pz = 1. / 3. * (this -> mZ[P1_index] + this  -> mZ[P2_index]
                           + this -> mZ[P3_index]);

    Vect Xc(3);
    Xc[0] = Px;
    Xc[1] = Py;
    Xc[2] = Pz;


    Vect acc = this -> PolyGrav(Xc, true);
    this -> surface_grav[facet][0] = acc[0];
    this -> surface_grav[facet][1] = acc[1];
    this -> surface_grav[facet][2] = acc[2];
    std::cout << facet << " /" << this -> mNOF << std::endl;


  }
}

void Asteroid::write_to_obj(std::string filename) {
  std::ofstream output_filestream (filename, std::ofstream::out);
  assert(output_filestream.is_open());

  for (unsigned int vertex = 0; vertex < this -> mNOV; ++vertex) {
    output_filestream << "v" << " " << this -> mX[vertex] << " " << this -> mY[vertex] << " " << this -> mZ[vertex] << std::endl;
  }

  for (unsigned int facet = 0; facet < this -> mNOF; ++facet) {
    output_filestream << "f" << " " << this -> mListTri[facet][0] << " " << this -> mListTri[facet][1] << " " << this -> mListTri[facet][2] << std::endl;
  }

  output_filestream.close();

}

double ** Asteroid::get_surface_grav() {
  return this -> surface_grav;
}
