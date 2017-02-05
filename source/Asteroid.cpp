// Adapted from the original work
// of Nicola Baresi, CSML, University of Colorado Boulder (2014)

#include "Asteroid.hpp"


Asteroid::Asteroid(vtkSmartPointer<vtkPolyData> polydata, double density) {

  this -> mGs = arma::datum::G * density;
  this -> polydata = polydata;


  this -> surface_grav = new double * [this -> polydata -> GetNumberOfPolys()];

  for (unsigned int facet = 0;
       facet < this -> polydata -> GetNumberOfPolys(); ++facet) {

    this -> surface_grav[facet] = new double [3];
    std::set< std::set <unsigned int > > empty_set;
    
    this -> facet_edge_point_ids.push_back(empty_set);

    // The facet is cast to a vtkTriangle to enable edge extraction
    vtkSmartPointer<vtkTriangle> trig = vtkTriangle::SafeDownCast(this -> polydata -> GetCell(facet));

    for (vtkIdType i = 0; i < trig -> GetNumberOfEdges(); i++) {

      vtkSmartPointer<vtkCell> edge = trig -> GetEdge(i);
      vtkSmartPointer<vtkIdList> pointIdList = edge -> GetPointIds();

      std::set<unsigned int> ids_points_in_edge;

      for (vtkIdType p = 0; p < pointIdList -> GetNumberOfIds(); p++) {
        ids_points_in_edge.insert(pointIdList -> GetId(p));
      }

      // This facet is augmented with three sets containing the indices
      // of the vertices held by its three edges
      this -> facet_edge_point_ids[facet].insert(ids_points_in_edge);

      // The reverse map is constructed to link the edges to the facets
      // they belong to
      this -> edge_point_facet_ids[ids_points_in_edge].insert(facet);


    }

  }


  this -> spin_rate = 0;
  this -> spin_axis = {0, 0, 1};
}

Asteroid::~Asteroid() {

  // delete[] mX;
  // delete[] mY;
  // delete[] mZ;

  for (int i = 0; i < this -> polydata -> GetNumberOfPolys(); i++) {
    //   delete[] mListTri[i];
    //   delete[] mListN[i];
    //   delete[] mF[i];
    delete[] surface_grav[i];

  }

  // delete[] mListTri;
  // delete[] mListN;
  // delete[] mF;
  delete[] surface_grav;


  // for (int i = 0; i < mNOE; i++) {
  //   delete[] mListE[i];
  //   delete[] mE[i];
  // }

  // delete[] mListE;
  // delete[] mE;

}

// Get G*sigma value
double Asteroid::GetGs() const {
  return mGs;
}

void Asteroid::setmGs(double mGs) {
  assert(mGs > 0);
  this -> mGs = mGs;
}

// Get No. of Vertices
unsigned int Asteroid::GetNOV() const {
  return this -> polydata -> GetNumberOfPoints();
}

// Get No. of Facets
unsigned int Asteroid::GetNOF() const {
  return this -> polydata -> GetNumberOfPolys();
}

// Get No. of Edges
unsigned int Asteroid::GetNOE() const {
  return mNOE;
}



Vect Asteroid::GetListTri() const {

  Vect ListTri(this -> polydata -> GetNumberOfPolys() * 3);

  for (int i = 0; i < this -> polydata -> GetNumberOfPolys(); i++) {
    for (int j = 0; j < 3; j++) {
      ListTri[3 * i + j] = (double) mListTri[i][j];
    }
  }

  return ListTri;
}


Vect Asteroid::GetListN() const {

  Vect ListN(this -> polydata -> GetNumberOfPolys() * 3);

  for (int i = 0; i < this -> polydata -> GetNumberOfPolys(); i++) {
    for (int j = 0; j < 3; j++) {
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

Vect Asteroid::GetF() const {

  Vect F(this -> polydata -> GetNumberOfPolys() * 9);

  for (int i = 0; i < this -> polydata -> GetNumberOfPolys(); i++) {
    for (int j = 0; j < 9; j++) {
      F[9 * i + j] = mF[i][j];
    }
  }

  return F;
}


Vect Asteroid::GetListE() const {

  Vect ListE(this -> mNOE * 2);

  for (int i = 0; i < this -> mNOE; i++) {
    for (int j = 0; j < 2; j++) {
      ListE[2 * i + j] = mListE[i][j];
    }
  }

  return ListE;
}


Vect Asteroid::GetE() const {

  Vect E(mNOE * 9);

  for (int i = 0; i < mNOE; i++) {
    for (int j = 0; j < 9; j++) {
      E[9 * i + j] = mE[i][j];
    }
  }

  return E;
}



Vect Asteroid::PolyGrav(Vect & Xsc) {
  Vect a_grav(3);


  int v1, v2, v3;
  Vect r1(3), r2(3), r3(3);
  double R1, R2, R3;
  Vect  nf(3);
  Matrix F(3, 3);
  int l;
  double wf;

  //Sum over Polyhedron Facets
  for (int i = 0; i < this -> polydata -> GetNumberOfPolys(); i++) {

    // Vertex no. 1
    v1 = this -> mListTri[i][0];
    r1[0] = this -> mX[v1];
    r1[1] = this -> mY[v1];
    r1[2] = this -> mZ[v1];

    // Vertex no. 2
    v2 = this -> mListTri[i][1];
    r2[0] = this -> mX[v2];
    r2[1] = this -> mY[v2];
    r2[2] = this -> mZ[v2];

    // Vertex no. 3
    v3 = this -> mListTri[i][2];
    r3[0] = this -> mX[v3];
    r3[1] = this -> mY[v3];
    r3[2] = this -> mZ[v3];

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
    r1[0] = this -> mX[e1];
    r1[1] = this -> mY[e1];
    r1[2] = this -> mZ[e1];
    r1 = r1 - Xsc;

    R1 = norm(r1);

    //Vertex 2 pos vector from s/c
    r2[0] = this -> mX[e2];
    r2[1] = this -> mY[e2];
    r2[2] = this -> mZ[e2];
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


arma::vec Asteroid::polygrav_vtk(arma::vec & Xsc) {
  arma::vec a_grav = {0, 0, 0};

  std::set<std::set<unsigned int> > used_edges;



  vtkSmartPointer<vtkIdTypeArray> polys_ids = this -> polydata -> GetPolys () -> GetData ();

  // Sum over Polyhedron Facets
  for (int facet = 0; facet < this -> polydata -> GetNumberOfPolys(); facet++) {

    double p[3];
    arma::vec r1 ;
    arma::vec r2 ;
    arma::vec r3 ;

    // Global indices of the three vertices in the facet
    unsigned int facet_vertices_ids[3];
    facet_vertices_ids[0] = * (polys_ids -> GetTuple (4 * facet + 1));
    facet_vertices_ids[1] = * (polys_ids -> GetTuple (4 * facet + 2));
    facet_vertices_ids[2] = * (polys_ids -> GetTuple (4 * facet + 3));

    this -> polydata -> GetPoint(facet_vertices_ids[0], p);
    r1 = {p[0], p[1], p[2]};

    this -> polydata -> GetPoint(facet_vertices_ids[1], p);
    r2 = {p[0], p[1], p[2]};

    this -> polydata -> GetPoint(facet_vertices_ids[2], p);
    r3 = {p[0], p[1], p[2]};

    //Pos vector wrt S/C
    r1 = r1 - Xsc;
    r2 = r2 - Xsc;
    r3 = r3 - Xsc;

    double R1 = arma::norm(r1);
    double R2 = arma::norm(r2);
    double R3 = arma::norm(r3);

    double wf = 2 * std::atan2(arma::dot(r1, arma::cross(r2, r3)),
                               R1 * R2 * R3 + R1 * arma::dot(r2, r3)
                               + R2 * arma::dot(r3, r1)
                               + R3 * arma::dot(r1, r2));

    // Facet Dyad
    double n_array[3];

    this -> polydata -> GetCellData() -> GetNormals() -> GetTuple(facet,
        n_array);

    arma::vec n = {n_array[0], n_array[1], n_array[2]};

    arma::mat F = n * n.t();

    a_grav = a_grav + F * r1 * wf;


    std::cout << facet << std::endl;
    
    // The edges owned by this facet are traversed
    for (auto edge_point_ids : this -> facet_edge_point_ids[facet]) {
      arma::vec r1 ;
      arma::vec r2 ;

      // If this edge has not been used yet, the gravity acceleration
      // is augmented with the corresponding term
      if (used_edges.find(edge_point_ids) == used_edges.end()) {

        unsigned int e1 = *edge_point_ids.begin();
        unsigned int e2 = *std::next(edge_point_ids.begin());
        double p[3];

        //Vertex 1 pos vector from spacecraft
        this -> polydata -> GetPoint(e1, p);
        r1 = {p[0], p[1], p[2]};
        r1 = r1 - Xsc;
        double R1 = norm(r1);

        //Vertex 2 pos vector from s/c
        this -> polydata -> GetPoint(e2, p);
        r2 = {p[0], p[1], p[2]};
        r2 = r2 - Xsc;
        R2 = arma::norm(r2);

        double Re = arma::norm(r2 - r1);
        double Le = std::log((R1 + R2 + Re) / (R1 + R2 - Re));

        // Edge dyad
        // Unit vector directing the edge
        arma::vec e = arma::normalise(r2 - r1);

        // Facet-tangent edge normals
        double n_array[3];
        unsigned int facet_A_id = *edge_point_facet_ids[edge_point_ids].begin();
        unsigned int facet_B_id = *std::next(edge_point_facet_ids[edge_point_ids].begin());

        this -> polydata -> GetCellData() -> GetNormals() -> GetTuple(facet_A_id,
            n_array);
        arma::vec n_A = {n_array[0], n_array[1], n_array[2]};

        this -> polydata -> GetCellData() -> GetNormals() -> GetTuple(facet_B_id,
            n_array);
        arma::vec n_B = {n_array[0], n_array[1], n_array[2]};

        arma::vec n_AA = arma::normalise(arma::cross(e, n_A));
        arma::vec n_BB = arma::normalise(arma::cross(e, n_B));

        // At this point, it is not guaranteed that n_AA and
        // n_BB are properly oriented

        unsigned int third_vertex_id_A = * (polys_ids -> GetTuple (4 * facet_A_id + 1));
        if (edge_point_ids.find(third_vertex_id_A) != edge_point_ids.end()) {
          third_vertex_id_A = * (polys_ids -> GetTuple (4 * facet_A_id + 2));
          if (edge_point_ids.find(third_vertex_id_A) != edge_point_ids.end()) {
            third_vertex_id_A = * (polys_ids -> GetTuple (4 * facet_A_id + 3));
          }
        }
        this -> polydata -> GetPoint(third_vertex_id_A, p);

        arma::vec p_A = {p[0], p[1], p[2]};

        unsigned int third_vertex_id_B = * (polys_ids -> GetTuple (4 * facet_B_id + 1));
        if (edge_point_ids.find(third_vertex_id_B) != edge_point_ids.end()) {
          third_vertex_id_B = * (polys_ids -> GetTuple (4 * facet_B_id + 2));
          if (edge_point_ids.find(third_vertex_id_B) != edge_point_ids.end()) {
            third_vertex_id_B = * (polys_ids -> GetTuple (4 * facet_B_id + 3));
          }
        }

        this -> polydata -> GetPoint(third_vertex_id_B, p);

        arma::vec p_B = {p[0], p[1], p[2]};


        // Check if facet-tangent edge normals are correctly
        // oriented
        if (arma::dot (p_A - r1 + Xsc, n_AA) > 0) {
          n_AA = - n_AA;
        }

        if (arma::dot (p_B - r2 + Xsc, n_BB) > 0) {
          n_BB = - n_BB;
        }

        // Edge dyad
        arma::mat E = n_A * n_AA.t() + n_B * n_BB.t();

        // Gravitational Acceleration
        a_grav = a_grav - E * r1 * Le;

        used_edges.insert(edge_point_ids);
      }

    }

  }

  return this -> mGs * a_grav;

}

void Asteroid::compute_global_pgm() {

  vtkSmartPointer<vtkIdTypeArray> polys_ids = this -> polydata -> GetPolys () -> GetData ();

  for (unsigned int facet = 0; facet < this -> polydata -> GetNumberOfPolys(); ++facet) {


    unsigned int facet_vertices_ids[3];
    facet_vertices_ids[0] = * (polys_ids -> GetTuple (4 * facet + 1));
    facet_vertices_ids[1] = * (polys_ids -> GetTuple (4 * facet + 2));
    facet_vertices_ids[2] = * (polys_ids -> GetTuple (4 * facet + 3));

    double * p0 = this -> polydata -> GetPoint(facet_vertices_ids[0]);
    double * p1 = this -> polydata -> GetPoint(facet_vertices_ids[1]);
    double * p2 = this -> polydata -> GetPoint(facet_vertices_ids[2]);


    double Px = 1. / 3. * (p0[0] + p1[0] + p2[0]);
    double Py = 1. / 3. * (p0[1] + p1[1] + p2[1]);
    double Pz = 1. / 3. * (p0[2] + p1[2] + p2[2]);

    arma::vec Xc = {Px, Py, Pz};


    arma::vec acc = this -> polygrav_vtk(Xc);
    this -> surface_grav[facet][0] = acc(0);
    this -> surface_grav[facet][1] = acc(1);
    this -> surface_grav[facet][2] = acc(2);
    std::cout << std::to_string(facet + 1) << " /" << this -> polydata -> GetNumberOfPolys() << std::endl;


  }
}


int Asteroid::load_surface_acceleration(std::string filename) {

  std::ifstream ifs (filename, std::ifstream::in);
  std::string line;

  // The number of facets is read
  std::getline (ifs, line);
  int nof = std::stoi(line);

  if (nof != this -> polydata -> GetNumberOfPolys()) {
    ifs.close();
    return 0;
  }

  // Density used for the computation of the PGM is read
  std::getline (ifs, line);
  double density = std::stod(line);
  this -> mGs = density * arma::datum::G;

  // Direction of the spin axis in the body frame
  std::getline (ifs, line);

  std::size_t first_space_loc = line.find_first_of(" ");
  double spin_x = std::stod(line.substr (0, first_space_loc));

  line.erase (line.begin() , line.begin() + first_space_loc + 1);

  std::size_t second_space_loc = line.find_first_of(" ");
  double spin_y = std::stod(line.substr (0, second_space_loc));

  line.erase (line.begin() , line.begin() + second_space_loc + 1);

  double spin_z = std::stod(line);

  this -> spin_axis = {spin_x, spin_y, spin_z};
  assert(arma::norm(this -> spin_axis) > 0);

  this -> spin_axis = this -> spin_axis / arma::norm(this -> spin_axis);

  // Spin rate (assumed constant)
  std::getline (ifs, line);
  this -> spin_rate = std::stod(line);

  // The surface gravity is computed
  for (unsigned int facet = 0; facet < this -> polydata -> GetNumberOfPolys(); ++facet) {
    std::getline (ifs, line);

    std::size_t first_space_loc = line.find_first_of(" ");
    double acc_x = std::stod(line.substr (0, first_space_loc));

    line.erase (line.begin() , line.begin() + first_space_loc + 1);

    std::size_t second_space_loc = line.find_first_of(" ");
    double acc_y = std::stod(line.substr (0, second_space_loc));

    line.erase (line.begin() , line.begin() + second_space_loc + 1);

    double acc_z = std::stod(line);

    this -> surface_grav[facet][0] = acc_x;
    this -> surface_grav[facet][1] = acc_y;
    this -> surface_grav[facet][2] = acc_z;

  }

  ifs.close();
  return 1;

}

void Asteroid::set_density(double density) {
  this -> mGs = density * arma::datum::G;
}

double Asteroid::get_density() {
  return this -> mGs / arma::datum::G;
}

void Asteroid::set_spin_axis(arma::vec & spin_axis) {
  assert(arma::norm(spin_axis) > 0);
  this -> spin_axis = spin_axis / arma::norm(spin_axis);
}


void Asteroid::set_spin_rate(double spin_rate) {
  this -> spin_rate = spin_rate;
}

vtkSmartPointer<vtkPolyData> Asteroid::get_polydata() {
  return this -> polydata;
}

void Asteroid::write_surface_acceleration(std::string filename) {

  std::ofstream output_filestream (filename, std::ofstream::out);
  assert(output_filestream.is_open());

  // The number of facets is written first
  output_filestream  << this -> polydata -> GetNumberOfPolys() << std::endl;

  // The density of the asteroid is written next
  output_filestream << this -> mGs / arma::datum::G << std::endl;

  // The direction of the spin axis is written next
  output_filestream << this -> spin_axis[0] << " " << this -> spin_axis[1] << " " << this -> spin_axis[2] << std::endl;

  // The spin rate is written next
  output_filestream << this -> spin_rate << std::endl;

  // Then, the gravity field is fetched in
  for (unsigned int facet = 0; facet < this -> polydata -> GetNumberOfPolys(); ++facet) {
    output_filestream  << this -> surface_grav[facet][0] << " " << this -> surface_grav[facet][1] << " " << this -> surface_grav[facet][2] << std::endl;
  }

  output_filestream.close();

}

double ** Asteroid::get_surface_grav() {
  return this -> surface_grav;
}


double Asteroid::get_spin_rate() const {
  return this -> spin_rate;
}

arma::vec Asteroid::get_spin_axis() const {
  return this -> spin_axis;
}
