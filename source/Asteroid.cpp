// Adapted from the original work
// of Nicola Baresi, CSML, University of Colorado Boulder (2014)

#include "Asteroid.hpp"


Asteroid::Asteroid(vtkSmartPointer<vtkPolyData> polydata, double density) {

  this -> mGs = arma::datum::G * density;
  this -> polydata = polydata;
  vtkSmartPointer<vtkIdTypeArray> polys_ids = this -> polydata -> GetPolys () -> GetData ();

  this -> surface_grav = new double * [this -> polydata -> GetNumberOfPolys()];
  this -> polydata -> GetPolys() -> InitTraversal();

  for (unsigned int facet = 0;
       facet < this -> polydata -> GetNumberOfPolys(); ++facet) {

    // Facet Dyad
    double n_array[3];

    this -> polydata -> GetCellData() -> GetNormals() -> GetTuple(facet,
        n_array);

    arma::vec n = {
      n_array[0],
      n_array[1],
      n_array[2]
    };

    // The facet dyads are stored
    arma::mat F = n * n.t();
    this -> F_vector.push_back(F);

    // Global indices of the three vertices in the facet
    std::array<unsigned int, 3> ids_array{ {
        (unsigned int) * (polys_ids -> GetTuple (4 * facet + 1)),
        (unsigned int) * (polys_ids -> GetTuple (4 * facet + 2)),
        (unsigned int) * (polys_ids -> GetTuple (4 * facet + 3))
      } };

    this -> facet_vertices_ids.push_back(ids_array);

    this -> surface_grav[facet] = new double [3];

    std::set< std::set <unsigned int > > empty_set;

    this -> facet_edge_point_ids.push_back(empty_set);

    vtkSmartPointer<vtkIdList> facet_ids = vtkSmartPointer<vtkIdList>::New();

    this -> polydata -> GetPolys() -> GetNextCell(facet_ids);

    std::set<unsigned int> ids_points_in_edge;

    std::set<unsigned int> edge_0;
    std::set<unsigned int> edge_1;
    std::set<unsigned int> edge_2;

    edge_0.insert(facet_ids -> GetId(0));
    edge_0.insert(facet_ids -> GetId(1));

    edge_1.insert(facet_ids -> GetId(1));
    edge_1.insert(facet_ids -> GetId(2));

    edge_2.insert(facet_ids -> GetId(0));
    edge_2.insert(facet_ids -> GetId(2));

    // This facet is augmented with three sets containing the indices
    // of the vertices held by its three edges
    this -> facet_edge_point_ids[facet].insert(edge_0);
    this -> facet_edge_point_ids[facet].insert(edge_1);
    this -> facet_edge_point_ids[facet].insert(edge_2);

    // The reverse map is constructed to link the edges to the facets
    // they belong to
    this -> edge_point_facet_ids[edge_0].insert(facet);
    this -> edge_point_facet_ids[edge_1].insert(facet);
    this -> edge_point_facet_ids[edge_2].insert(facet);
  }

  // The edge dyads are stored
  // The edges owned by this facet are traversed




  for (const auto & edge_point_ids : this -> edge_point_facet_ids) {

    arma::vec r1 ;
    arma::vec r2 ;

    unsigned int e1 = *edge_point_ids.first.begin();
    unsigned int e2 = *std::next(edge_point_ids.first.begin());
    double p[3];

    //Vertex 1 pos vector from spacecraft
    this -> polydata -> GetPoint(e1, p);
    r1 = {p[0], p[1], p[2]};

    //Vertex 2 pos vector from s/c
    this -> polydata -> GetPoint(e2, p);
    r2 = {p[0], p[1], p[2]};

    // Edge dyad
    // Unit vector directing the edge
    arma::vec e = arma::normalise(r2 - r1);

    // Facet-tangent edge normals
    double n_array[3];
    unsigned int facet_A_id = *edge_point_ids.second.begin();
    unsigned int facet_B_id = *std::next(edge_point_ids.second.begin());

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
    if (edge_point_ids.first.find(third_vertex_id_A) != edge_point_ids.first.end()) {
      third_vertex_id_A = * (polys_ids -> GetTuple (4 * facet_A_id + 2));
      if (edge_point_ids.first.find(third_vertex_id_A) != edge_point_ids.first.end()) {
        third_vertex_id_A = * (polys_ids -> GetTuple (4 * facet_A_id + 3));
      }
    }
    this -> polydata -> GetPoint(third_vertex_id_A, p);

    arma::vec p_A = {p[0], p[1], p[2]};

    unsigned int third_vertex_id_B = * (polys_ids -> GetTuple (4 * facet_B_id + 1));
    if (edge_point_ids.first.find(third_vertex_id_B) != edge_point_ids.first.end()) {
      third_vertex_id_B = * (polys_ids -> GetTuple (4 * facet_B_id + 2));
      if (edge_point_ids.first.find(third_vertex_id_B) != edge_point_ids.first.end()) {
        third_vertex_id_B = * (polys_ids -> GetTuple (4 * facet_B_id + 3));
      }
    }

    this -> polydata -> GetPoint(third_vertex_id_B, p);

    arma::vec p_B = {p[0], p[1], p[2]};


    // Check if facet-tangent edge normals are correctly
    // oriented
    if (arma::dot (p_A - r1, n_AA) > 0) {
      n_AA = - n_AA;
    }

    if (arma::dot (p_B - r1, n_BB) > 0) {
      n_BB = - n_BB;
    }

    // Edge dyad
    arma::mat E = n_A * n_AA.t() + n_B * n_BB.t();

    this -> edge_dyads[edge_point_ids.first] = E;

  }

  this -> spin_rate = 0;
  this -> spin_axis = {0, 0, 1};
}

Asteroid::~Asteroid() {

  for (int facet = 0; facet < this -> polydata -> GetNumberOfPolys(); facet++) {
    delete[] surface_grav[facet];

  }

  delete[] surface_grav;

}

// Get G*sigma value
double Asteroid::GetGs() const {
  return mGs;
}

void Asteroid::setmGs(double mGs) {
  assert(mGs > 0);
  this -> mGs = mGs;
}





arma::vec Asteroid::polygrav_vtk(arma::vec & Xsc) {
  arma::vec a_grav = {0, 0, 0};

  std::set<std::set<unsigned int > > used_edges;

  // Sum over Polyhedron Facets
  for (unsigned int facet = 0;
       facet < this -> polydata -> GetNumberOfPolys();
       ++facet) {

    double p[3];
    arma::vec r1 ;
    arma::vec r2 ;
    arma::vec r3 ;

    // Global indices of the three vertices in the facet
    this -> polydata -> GetPoint(this -> facet_vertices_ids[facet][0], p);
    r1 = {p[0], p[1], p[2]};

    this -> polydata -> GetPoint(this -> facet_vertices_ids[facet][1], p);
    r2 = {p[0], p[1], p[2]};

    this -> polydata -> GetPoint(this -> facet_vertices_ids[facet][2], p);
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

    a_grav = a_grav + this -> F_vector[facet] * r1 * wf;

    for (const auto & edge : this -> facet_edge_point_ids[facet] ) {

      // For all the edges in this facet
      if (used_edges.find(edge) == used_edges.end()) {

        arma::vec r1 ;
        arma::vec r2 ;

        unsigned int e1 = *edge.begin();
        unsigned int e2 = *std::next(edge.begin());
        double p[3];

        //Vertex 1 pos vector from spacecraft
        this -> polydata -> GetPoint(e1, p);
        r1 = {p[0], p[1], p[2]};
        r1 = r1 - Xsc;
        double R1 = arma::norm(r1);

        //Vertex 2 pos vector from s/c
        this -> polydata -> GetPoint(e2, p);
        r2 = {p[0], p[1], p[2]};
        r2 = r2 - Xsc;

        double R2 = arma::norm(r2);
        double Re = arma::norm(r2 - r1);
        double Le = std::log((R1 + R2 + Re) / (R1 + R2 - Re));

        // Gravitational Acceleration
        a_grav = a_grav - this -> edge_dyads[edge] * r1 * Le;
        used_edges.insert(edge);

      }
    }
  }

  return this -> mGs * a_grav;

}

void Asteroid::compute_global_pgm() {


  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  for (unsigned int facet = 0; facet < this -> polydata -> GetNumberOfPolys(); ++facet) {

    double  p [3];
    this -> polydata -> GetPoint( this -> facet_vertices_ids[facet][0], p);
    arma::vec P0 = {p[0], p[1], p[2]};

    this -> polydata -> GetPoint( this -> facet_vertices_ids[facet][1], p);
    arma::vec P1 = {p[0], p[1], p[2]};

    this -> polydata -> GetPoint( this -> facet_vertices_ids[facet][2], p);
    arma::vec P2 = {p[0], p[1], p[2]};


    arma::vec Xc = (P0 + P1 + P2) / 3;

    arma::vec acc = this -> polygrav_vtk(Xc);

    this -> surface_grav[facet][0] = acc(0);
    this -> surface_grav[facet][1] = acc(1);
    this -> surface_grav[facet][2] = acc(2);

    std::cout << std::to_string(facet + 1) << " /" << this -> polydata -> GetNumberOfPolys() << std::endl;

  }

  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  double dif = std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count();
  std::printf ("Elasped time is %lf seconds.\n", dif );


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
