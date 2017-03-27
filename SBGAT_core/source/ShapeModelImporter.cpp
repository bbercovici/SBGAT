#include "ShapeModelImporter.hpp"

ShapeModelImporter::ShapeModelImporter(std::string filename, double unit_factor) {
	this -> filename = filename;
	this -> unit_factor = unit_factor;
}

void ShapeModelImporter::load_shape_model(ShapeModel * shape_model ) const {

	std::ifstream ifs(this -> filename);

	if (!ifs.is_open()) {
		std::cout << "There was a problem opening the input file!\n";
		throw;
	}

	std::string line;
	std::vector<arma::vec> vertices;
	std::vector<arma::uvec> facet_vertices;

	std::set<std::set<unsigned int> > edge_vertices_indices;


	std::cout << " Reading " << this -> filename << std::endl;
	while (std::getline(ifs, line)) {

		std::stringstream linestream(line);

		std::string word1, word2, word3;

		char type;
		linestream >> type;

		if (type == '#' || type == 's'  || type == 'o' || type == 'm' || type == 'u') {
			continue;
		}

		else if (type == 'v') {
			double vx, vy, vz;
			linestream >> vx >> vy >> vz;
			arma::vec vertex = {vx, vy, vz};
			vertices.push_back(this -> unit_factor * vertex);

		}

		else if (type == 'f') {
			unsigned int v0, v1, v2;
			linestream >> v0 >> v1 >> v2;

			arma::uvec vertices_in_facet = {v0 - 1, v1 - 1, v2 - 1};
			facet_vertices.push_back(vertices_in_facet);

			std::set<unsigned int> edge_0_vertex_indices;
			edge_0_vertex_indices.insert(v0 - 1);
			edge_0_vertex_indices.insert(v1 - 1);

			std::set<unsigned int> edge_1_vertex_indices;
			edge_1_vertex_indices.insert(v0 - 1);
			edge_1_vertex_indices.insert(v2 - 1);

			std::set<unsigned int> edge_2_vertex_indices;
			edge_2_vertex_indices.insert(v1 - 1);
			edge_2_vertex_indices.insert(v2 - 1);

			edge_vertices_indices.insert(edge_0_vertex_indices);
			edge_vertices_indices.insert(edge_1_vertex_indices);
			edge_vertices_indices.insert(edge_2_vertex_indices);

		}

		else {
			std::cout << " unrecognized type: " << type << std::endl;
			throw;
		}

	}

	std::cout << " Number of vertices: " << vertices.size() << std::endl;
	std::cout << " Number of facets: " << facet_vertices.size() << std::endl;
	std::cout << " Number of edges: " << edge_vertices_indices.size() << std::endl;


	// Vertices are added to the shape model
	std::vector<std::shared_ptr<Vertex>> vertex_index_to_ptr(vertices.size(), nullptr);

	std::cout << std::endl << " Constructing Vertices " << std::endl  ;
	boost::progress_display progress_vertices(vertices.size()) ;

	for (unsigned int vertex_index = 0; vertex_index < vertices.size(); ++vertex_index) {

		std::shared_ptr<arma::vec> coordinates = std::make_shared<arma::vec>(vertices[vertex_index]);

		std::shared_ptr<Vertex> vertex = std::make_shared<Vertex>(Vertex());
		vertex -> set_coordinates(coordinates);

		vertex_index_to_ptr[vertex_index] = vertex;
		shape_model -> add_vertex(vertex);
		++progress_vertices;

	}

	std::cout << std::endl << " Constructing Facets " << std::endl ;

	boost::progress_display progress_facets(facet_vertices.size()) ;

	// Facets are added to the shape model
	for (unsigned int facet_index = 0; facet_index < facet_vertices.size(); ++facet_index) {

		// The vertices stored in this facet are pulled.
		std::shared_ptr<Vertex> v0 = vertex_index_to_ptr[facet_vertices[facet_index][0]];
		std::shared_ptr<Vertex> v1 = vertex_index_to_ptr[facet_vertices[facet_index][1]];
		std::shared_ptr<Vertex> v2 = vertex_index_to_ptr[facet_vertices[facet_index][2]];

		std::vector<std::shared_ptr<Vertex>> vertices;
		vertices.push_back(v0);
		vertices.push_back(v1);
		vertices.push_back(v2);


		// Was invariably getting the same memory address if using
		// std::make_shared. The destructor of ShapeModel will take care
		// of those
		Facet * facet = new Facet(std::make_shared<std::vector<std::shared_ptr<Vertex>>>(vertices));

		shape_model -> add_facet(facet);
		++progress_facets;
	}


	// Edges are added to the shape model
	std::cout << std::endl << " Constructing Edges " << std::endl ;

	boost::progress_display progress_edges(edge_vertices_indices.size()) ;
	for (auto edge_iter = edge_vertices_indices.begin(); edge_iter != edge_vertices_indices.end(); ++edge_iter) {

		std::shared_ptr<Vertex> v0 = vertex_index_to_ptr[*edge_iter -> begin()];
		std::shared_ptr<Vertex> v1 = vertex_index_to_ptr[*std::next(edge_iter -> begin())];

		std::shared_ptr<Edge> edge = std::make_shared<Edge>(v0, v1);

		shape_model -> add_edge(edge);
		++progress_edges;

	}



	// The consistency of the surface normals is checked
	shape_model -> check_normals_consistency();


	// The surface area, volume, center of mass of the shape model
	// are computed
	shape_model -> update_mass_properties();

	// The shape model is shifted so as to have its coordinates
	// expressed in its barycentric frame
	shape_model -> shift(-(*shape_model -> get_center_of_mass()));



}
