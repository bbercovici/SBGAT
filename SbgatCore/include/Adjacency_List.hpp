/**
\Adjacency_List.hh
\author Andrew French
\date Jan 8, 2014
 \brief Adjacency_List class template. This is a template for a graph of vertices
        and edges of types T1 and T2, respectively.

\details The Adjacency_List class is a container class that stores information
 in the form of a graph, either directed or undirected, through locally defined
 vertices and edges. This class supports dynamic insertion and removal of both
 vertices and edges, contrary to popular implementations. A basic
 depth-firsth-algorithm is also included, providing simple path-finding
 between vertices, so long as the graph is connected.

*/

#ifndef ADJACENCYLIST_HPP
#define ADJACENCYLIST_HPP

#include <map>
#include <list>
#include <queue>
#include <algorithm>
#include <iostream>
#include <iomanip>
// #include "TPint.hh"
#include <stdexcept>
#include <set>

/*
 * Template type declaration:
 *
 * T1: vertex type (node)
 * T2: edge type (link)
 *
 * ----------------------------- Rule of Three! --------------------------------
 * Type T1 (not T2) MUST satisfy the conditions required for C++ STL constainers.
 * i.e. it must have a destructor, copy constructor, and copy assignment operator.
 *
 */
template<class T1, class T2>
class Adjacency_List {

public:

	/*-------------------- General Types and Constants -----------------------*/

	/**
	 * \brief Vertex type definition for the Adjacency_List graph class
	 *
	 * \details
	 * Each vertex of type T1 contains a list of pointers to the vertices adjacent
	 * to it, paired with an object of type T2. Each vertex also contains a boolean
	 * that is used in the search algorithms.
	 *
	 */
	typedef struct vertex {
		typedef std::pair<T2, vertex*>    edge;
		std::list<edge>                   adj;
		T1                                data;
		bool                              indentified;
		vertex(T1 d) {
			data  = d;
			indentified  = false;
		}
	} vertex_type;

	/**
	 * \brief edge type definition:
	 *
	 * \details
	 * pair of type T2 object and pointer to a vertex
	 *
	 */
	typedef std::pair<T2, vertex*>      edge;

	/**
	 * \brief vertexmap type definition:
	 *
	 * \details
	 * This is the all encompassing data structure. A pointer to each vertex is
	 * stored in a map. A map was chosen for its fast (O(1)) access and (O(log(n)))
	 * search methods (i.e. std::find()).
	 */
	typedef std::map<T1, vertex*>       vertexmap;

	/**
	 * \brief Constructor
	 *
	 * \details Default constructor: creates empty undirected graph.
	 */
	Adjacency_List() {
		this->isdirected = false;
	}

	/**
	 * \brief Destructor
	 */
	~Adjacency_List() {
		typename vertexmap::iterator it;

		// free memory allocated to vertices
		for (it = this->vmap.begin(); it != this->vmap.end(); ++it) {
			delete it->second;
		}
	}

	/*-------------------------- Access Methods ------------------------------*/

	/**
	 * \brief uint16_t getnumv()
	 *
	 * \details Return the number of vertices.
	 */
	uint16_t getnumv() {
		return (this->vmap.size());
	}

	/**
	 * \brief T2 getedge(T1 src, T1 dest)
	 *
	 * \details Return edge between two vertices, if one exists.
	 */
	T2 getedge(T1 src, T1 dest) {

		// local type definitions
		typename vertexmap::iterator it  = (this->vmap.find(src));
		typename std::list<edge>::iterator adjIt;
		bool exist = false;

		// search source's adjacency list for edge
		for (adjIt = it->second->adj.begin(); adjIt != it->second->adj.end(); ++adjIt) {
			if (adjIt->second->data == dest) {
				exist = true;
				break;
			}
		}

		// throw error if edge does not exist
		if (!exist) {
			std::cout << "Error: Adjacency_List::getedge(). Edge between "
			          << src << " and " << dest << " does not exist.\n";
			throw (std::runtime_error(""));
		}

		return (adjIt->first);
	}



	/**
	\brief std::set<T1> getneighbors(T1 src)
	 *
	 * \details Return set of src's neighbors, if src exists.
	 */
	std::set<T1> getneighbors(T1 src) {
		std::set<T1> neighbors;
		if (vertexexists(src) == true) {

			// local type definitions
			typename vertexmap::iterator it  = (this->vmap.find(src));
			typename std::list<edge>::iterator adjIt;




			// search source's adjacency list for edge
			for (adjIt = it->second->adj.begin(); adjIt != it->second->adj.end(); ++adjIt) {
				neighbors.insert(adjIt -> second -> data);
			}
		}
		else {
			std::cout << "Error Adjacency_List::getneighbors(). vertex "
			          << src << " does not exist!\n";
			throw (std::runtime_error(""));

		}

		return neighbors;
	}

	/**
	 * \brief bool edgeexist(T1 src, T1 dest)
	 *
	 * \details Determine if an edge exists between two vertices. To prevent runtime errors
	 * this should ALWAYS be called before getedge().
	 */
	bool edgeexist(T1 src, T1 dest) {

		// local type definitions
		typename vertexmap::iterator it  = (this->vmap.find(src));
		typename std::list<edge>::iterator adjIt;
		bool exist = false;

		// search source's adjacency list for edge
		for (adjIt = it->second->adj.begin(); adjIt != it->second->adj.end(); ++adjIt) {
			if (adjIt->second->data == dest) {
				exist = true;
				break;
			}
		}

		return (exist);
	}


	/**
	 * \brief bool vertexexists(T1 node)
	 *
	 * \details Check if vertex is in graph
	 */
	bool vertexexists(T1 node) {

		// local type definitions
		typename vertexmap::iterator it  = (this->vmap.find(node));

		if (it == vmap.end()) {
			return (false);
		}
		return (true);
	}

	/**
	 * \brief T1 getvertex(T1 node)
	 *
	 * \details Return vertex from graph, if it exists
	 */
	T1 getvertex(T1 node) {

		// local type definitions
		typename vertexmap::iterator it  = (this->vmap.find(node));

		return (it->second->data);
	}

	/**
	 * \brief void displaygraph(uint16_t w)
	 *
	 * \details Display graph to console. Parameter w determines setw formatting.
	 * Could overload << operator
	 */
	void displaygraph(uint16_t w) {

		// local type definitions
		typename vertexmap::iterator        it;
		typename std::list<edge>::iterator  adjIt;

		std::cout << std::endl << std::endl;

		// for all vertices
		for (it = this->vmap.begin(); it != this->vmap.end(); ++it) {

			// display vertex
			std::cout << std::setw(w) << it->second->data << " |";

			// for all edges in adjacency list
			for (adjIt = it->second->adj.begin(); adjIt != it->second->adj.end();
			        ++adjIt) {

				// display edge
				std::cout << std::setw(w) << adjIt->first << " ";

			}
			std::cout << std::endl;
		}
		std::cout << std::endl << std::endl;

	}

	/**
	\brief std::set<T2> get_edges()
	\details Returns a set storing all the existing edges
	*/

	std::set<T2> get_edges() {
		std::set<T2> edges;
		// local type definitions
		typename vertexmap::iterator        it;
		typename std::list<edge>::iterator  adjIt;

		// for all vertices
		for (it = this->vmap.begin(); it != this->vmap.end(); ++it) {

			// for all edges in adjacency list
			for (adjIt = it->second->adj.begin(); adjIt != it->second->adj.end();
			        ++adjIt) {

				// insert edge
				edges.insert(adjIt->first);

			}
		}
		return edges;
	}

	/**
	\brief std::set<T1> get_vertices()
	\details Returns a set storing all the existing vertices
	*/

	std::set<T1> get_vertices() {
		std::set<T1> vertices;
		// local type definitions
		typename vertexmap::iterator        it;

		// for all vertices
		for (it = this -> vmap.begin(); it != this -> vmap.end(); ++it) {
			// insert vertex
			vertices.insert(it -> second -> data);

		}



		return vertices;
	}



	/**
	 * \brief std::deque<T1> dfs(T1 src, T1 dest)
	 *
	 * \details Wrapper around _dfs (depth-first-search) algorithm that provides a
	 * cleaner user interface. See _dfs in the private section to view the
	 * algorithm implementation.
	 *
	 * This function returns a path between two vertices if one exists. If no path
	 * exists, the dq output will be empty, which can be checked externally using
	 * dq.empty(). If more than one path exists, there is NO guarantee that the
	 * path found will be the shorted path.
	 *
	 * Note: std::stack can be used in place of std::deque for less overhead
	 *       however the deque comes with an iterator. this is an easy fix if
	 *       needed later on.
	 *
	 */
	std::deque<T1> dfs(T1 src, T1 dest) {

		// local type definitions
		typename vertexmap::iterator it, src_it, dest_it;
		src_it  = (this->vmap.find(src));
		dest_it = (this->vmap.find(dest));

		// containter
		std::deque<T1> dq;
		// flag that indicates whether or not dest has been found
		// (ie a path has been found)
		bool found = false;

		// check that both vertices exist
		if (src_it == vmap.end()) {
			std::cout << "Error Adjacency_List::dfs(). Source vertex "
			          << src << " does not exist!\n";
			throw (std::runtime_error(""));
		} else if (dest_it == vmap.end()) {
			std::cout << "Error Adjacency_List::dfs(). Destination vertex "
			          << dest << " does not exist!\n";
			throw (std::runtime_error(""));
		}

		// set all vertices to unidentified
		for (it = this->vmap.begin(); it != this->vmap.end(); ++it) {
			it->second->indentified = false;
		}

		// call the depth-first-search algorithm
		_dfs(src, dest, dq, found);
		return (dq);
	}


	/**
	 * \brief void ddfs
	 *
	 * \details double ended (2-way) depth-first-search. Not yet implemented.
	 *
	 * Advantages: faster searches for large graphs
	 */
	void ddfs(T1 src, T1 dest) {
		std::cout << "Error Adjacency_List::ddfs(). This function has not "
		          << "yet been implemented. Sorry!\n";
		throw (std::runtime_error(""));
	}

	/*------------------------- Mutator Methods ------------------------------*/

	/**
	 * \brief void setdirected(bool isdirected)
	 *
	 * Set the directedness of the graph. Default is undirected (false)
	 */
	void setdirected(bool isdirected) {
		this->isdirected = isdirected;
	}

	/**
	 * \brief void addedge(T1 src, T1 dest, T2 link)
	 *
	 * Add edge between two vertices, if both exist.
	 *
	 */
	void addedge(T1 src, T1 dest, T2 link) {

		// local type definitions
		typename vertexmap::iterator src_it, dest_it;
		src_it  = (this->vmap.find(src));
		dest_it = (this->vmap.find(dest));

		// could add check for edge existence here if worried about adding
		// multiple edges between vertices
		// could also check for edge uniqueness to allow for multiple edges
		// between vertices so long as they contain differet values

		// check that both vertices exist
		if (src_it == vmap.end()) {
			std::cout << "Warning Adjacency_List::addedge(). Source vertex "
			          << src << " does not exist! Edge was not added!\n";
		} else if (dest_it == vmap.end()) {
			std::cout << "Warning Adjacency_List::addedge(). Destination vertex "
			          << dest << " does not exist! Edge was not added!\n";
		} else {

			// add edge
			vertex *from = (this->vmap.find(src)->second);
			vertex *to   = (this->vmap.find(dest)->second);

			edge e1 = std::make_pair(link, to);
			from->adj.push_back(e1);

			// if undirected, add 'other' edge
			if (!this->isdirected) {
				edge e2 = std::make_pair(link, from);
				to->adj.push_back(e2);
			}
		}
	}

	/**
	 * \brief void addvertex(T1 node)
	 *
	 * \details Add vertex to graph.
	 *
	 */
	void addvertex(T1 node) {

		// local type definitions
		typename vertexmap::iterator it = this->vmap.begin();
		it = this->vmap.find(node);

		// if vertex does not already exist, add it
		if (it == vmap.end()) {
			vertex *v;
			v = new vertex(node);
			vmap[node] = v;
			return;
		}
		std::cout << "Warning Adjacency_List::addvertex(). Vertex "
		          << node << " already exists!\n";
	}

	/**
	 * \brief void removeedge(T1 src, T1 dest)
	 *
	 * \details Remove edge from graph.
	 *
	 */
	void removeedge(T1 src, T1 dest) {

		// local type definitions
		typename vertexmap::iterator it  = (this->vmap.find(src));
		typename std::list<edge>::iterator adjIt;

		// search for and delete edge
		for (adjIt = it->second->adj.begin(); adjIt != it->second->adj.end();
		        ++adjIt) {
			if (adjIt->second->data == dest) {
				it->second->adj.remove(*adjIt);
				break;
			}
		}

		// if this graph is undirected, delete the 'other' edge
		if (!this->isdirected) {
			it  = (this->vmap.find(dest));
			for (adjIt = it->second->adj.begin(); adjIt != it->second->adj.end();
			        ++adjIt) {
				if (adjIt->second->data == src) {
					it->second->adj.remove(*adjIt);
					break;
				}
			}
		}
	}

	/**
	 * \brief void removevertex(T1 node)
	 *
	 * \details Remove vertex, and all associated edges, from graph.
	 *
	 */
	void removevertex(T1 node) {

		// local type definitions
		typename vertexmap::iterator it;

		// remove all edges from node
		it = (this->vmap.find(node));
		it->second->adj.clear();

		// remove all edges to node
		for (it = this->vmap.begin(); it != this->vmap.end(); ++it) {
			removeedge(it->second->data, node);
		}

		// remove node
		delete vmap[node];
		this->vmap.erase(node);

	}


private:

	vertexmap vmap;
	bool isdirected;

	/**
	 * \brief void _dfs(T1 src, T1 dest, std::deque<T1> &dq, bool &found)
	 *
	 * \details Recursively defined depth-first-search algorithm. The path found is stored
	 * in the deque. If no path exists, the function will not modify the dq
	 * (i.e. it will remain empty if it was passed in empty).
	 *
	 * Neccesary conditions for path retrieval:
	 *
	 * 1) Graph must be connected (path must exist between src & dest)
	 *
	 * Shortest path gaurenteed when:
	 *
	 * 1) Only one path exists -- different algorithm can be implemented to
	 *                            accomodate multiple paths
	 *                            -- a variation of breadth-first-search
	 *
	 * ** This will work for all trees!
	 */
	void _dfs(T1 src, T1 dest, std::deque<T1> &dq, bool &found) {

		// local type definitions
		typename vertexmap::iterator it;
		typename std::list<edge>::iterator adjIt;

		it  = (this->vmap.find(src));

		// mark source as identified
		it->second->indentified = true;

		// for all edges in adjacency list
		for (adjIt = it->second->adj.begin(); adjIt != it->second->adj.end();
		        ++adjIt) {

			// if destination, push dest and src into queue and return
			// path found!
			if (adjIt->second->data == dest) {
				dq.push_front(dest);
				dq.push_front(src);
				found = true;
				return;

				// skip edges back to self
			} else if (adjIt->second->data == src) {
				//... do nothing ...

				// skip previously indentified edges
			} else if (adjIt->second->indentified == true) {
				//... do nothing ...

				// recursively call _dfs on adjacent vertices
			} else {
				_dfs(adjIt->second->data, dest, dq, found);
				// if dest has been found, push entire path into queue
				if (found == true) {
					dq.push_front(src);
					return;
				}
			}
		}
		return;
	}

	/*
	 * Disabled methods: copy constructor, assignment operator
	 */
	Adjacency_List( const Adjacency_List &rhs );
	void operator=( const Adjacency_List &rhs );
};

#endif /* ADJACENCYLIST_HPP */
