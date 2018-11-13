/**
\file Constants.hpp
\author Benjamin Bercovici
\date October 2018
\brief A collection of common constants
\details Holds a few useful common constants like density expressed in SI units
\copyright MIT License, Benjamin Bercovici and Jay McMahon
*/

#ifndef HEADER_CONSTANTS
#define HEADER_CONSTANTS


namespace SBGATConstants {

// Enumeration storing the average densities of a number of bodies of reference
// All densities are expressed in kg/m^3
enum density {

	earth = 5510,
	moon = 3340,
	kw4_alpha = 1970,
	itokawa = 1900,
	bennu = 1260,
	eros = 2670

};

}


#endif