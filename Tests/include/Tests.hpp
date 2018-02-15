/**
 @file   Test.hpp
 @Author Benjamin Bercovici (bebe0705@colorado.edu)
 @date   August, 2017
 @brief  Declaration of a series of test validating the implementation of the 
 Polyhedron Gravity Model computation
*/


#ifndef HEADER_TESTS
#define HEADER_TESTS

namespace TestsSBCore {

void run();
void test_loading_shape();
void test_sbgat_mass_properties();
void test_sbgat_pgm();
void test_pgm_consistency_cube();
void test_pgm_consistency_ellipsoid();
void test_spherical_harmonics_consistency();
void test_spherical_harmonics_invariance();
void test_geometrical_measures();

}

#endif