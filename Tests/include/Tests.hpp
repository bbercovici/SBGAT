/**
 @file   Test.hpp
 @Author Benjamin Bercovici (bebe0705@colorado.edu)
 @date   August, 2017
 @brief  Tests validating SBGAT implementation
*/


#ifndef HEADER_TESTS
#define HEADER_TESTS

namespace TestsSBCore {

void run();
void test_sbgat_mass_properties();
void test_sbgat_pgm_speed();
void test_sbgat_pgm();

void test_spherical_harmonics_consistency();
void test_spherical_harmonics_invariance();

}

#endif