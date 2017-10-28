/**
 @file   Test.hpp
 @Author Benjamin Bercovici (bebe0705@colorado.edu)
 @date   August, 2017
 @brief  Declaration of a series of test validating the implementation of the 
 Polyhedron Gravity Model computation
*/


#ifndef HEADER_TESTS
#define HEADER_TESTS


#include <ShapeModelImporter.hpp>
#include <ShapeModel.hpp>
#include <DynamicAnalyses.hpp>
#include <Constants.hpp>
#include <RigidBodyKinematics.hpp>

#include <assert.h>

namespace TestsSBCore {

void run();
void test_loading_shape();
void test_pgm_consistency_cube();
void test_pgm_consistency_ellipsoid();
void test_spherical_harmonics_consistency();


}

#endif