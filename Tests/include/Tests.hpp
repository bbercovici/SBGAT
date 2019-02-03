/** MIT License

Copyright (c) 2018 Benjamin Bercovici and Jay McMahon

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


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
void test_sbgat_transform_shape();

void test_spherical_harmonics_coefs_consistency();
void test_spherical_harmonics_partials_consistency();
void test_sbgat_shape_uq();


void test_radar_obs();
void test_lightcurve_obs();
void test_frame_conversion();



}

#endif