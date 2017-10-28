#include    "mex.h"
#include  <math.h>
#include  <stdarg.h>
#include  <stdio.h>
#include  <stdlib.h>
#include   <string.h>


//////////////////////////////////////////////////////////////////////////
// File name: GravAccelPartialsExterior_mex.c
//
// Author: Siamak Hesar             3/22/2015
//
// This function is adapted from the original code written by Yu Takahashi 
// called AccelInteriorPotential_mex.c. However I implemented Cunningham's approach
// in computing first and second partial derivatives of the potential w.r.t
// the Cartesian coordinates. The original code used Hotine's method.
// 
///////////////////////////////////////////////////////////////////////////
//
// Description:
//
//  This function computes the acceleration, and dynamics 
//  matrix for the exterior potential. It returns the full partial matrix 
//  of acceleration w.r.t state and gravity harmonics partials.
//
// Inputs:
//
//     n_degree   [n.d.]     : Degree of the spherical harmonics
//
//     ref_radius      [km]       : Reference distance. usually the reference radius
//
//     mu      	  [km^3/sec^2] : Gravitational Parameter.
//
//     r_vec      [km]       : Field point (spacecraft) vector = [x_sat, y_sat, z_sat]
//
//     Cbar       [n.d.]     : Normalized C spherical harmonics
//
//     Sbar       [n.d.]     : Normalized S spherical harmonics
//
// Outputs:
//
//     Accel              : Acceleration by basis \bar{K}_{nm}^i (normalized)
//
//     Jacobian               : STM by \bar{b}_{nm}^i
//
//     STM_C_bar_1              : STM of Cbar spherical harmonics by \bar{b}_{nm}^i
//
//     STM_S_bar_1              : STM of Sbar spherical harmonics by \bar{b}_{nm}^i
//
// Assumptions/References:
//	- 1: S. V. Bettadpur, "Hotine's geopotential formulation: revisited", Bulletin Geodesique (1995) 69:i35-142
//  - 2: R. A. Werner, "Evaluating Descent and Ascent Trajectories Near Non-Spherical Bodies", Technical Support Package
//	- 3: L. E. Cunningham, "On the computation of the spherical harmonic terms needed during the numerical integration of the orbital motion of an artificial satellite"
//	
//  - Mathematical Formulation
//
//    ~ Note the following definitions
//
//       (1) b_{n,m}^e       = (ref_radius/r)^{n+1} * Pnm * [cos(m*lambda); sin(m*lambda)]
//
//       (2) \bar{b}_{n,m}^e = \sqrt{(2 - \delta_{0,m}) * (2n + 1) * (n - m)! / (n + m)!} * (ref_radius/r)^{n+1} * Pnm * [cos(m*lambda); sin(m*lambda)]
//
//       (3) c_{n,m}^e       = ( 2 - \delta_{0,m} )*(n - m)!/(n + m)!*(r'/ref_radius)^n * Pnm * [cos(m*lambda'); sin(m*lambda')]
//
//       (4) \bar{c}_{n,m}^e = \sqrt{ (2 - \delta_{0,m})*(n - m)!/( (2n + 1)*(n + m)! ) } *(r'/ref_radius)^n * Pnm * [cos(m*lambda'); sin(m*lambda')]
//
//      where ' indicates the parameters of the differential mass.
//
//    ~ Note that these expressions can be considered as imaginary numbers. That is,
//
//       (1) b_{n,m}^i       = (ref_radius/r)^{n+1} * Pnm * e^{i*m*lambda}
//
//       (2) \bar{b}_{n,m}^i = \sqrt{(2 - \delta_{0,m}) * (2n + 1) * (n - m)! / (n + m)!} * (ref_radius/r)^{n+1} * Pnm * e^{i*m*lambda}
//
//       (3) c_{n,m}^i       = ( 2 - \delta_{0,m} )*(n - m)!/(n + m)! * (r'/ref_radius)^n * Pnm * e^{i*m*lambda'}
//
//       (4) \bar{c}_{n,m}^i = \sqrt{ (2 - \delta_{0,m})*(n - m)!/( (2n + 1)*(n + m)! ) } * (r'/ref_radius)^n * Pnm * e^{i*m*lambda'}
//
//    ~ The normalization factor is defined as N, where
//
//        \bar{P}_{n,m} = N*P_{n,m}
//
//      Thus, N is immediately recognized as
//
//         N        = \sqrt{(2 - \delta_{0,m}) * (2n + 1) * (n - m)! / (n + m)!}
//
//    ~ The external potential of the field point is computed as
//
//        (1) U^e       = \frac{G*M_ref*}{ref_radius} * \sum^{\infty}_{n = 0} \sum^n_{m = 0} b_{n,m}^i * (1/M_ref) * \int_M c_{n,m}^i dm'
//
//        (2) U^e       = \frac{G*M_ref*}{ref_radius} * \sum^{\infty}_{n = 0} \sum^n_{m = 0} \bar{b}_{n,m}^i * (1/M_ref) * \int_M \bar{c}_{n,m}^i dm'
//
//    ~ This function has the following recursive formulae:
//
//      -------------------------------------------------------------------
//
//      (1) Basis function          : b_{0,0}^e   = (ref_radius/r) * [1; 0]
//
//      (2) Diagonal recurrences    : b_{n,n}^e   = (2*n - 1) * (ref_radius/r) * [x/r, -y/r; y/r, x/r]*b_{n-1,n-1}^e
//
//      (3) Subdiagonal recurrences : b_{n,n-1}^e = (2*n - 1) * (ref_radius/r) * (z/r) * b_{n-1,n-1}^e
//
//      (4) Vertical recurrences    : b_{n,m}^e   = \frac{2*n - 1}{n - m}* (ref_radius/r) * (z/r) * b_{n-1,m}^e - \frac{n + m - 1}{n - m}*(ref_radius/r)^2*b_{n-2,m}^e
//      -------------------------------------------------------------------
//
//      (1) Basis function          : \bar{b}_{0,0}^e   = (ref_radius/r) * [1; 0]
//
//      (2) Diagonal recurrences    : \bar{b}_{n,n}^e   = \sqrt{(1 + \delta_{1,n})*(2n + 1)/(2n)}* (ref_radius/r) * [x/r, -y/r; y/r, x/r]*\bar{b}_{n-1,n-1}^e
//
//      (3) Subdiagonal recurrences : \bar{b}_{n,n-1}^e = \sqrt{2*n - 1}* (ref_radius/r) * (z/r)*\bar{b}_{n-1,n-1}^e
//
//      (4) Vertical recurrences    : \bar{b}_{n,m}^e   = \frac{4n^2 - 1}{n^2 - m^2}* (ref_radius/r) * (z/r)*\bar{b}_{n-1,m}^e - \frac{(2n + 1)*( (n - 1)^2 - m^2 )}{(2n - 3)*(n^2 - m^2)}*(ref_radius/r)^2*\bar{b}_{n-2,m}^e
//
//      -------------------------------------------------------------------
//
//      (1) Basis function          : c_{0,0}^e   = [1; 0]
//
//      (2) Diagonal recurrences    : c_{n,n}^e   = (1 + \delta_{1,n})/(2n)*[x'/ref_radius, -y'/ref_radius; y'/ref_radius, x'/ref_radius]*c_{n-1,n-1}^e
//
//      (3) Subdiagonal recurrences : c_{n,n-1}^e = (z'/ref_radius)*c_{n-1,n-1}^e
//
//      (4) Vertical recurrences    : c_{n,m}^e   = \frac{2n - 1}{n + m}*(z'/ref_radius)*c_{n-1,m}^e - \frac{n - m - 1}{n + m}*(r'/ref_radius)^2*c_{n-2,m}^e
//
//      -------------------------------------------------------------------
//
//      (1) Basis function          : \bar{c}_{0,0}^e   = [1; 0]
//
//      (2) Diagonal recurrences    : \bar{c}_{n,n}^e   = (2n - 1)*\sqrt{ (1 + \delta_{1,n})/( (2n)*(2n + 1) ) }*[x'/ref_radius, -y'/ref_radius; y'/ref_radius, x'/ref_radius]*\bar{c}_{n-1,n-1}^e
//
//      (3) Subdiagonal recurrences : \bar{c}_{n,n-1}^e = \frac{(2n - 1)}{\sqrt{2n + 1}}*(z'/ref_radius)*\bar{c}_{n-1,n-1}^e
//
//      (4) Vertical recurrences    : \bar{c}_{n,m}^e   = (2n - 1)*\sqrt{\frac{2n - 1}{(2n + 1)*(n^2 - m^2)}}*(z'/ref_radius)*\bar{c}_{n-1,m}^e - \sqrt{\frac{(2n - 3)*( (n - 1)^2 - m^2 )}{(2n + 1)*(n^2 - m^2)}}*(r'/ref_radius)^2*\bar{c}_{n-2,m}^e
//
//      -------------------------------------------------------------------
//
//      where subdiagonal recurrences are determined by m = n - 1 through the vertical recurrences.
//
//      -------------------------------------------------------------------
//
// Note:
//
//  - Modified from the original version 
//
// Dependencies:
//
//  - None
//
// Call
//
//  - None
//
// Called by
//
// - TBD
//
// Modification History:
// 10/27/2017 Benjamin Bercovici Added to SBGAT
// 
//  3/22/2015   Siamak Hesar    Modified from the original code written by Yu Takahashi (AccelInteriorPotential_mex.c)
//                              to compute the normalized accelerations and full partials matrix for an exterior gravity field.
// 
//  3/22/2015   Siamak Hesar    Added statements for validating the mex function input types and sizes.
//////////////////////////////////////////////////////////////////////////


#define ABS(x) ((x) < 0) ? -(x) : (x)

#define G 6.67384E-20

#define n_max  50

#define num_C_max  1325
#define num_S_max  1275

/*//////////////////
// -- Outputs -- //
//////////////////*/

// double *Accel_ptr, *Jacobian_ptr, *Partial_C_ptr, *Partial_S_ptr;

/*/////////////////
// -- Inputs -- //
/////////////////*/

int    n_degree;
// double ref_radius, mu;
// double *r_vec_ptr, x_sat, y_sat, z_sat, r_sat;
// double *Cbar_ptr, *Sbar_ptr;
double Cbar[n_max + 1][n_max + 1], Sbar[n_max + 1][n_max + 1];

/*////////////////
// -- Index -- //
////////////////*/

int ii, jj, mm, nn;
int Index_C, Index_S;  

/*///////////////////////////
// -- Degree and Order -- //
///////////////////////////*/

double n, m;
double delta_0_m, delta_1_n, delta_1_m, delta_2_m;
double g_bar_int_1, g_bar_int_2, g_bar_int_3;
double s_bar_int_1, s_bar_int_2, s_bar_int_3, s_bar_int_4, s_bar_int_5, s_bar_int_6;

/*///////////////////////////////
// -- Satellite Parameters -- //
///////////////////////////////*/

double x_ddot, y_ddot, z_ddot;
double ddU_dxdx, ddU_dxdy, ddU_dxdz;
double ddU_dydy, ddU_dydz, ddU_dzdz;
double K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10;

double dxddot_dx_bbar_1, dxddot_dy_bbar_1, dxddot_dz_bbar_1, dyddot_dz_bbar_1, dzddot_dz_bbar_1;
double dxddot_dC_bbar_1[num_C_max][num_C_max], dxddot_dS_bbar_1[num_S_max][num_S_max];
double dyddot_dC_bbar_1[num_C_max][num_C_max], dyddot_dS_bbar_1[num_S_max][num_S_max];
double dzddot_dC_bbar_1[num_C_max][num_C_max], dzddot_dS_bbar_1[num_S_max][num_S_max];

/*/////////////////
// -- Output -- //
/////////////////*/

/*////////////////////
// -- Functions -- //
////////////////////*/

void GetBnmNormalizedExterior(void);

/***************************************************************
 * main program
 ***************************************************************/


/*///////////////////////////////////////////////////////////////////////////*/


/*///////////////////////////////////////////////////////////////////////////*/