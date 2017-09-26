#ifndef HEADER_COMPUTEPOLYHEDRALCSBAR
#define HEADER_COMPUTEPOLYHEDRALCSBAR


#include <armadillo>



#define nmDim 70
#define dimMax (nmDim+1)*(nmDim+2)/2

/* macro to access trinomial coefficients */
/* coded as a macro so the compiler can optimize loop invariants */

#define Index(i,j,k,n) ( ((n)-(i)) * ((n)-(i)+1) / 2 + (k) )





/**
structure to store trinomial coefficients
**/
typedef struct {
    int degree;
    double data [dimMax]; /* worst-case dimension */
} Trinomial;


void    CalculateBasicTables (int n_degree,
                              double (&mixingFactors) [nmDim + 1] [nmDim + 1] [nmDim + 1],
                              int (&trinomialCoefficientCount)[nmDim + 1]);
void    CalculateFullyNormalizedTables (int n_degree,
                                        double (&diagonalFactors) [nmDim + 1],
                                        double (&subdiagonalFactors) [nmDim + 1],
                                        double (&vertical1Factors) [nmDim + 1] [nmDim + 1],
                                        double (&vertical2Factors) [nmDim + 1] [nmDim + 1]);
void    CalculateUnnormalizedTables (int n_degree,
                                     double (&diagonalFactors) [nmDim + 1],
                                     double (&subdiagonalFactors) [nmDim + 1],
                                     double (&vertical1Factors) [nmDim + 1] [nmDim + 1],
                                     double (&vertical2Factors) [nmDim + 1] [nmDim + 1]);

double  IntegrateOneSimplex (const Trinomial *tri,
                             double (&mixingFactors) [nmDim + 1] [nmDim + 1] [nmDim + 1]);
void    TriAdd  (Trinomial * const result,
                 const Trinomial * const left,
                 const Trinomial * const right,
                 int (&trinomialCoefficientCount)[nmDim + 1]);
void    TriCopy (Trinomial * const target, const Trinomial * const source,
                 int (&trinomialCoefficientCount) [nmDim + 1]);

void    TriMult (Trinomial * const result, const Trinomial * const left,
                 const Trinomial * const right,
                 int (&trinomialCoefficientCount) [nmDim + 1]);

void    TriMultScalar (Trinomial * const result, const double scalar,
                       int (&trinomialCoefficientCount) [nmDim + 1]);
void    TriPrint (const Trinomial * const tri, const char * words);
void    TriSub  (Trinomial * const result, const Trinomial * const left,
                 const Trinomial * const right,
                 int (&trinomialCoefficientCount) [nmDim + 1]);



/**
 * Accumulate contribution of one simplex to potential coefficients Cnm and Snm.
 *
 * Inputs:  Simplex Cartesian coordinates (x0,y0,z0), (x1,y1,z1), (x2,y2,z2).
 * Outputs: Results for this simplex are added to array arguments Cnm and Snm.
 * Globals: totalMass, ref_radius, density,
 *          diagonalFactors, subdiagonalFactors, vertical1factors, vertical2factors
 *******************************************************************/
void AccumulateOneSimplex (
    int n_degree,
    double Cnm [nmDim + 1] [nmDim + 1],
    double Snm [nmDim + 1] [nmDim + 1],
    double x0, double y0, double z0,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    int (&trinomialCoefficientCount) [nmDim + 1],
    double (&diagonalFactors) [nmDim + 1],
    double (&subdiagonalFactors) [nmDim + 1],
    double (&vertical1Factors) [nmDim + 1] [nmDim + 1],
    double (&vertical2Factors) [nmDim + 1] [nmDim + 1],
    double (&mixingFactors) [nmDim + 1] [nmDim + 1] [nmDim + 1],
    double density,
    double polygon_mass,
    double ref_radius
);

void ComputePolyhedralCS(
    arma::mat & Cnm2f,
    arma::mat & Snm2f,
    int n_degree, 
    double ref_radius,
    double polygon_mass, 
    double density,
    double * r0,
    double * r1,
    double * r2,
    double total_mass,
    bool normalized) ;

#endif