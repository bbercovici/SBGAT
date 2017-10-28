#include "PolyhedralSphericalHarmo.hpp"


void ComputePolyhedralCS(arma::mat & Cnm2f,
                         arma::mat & Snm2f,
                         int n_degree,
                         double ref_radius,
                         double polygon_mass,
                         double density,
                         double * r0,
                         double * r1,
                         double * r2,
                         double total_mass,
                         bool normalized) {

	/* result accumulators */
	double Cnm [nmDim + 1] [nmDim + 1], Snm [nmDim + 1] [nmDim + 1];
	double	mixingFactors [nmDim + 1] [nmDim + 1] [nmDim + 1]; /* (i!) (j!) (k!) / (n+3)! */
	int		trinomialCoefficientCount [nmDim + 1]; /* number of coefficients in trinomial for any degree */
	double diagonalFactors [nmDim + 1];
	double subdiagonalFactors [nmDim + 1];
	double vertical1Factors [nmDim + 1] [nmDim + 1];
	double vertical2Factors [nmDim + 1] [nmDim + 1];


	double x0 = r0[0];
	double y0 = r0[1];
	double z0 = r0[2];

	double x1 = r1[0];
	double y1 = r1[1];
	double z1 = r1[2];

	double x2 = r2[0];
	double y2 = r2[1];
	double z2 = r2[2];



	/*** set overall problem size ***/
	if (n_degree < 0 || n_degree > nmDim) {
		throw (std::runtime_error("n_degree " + std::to_string(n_degree) + " outside range [0, " + std::to_string(nmDim) + " ]\n"));
	}

	/*** calculate tables needed by AccumulateOneSimplex ***/

	CalculateBasicTables (n_degree, mixingFactors,
	                      trinomialCoefficientCount);
	if (normalized) {
		CalculateFullyNormalizedTables (n_degree,
		                                diagonalFactors,
		                                subdiagonalFactors,
		                                vertical1Factors,
		                                vertical2Factors);
	}
	else {

		CalculateUnnormalizedTables (n_degree,
		                             diagonalFactors,
		                             subdiagonalFactors,
		                             vertical1Factors,
		                             vertical2Factors);
	}




	// /* or, you might invoke CalculateUnnormalizedTables if you want unnormalized coefficients */

	// /** zero the result accumulators **/

	for ( int n = 0; n <= nmDim; ++n) {
		for ( int m = 0; m <= n; ++m) {
			Cnm [n][m] = 0;
			Snm [n][m] = 0;
		}
	}

	/** loop per simplex of polyhedron **/

	int simplexCount = 1; /* whatever */
	for ( int s = 0; s < simplexCount; ++s) {

		AccumulateOneSimplex (n_degree,
		                      Cnm,
		                      Snm,
		                      x0, y0, z0,  x1, y1, z1,  x2, y2, z2,
		                      trinomialCoefficientCount,
		                      diagonalFactors,
		                      subdiagonalFactors,
		                      vertical1Factors,
		                      vertical2Factors,
		                      mixingFactors,
		                      density,
		                      polygon_mass,
		                      ref_radius);

	} /* for s */

	// /** print the results **/

	// /*// 	printf ("\nPotential coefficients of polyhedron\n"); */
	// /*// 	printf ("total mass %g, reference distance %g, polyhedron density %g\n", */
	// /*// 		polygon_mass, ref_radius, density); */
	// /*//PrintCoefficients (Cnm, Snm); */


	Cnm2f = arma::zeros<arma::mat>(n_degree + 1, n_degree + 1);
	Snm2f = arma::zeros<arma::mat>(n_degree + 1, n_degree + 1);

	for (int m = 0; m < n_degree + 1; ++m) {

		for (int n = m; n < n_degree + 1; ++n) {

			Cnm2f.row(n)(m) = Cnm [n][m];
			Snm2f.row(n)(m) = Snm [n][m];

		} /*// For n */

	} /*// For m */



}

void	AccumulateOneSimplex (
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
) {


	Trinomial xTri, yTri, zTri;

	Trinomial rSquared;
	Trinomial diagonalC, diagonalS;

	/* arrays for vertical recurrence; store (n,m), (n-1,m), and (n-2,m) trinomials */
	Trinomial verticalC [3], verticalS [3];

	/* subscripts which index verticalC, verticalS arrays */
	int prior2  = 0; /* n-2,m */
	int prior1  = 1; /* n-1,m */
	int present = 2; /* n,m */

	/* calculate Jacobian determinant *before* normalizing the coordinates */

	const double det = x0 * (y1 * z2 - y2 * z1) + x1 * (y2 * z0 - z2 * y0) + x2 * (y0 * z1 - y1 * z0);

	/* prepare dimensionless overall scale factor for this simplex */

	const double overallFactor = density * det / polygon_mass;




	if (overallFactor == 0.0) return; /* nothing to do */

	/* normalize coordinates */

	x0 /= ref_radius; y0 /= ref_radius; z0 /= ref_radius;
	x1 /= ref_radius; y1 /= ref_radius; z1 /= ref_radius;
	x2 /= ref_radius; y2 /= ref_radius; z2 /= ref_radius;

	/* initialize trinomial x = x0*X + x1*Y + x2*Z */
	xTri.degree = 1;
	xTri.data [Index (1, 0, 0, 1)] = x0;
	xTri.data [Index (0, 1, 0, 1)] = x1;
	xTri.data [Index (0, 0, 1, 1)] = x2;

	/* initialize trinomial y = y0*X + y1*Y + y2*Z */
	yTri.degree = 1;
	yTri.data [Index (1, 0, 0, 1)] = y0;
	yTri.data [Index (0, 1, 0, 1)] = y1;
	yTri.data [Index (0, 0, 1, 1)] = y2;

	/* initialize trinomial z = z0*X + z1*Y + z2*Z */
	zTri.degree = 1;
	zTri.data [Index (1, 0, 0, 1)] = z0;
	zTri.data [Index (0, 1, 0, 1)] = z1;
	zTri.data [Index (0, 0, 1, 1)] = z2;

	/* calculate trinomial rSquared = x*x + y*y + z*z */
	rSquared.degree = 2;
	rSquared.data [Index (2, 0, 0, 2)] =  x0 * x0 + y0 * y0 + z0 * z0;
	rSquared.data [Index (0, 2, 0, 2)] =  x1 * x1 + y1 * y1 + z1 * z1;
	rSquared.data [Index (0, 0, 2, 2)] =  x2 * x2 + y2 * y2 + z2 * z2;
	rSquared.data [Index (1, 1, 0, 2)] = (x0 * x1 + y0 * y1 + z0 * z1) * 2;
	rSquared.data [Index (0, 1, 1, 2)] = (x1 * x2 + y1 * y2 + z1 * z2) * 2;
	rSquared.data [Index (1, 0, 1, 2)] = (x2 * x0 + y2 * y0 + z2 * z0) * 2;

	for ( int m = 0; m <= n_degree; ++m) {	/** order m **/

		for ( int n = m; n <= n_degree; ++n) {	/** degree n **/


			if (n == m) {	/** diagonal **/

				if (m == 0) {	/* anchor */

					verticalC [prior2].degree = 0;
					verticalC [prior2].data [Index (0, 0, 0, 0)] =  overallFactor;

					verticalS [prior2].degree = 0;
					verticalS [prior2].data [Index (0, 0, 0, 0)] = 0.0;

				} else if (m == 1) {	/* special case */

					const double fac = overallFactor * diagonalFactors [1];

					verticalC [prior2].degree = 1;
					verticalC [prior2].data [Index (1, 0, 0, 1)] = x0 * fac;
					verticalC [prior2].data [Index (0, 1, 0, 1)] = x1 * fac;
					verticalC [prior2].data [Index (0, 0, 1, 1)] = x2 * fac;

					verticalS [prior2].degree = 1;
					verticalS [prior2].data [Index (1, 0, 0, 1)] = y0 * fac;
					verticalS [prior2].data [Index (0, 1, 0, 1)] = y1 * fac;
					verticalS [prior2].data [Index (0, 0, 1, 1)] = y2 * fac;



				} else { /* general diagonal shift using prior column's diagonal */

					Trinomial temp1, temp2;

					TriMult (&temp1, &xTri, &diagonalC, trinomialCoefficientCount);
					TriMult (&temp2, &yTri, &diagonalS, trinomialCoefficientCount);
					TriSub  (&verticalC [prior2], &temp1, &temp2, trinomialCoefficientCount);
					TriMultScalar (&verticalC [prior2], diagonalFactors [m], trinomialCoefficientCount);

					TriMult (&temp1, &yTri, &diagonalC, trinomialCoefficientCount);
					TriMult (&temp2, &xTri, &diagonalS, trinomialCoefficientCount);
					TriAdd  (&verticalS [prior2], &temp1, &temp2, trinomialCoefficientCount);
					TriMultScalar (&verticalS [prior2], diagonalFactors [m], trinomialCoefficientCount);

				} /* general diagonal shift */

				/* accumulate harmonic coefficients due to this simplex */



				Cnm [n][m] += IntegrateOneSimplex (&verticalC [prior2], mixingFactors);
				Snm [n][m] += IntegrateOneSimplex (&verticalS [prior2], mixingFactors);

				

				/* remember this diagonal element for initializing next column */

				TriCopy (&diagonalC, &verticalC [prior2], trinomialCoefficientCount);
				TriCopy (&diagonalS, &verticalS [prior2], trinomialCoefficientCount);

			} /* n == m */ else if (n == m + 1) { /** subdiagonal **/

				TriMult (&verticalC [prior1], &verticalC [prior2], &zTri, trinomialCoefficientCount);
				TriMult (&verticalS [prior1], &verticalS [prior2], &zTri, trinomialCoefficientCount);

				TriMultScalar (&verticalC [prior1], subdiagonalFactors [n], trinomialCoefficientCount);
				TriMultScalar (&verticalS [prior1], subdiagonalFactors [n], trinomialCoefficientCount);

				/* accumulate harmonic coefficients due to this simplex */


				

				Cnm [n][m] += IntegrateOneSimplex (&verticalC [prior1], mixingFactors);
				Snm [n][m] += IntegrateOneSimplex (&verticalS [prior1], mixingFactors);


			} else { /** ordinary vertical recurrence **/

				/* prior2 and prior1 index the (n-2,m) and (n-1,m) trinomials */

				Trinomial temp1, temp2;

				TriMult (&temp1, &verticalC [prior1], &zTri, trinomialCoefficientCount);
				TriMultScalar (&temp1, vertical1Factors [n][m], trinomialCoefficientCount);
				TriMult (&temp2, &verticalC [prior2], &rSquared, trinomialCoefficientCount);
				TriMultScalar (&temp2, vertical2Factors [n][m], trinomialCoefficientCount);
				TriSub (&verticalC [present], &temp1, &temp2, trinomialCoefficientCount);

				TriMult (&temp1, &verticalS [prior1], &zTri, trinomialCoefficientCount);
				TriMultScalar (&temp1, vertical1Factors [n][m], trinomialCoefficientCount);
				TriMult (&temp2, &verticalS [prior2], &rSquared, trinomialCoefficientCount);
				TriMultScalar (&temp2, vertical2Factors [n][m], trinomialCoefficientCount);
				TriSub (&verticalS [present], &temp1, &temp2, trinomialCoefficientCount);

				/* accumulate harmonic coefficients due to this simplex */



				

				Cnm [n][m] += IntegrateOneSimplex (&verticalC [present], mixingFactors);
				Snm [n][m] += IntegrateOneSimplex (&verticalS [present], mixingFactors);
				



				/** cycle by adjusting subscripts **/
				/* what was "prior1" becomes "prior2"; what was "present" becomes "prior1" */

				int i = prior2;
				prior2 = prior1;
				prior1 = present;
				present = i;

			} /* ordinary vertical recurrence */

		} /* for n (row) */

	} /* for m (column) */

}


/*******************************************************************
 * Calculate basic tables
 *******************************************************************/
void	CalculateBasicTables (int n_degree, 
	double (&mixingFactors) [nmDim + 1] [nmDim + 1] [nmDim + 1],
                              int (&trinomialCoefficientCount)[nmDim + 1]) {

	/* local table of factorials 0! to (n+3)! */
	double	factorials [nmDim + 4];
	factorials [0] = 1;
	for ( int i = 1; i <= n_degree + 3; ++i) factorials [i] = factorials [i - 1] * i;

	for ( int n = 0; n <= n_degree; ++n) {

		/* number of coefficients in a trinomial of degree n */
		trinomialCoefficientCount [n] = (n + 1) * (n + 2) / 2;

		/* mixing factors table (i!) (j!) (k!) / (n+3)! */
		for ( int i = 0; i <= n; ++i) {
			for ( int j = 0; j <= n - i; ++j) {
				int k = n - i - j;
				mixingFactors [i] [j] [k] = factorials [i] * factorials [j] * factorials [k]
				                            / factorials [n + 3];
			} /* for j */
		} /* for i */

	} /* for n */
}


/*******************************************************************
 * Calculate coefficient tables used by verticalS in AccumulateOneSimplex
 * for fully normalized verticalS
 *******************************************************************/
void	CalculateFullyNormalizedTables (
    int n_degree,
    double (&diagonalFactors) [nmDim + 1],
    double (&subdiagonalFactors) [nmDim + 1],
    double (&vertical1Factors) [nmDim + 1] [nmDim + 1],
    double (&vertical2Factors) [nmDim + 1] [nmDim + 1]
) {

	/* this diagonal from preceding diagonal */
	/* diagonalFactors [0] not used */
	diagonalFactors [1] = 1.0 / sqrt (3.0);
	for (int n = 2; n <= n_degree; ++n) {
		diagonalFactors [n] =  (2.0 * n - 1) / std::sqrt ( (2. * n * (2 * n + 1)) );

	}


	/* subdiagonal */
	for (int n = 0; n <= n_degree; ++n) {
		subdiagonalFactors [n] =  (2.0 * n - 1) / std::sqrt ( (2. * n + 1) );
	}

	/* vertical recurrence */
	/* n = 0 and n = 1 are handled by anchor and subdiagonal code */
	for (int n = 2; n <= n_degree; ++n) {
		/* m = n-1 and m = n are handled by subdiagonal and diagonal code */
		for (int m = 0; m <= n - 2; ++m) {
			const double denom = (2 * n + 1) * (n + m) * (n - m);
			vertical1Factors [n][m] = (2 * n - 1) * std::sqrt ( (2. * n - 1) / denom);
			vertical2Factors [n][m] = std::sqrt (((2. * n - 3) * (n + m - 1) * (n - m - 1)) / denom);
		} /* for m */
	} /* for n */
}


/*******************************************************************
 * Calculate coefficient tables used by verticalS in AccumulateOneSimplex
 * used by ordinary unnormalized verticalS
 *******************************************************************/
void	CalculateUnnormalizedTables (
    int n_degree,
    double (&diagonalFactors) [nmDim + 1],
    double (&subdiagonalFactors) [nmDim + 1],
    double (&vertical1Factors) [nmDim + 1] [nmDim + 1],
    double (&vertical2Factors) [nmDim + 1] [nmDim + 1]) {


	/* this diagonal from preceding diagonal */
	/* diagonalFactors [0] not used */
	diagonalFactors [1] = 1.0;
	for (int n = 2; n <= n_degree; ++n) {
		diagonalFactors [n] = 1.0 / ((double) (2 * n) );
	}

	/* subdiagonal */
	for (int n = 0; n <= n_degree; ++n) {
		subdiagonalFactors [n] = 1.0;
	}

	/* vertical recurrence */
	/* n = 0 and n = 1 are handled by anchor and subdiagonal code */
	for (int n = 2; n <= n_degree; ++n) {
		/* m = n-1 and m = n are handled by subdiagonal and diagonal code */
		for (int m = 0; m <= n - 2; ++m) {
			vertical1Factors [n][m] = ((double) (2 * n - 1)) / ((double) (n + m));
			vertical2Factors [n][m] = ((double) (n - m - 1)) / ((double) (n + m));
		} /* for m */
	} /* for n */
}


/***************************************************************
 * integrate a trinomial which corresponds to a simplex
 ***************************************************************/
double	IntegrateOneSimplex (const Trinomial *tri,
                             double (&mixingFactors) [nmDim + 1] [nmDim + 1] [nmDim + 1]) {
	const int n = tri->degree;


	double accum = 0;
	double const *pTri = tri->data;
	for (int i = n; i >= 0; --i) {
		for (int k = 0; k <= n - i; ++k) {
			int j = n - i - k;
			accum += mixingFactors [i] [j] [k] * (*pTri++);
		} /* for k */
	} /* for i */




	return accum;
}


/***************************************************************
 * trinomial addition
 * Globals: trinomialCoefficientCount
 ***************************************************************/
void	TriAdd (Trinomial * const result,
                const Trinomial * const left,
                const Trinomial * const right,
                int (&trinomialCoefficientCount)[nmDim + 1]) {
	int n = trinomialCoefficientCount [result->degree = left->degree];
	double *pResult = result->data;
	double const *pLeft = left->data;
	double const *pRight = right->data;
	while (n--) {
		*pResult++ = (*pLeft++) + (*pRight++);
	}
}


/***************************************************************
 * make a copy of a trinomial
 * Globals: trinomialCoefficientCount
 ***************************************************************/
void	TriCopy (Trinomial * const target, const Trinomial * const source,
                 int (&trinomialCoefficientCount) [nmDim + 1]) {
	int n = trinomialCoefficientCount [target->degree = source->degree];
	double *pTarget = target->data;
	double const *pSource = source->data;
	while (n--) {
		*pTarget++ = *pSource++;
	}
}


/***************************************************************
 * trinomial multiplication
 * Globals: trinomialCoefficientCount
 ***************************************************************/
void	TriMult (Trinomial * const result, const Trinomial * const left,
                 const Trinomial * const right,
                 int (&trinomialCoefficientCount) [nmDim + 1]) {

	const  double *pCoeffLeft, *pCoeffRight;
	double *pResult;

	result->degree = left->degree + right->degree;

	/* zero the result before accumulation */
	pResult = result->data;
	int i = trinomialCoefficientCount [result->degree];
	while (i--) *pResult++ = 0;

	pCoeffLeft = left->data;
	for (int iLeft = left->degree; iLeft >= 0; --iLeft) {
		for (int kLeft = 0; kLeft <= left->degree - iLeft; ++kLeft) {

			if (*pCoeffLeft != 0.0) {

				int jLeft = left->degree - iLeft - kLeft;

				pCoeffRight = right->data;
				for (int iRight = right->degree; iRight >= 0; --iRight) {
					for (int kRight = 0; kRight <= right->degree - iRight; ++kRight) {

						if (*pCoeffRight != 0.0) {
							int iResult, jResult, kResult, ndxResult;

							int jRight = right->degree - iRight - kRight;

							iResult = iLeft + iRight;
							jResult = jLeft + jRight;
							kResult = kLeft + kRight;
							ndxResult = Index (iResult, jResult, kResult, result->degree);

							result->data [ndxResult] += (*pCoeffLeft) * (*pCoeffRight);

						} /* coeffRight != 0 */

						++pCoeffRight;

					} /* for kRight */
				} /* for iRight */

			} /* coeffLeft != 0 */

			++pCoeffLeft;

		} /* for kLeft */
	} /* iLeft */
}


/***************************************************************
 * multiply a trinomial by a scalar
 * Globals: trinomialCoefficientCount
 ***************************************************************/
void	TriMultScalar (Trinomial * const result, const double scalar,
                       int (&trinomialCoefficientCount) [nmDim + 1]) {
	int n = trinomialCoefficientCount [result->degree];
	double * pResult = result->data;
	while (n--) {
		*pResult++ *= scalar;
	}
}


/***************************************************************
 * print a trinomial
 ***************************************************************/
void	TriPrint (const Trinomial * const tri, const char *words) {

	/**
	Disabled
	*/
}


/***************************************************************
 * trinomial subtraction
 * Globals: trinomialCoefficientCount
 ***************************************************************/
void	TriSub (Trinomial * const result, const Trinomial * const left,
                const Trinomial * const right,
                int (&trinomialCoefficientCount) [nmDim + 1]) {
	int n = trinomialCoefficientCount [result->degree = left->degree];
	double *pResult = result->data;
	double const *pLeft = left->data;
	double const *pRight = right->data;
	while (n--) {
		*pResult++ = (*pLeft++) - (*pRight++);
	}
}

