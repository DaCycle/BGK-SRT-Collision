# include "d2q9.h"

Matrix<double, 1, 9> eqm_d2q9(double Rho, const Matrix<double, 2, 1>& U)
{
	// Define lateral velocity and weight
	Matrix<double, 9, 2> Ksi;
	Ksi << 0,  0,
		   1,  0,
		   0,  1,
		  -1,  0,
		   0, -1,
		   1,  1,
		  -1,  1,
		  -1, -1,
		   1, -1;
	Matrix<double, 9, 1> w;
	w << 4.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36;

	// Calculate the equilibrium distribution function for D2Q9 lattice Boltzmann model
	Matrix<double, 1, 9> f_eq;
	f_eq = ((w * Rho).array() * 
		(1 + (3*Ksi*U).array() +
		4.5 * ((Ksi*U).array()).pow(2) -
		1.5 * (Array<double, 9, 1>::Constant(U.squaredNorm())))).transpose();

	// Return the equilibrium distribution function value
	return f_eq;
}