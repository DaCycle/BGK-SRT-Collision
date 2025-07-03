# include "d2q9.h"

tuple<double, Matrix<double, 2, 1>> moment_rho_U_d2q9(const Matrix<double, 1, 9>& f)
{
	// Define lattice velocity
	Matrix<double, 2, 9> Ksi;
	Ksi << 0, 1, 0, -1, 0, 1, -1, -1, 1,
		0, 0, 1, 0, -1, 1, 1, -1, -1;

	// Find Rho
	double Rho = f.sum();

	// Find U
	Matrix<double, 2, 1> U;
	U = Ksi * f.transpose() / Rho;

	// Return
	return make_tuple(Rho, U);
}