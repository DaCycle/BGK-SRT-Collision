# include "d2q9.h"

int main()
{
	// Define Initial f vecotr
	Matrix<double, 1, 9> f1;
	f1 << 1.63, 0.61, 0.41, 0.27, 0.41, 0.15, 0.07, 0.07, 0.16;

	// Call the moment_rho_U_d2q9 function
	auto [Rho1, U1] = moment_rho_U_d2q9(f1);
	cout << "Rho1: " << Rho1 << endl;
	cout << "U1: " << U1.transpose() << endl;

	// Call the eqm_d2q9 function
	Matrix<double, 1, 9> f_eq1;
	f_eq1 = eqm_d2q9(Rho1, U1);

	// Output the result
	cout << "f_eq1:" << endl << f_eq1 << endl;
}