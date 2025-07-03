#ifndef VECTOR
#define VECTOR
#include <vector>
#endif

#ifndef EIGEN
#define EIGEN
#include <Eigen/Dense>
using namespace Eigen;
#endif

#ifndef IOSSTREAM
#define IOSSTREAM
#include <iostream>
using namespace std;
#endif

#ifndef TUPLE
#define TUPLE
#include <tuple>
#endif


#ifndef D2Q9_H
#define D2Q9_H
Matrix<double, 1, 9> eqm_d2q9(double Rho, const Matrix<double, 2, 1>& U);
#endif

#ifndef MOMENT_RHO_U_D2Q9_H
#define MOMENT_RHO_U_D2Q9_H
tuple<double, Matrix<double, 2, 1>> moment_rho_U_d2q9(const Matrix<double, 1, 9>& f);
#endif