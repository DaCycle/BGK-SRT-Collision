#ifndef CUDA
#define CUDA
#include <cuda_runtime.h>
#endif

#ifndef CUBLAS
#define CUBLAS
#include <cublas_v2.h>
#pragma comment(lib, "cublas.lib")
#endif

#ifndef IOSSTREAM
#define IOSSTREAM
#include <iostream>
using namespace std;
#endif

#ifndef D2Q9_H
#define D2Q9_H
// For future
#endif

#ifndef MOMENT_RHO_U_D2Q9_H
#define MOMENT_RHO_U_D2Q9_H
__global__ void moment_rho_u_d2q9(double* f, double* rho, double* U, int g_x, int g_y);
#endif