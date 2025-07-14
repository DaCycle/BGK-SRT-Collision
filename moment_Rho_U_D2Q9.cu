#include "d2q9.cuh"

__global__
void moment_rho_u_d2q9(double* f, double* rho, double* U, int g_x, int g_y) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	// Calculate density and velocity
	for (int i = idx; i < g_x * g_y; i += stride)
	{
		rho[i] = f[i * 9] + f[i * 9 + 1] + f[i * 9 + 2] + f[i * 9 + 3] + f[i * 9 + 4] + f[i * 9 + 5] + f[i * 9 + 6] + f[i * 9 + 7] + f[i * 9 + 8];
		U[i] = (f[i * 9 + 1] - f[i * 9 + 3] + f[i * 9 + 5] - f[i * 9 + 6] - f[i * 9 + 7] + f[i * 9 + 8]) / rho[i]; // x
		U[i + g_x * g_y] = (f[i * 9 + 2] - f[i * 9 + 4] + f[i * 9 + 5] + f[i * 9 + 6] - f[i * 9 + 7] - f[i * 9 + 8]) / rho[i]; // y
	}
  }