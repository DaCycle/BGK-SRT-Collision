#include "d2q9.cuh"

__global__ 
void collision_d2q9(double* f_new, double* f_eq, double* f, int N_x, int N_y) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	if (idx >= N_x * N_y * 9) return; // Ensure we don't access out of bounds

	for (int i = idx; i < N_x * N_y * 9; i += stride) {
		f[i] = f_new[i] - (f_new[i] - f_eq[i]) / Tau;
	}
}