#include "d2q9.cuh"

__global__
void eqm_d2q9(double* rho, double* U, double* f_eq, int N_x, int N_y) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	if (idx >= N_x * N_y * 9) return; // Ensure we don't access out of bounds

	for (int i = idx; i < N_x * N_y * 9; i += stride) {
		int I = i / 9; // Index for rho and U
		int d = i % 9; // Direction index

		double U_x = U[I]; // x-component of velocity
		double U_y = U[I + N_x * N_y]; // y-component of velocity
		double Ksi_U = Ksi[d][0] * U_x + Ksi[d][1] * U_y; // lattice velocity dot product with fluid velocity
		double UU = U_x * U_x + U_y * U_y; // squared velocity magnitude

		f_eq[i] = w[d] * rho[I] * (1.0 + 3.0 * Ksi_U + 4.5 * Ksi_U * Ksi_U - 1.5 * UU);
	}
}
		/*switch (d) {
			case 0:
				f_eq[i] = 4.0 / 9.0 * rho[I] * (1.0 - 1.5 * (U_x * U_x + U_y * U_y));
				break;
			case 1:
				f_eq[i] = 1.0 / 9.0 * rho[I] * (1.0 + 3.0 * U_x + 4.5 * U_x * U_x - 3.0 /2.0 * (U_x * U_x + U_y * U_y));
				break;
			case 2:
				f_eq[i] = 1.0 / 9.0 * rho[I] * (1.0 + 3.0 * U_y + 4.5 * U_y * U_y - 1.5 * (U_x * U_x + U_y * U_y));
				break;
			case 3:
				f_eq[i] = 1.0 / 9.0 * rho[I] * (1.0 - 3.0 * U_x + 4.5 * U_x * U_x - 1.5 * (U_x * U_x + U_y * U_y));
				break;
			case 4:
				f_eq[i] = 1.0 / 9.0 * rho[I] * (1.0 - 3.0 * U_y + 4.5 * U_y * U_y - 1.5 * (U_x * U_x + U_y * U_y));
				break;
			case 5:
				f_eq[i] = 1.0 / 36.0 * rho[I] * (1.0 + 3.0 * (U_x + U_y) + 4.5 * (U_x + U_y) * (U_x + U_y) - 1.5 * (U_x * U_x + U_y * U_y));
				break;
			case 6:
				f_eq[i] = 1.0 / 36.0 * rho[I] * (1.0 + 3.0 * (-U_x + U_y) + 4.5 * (-U_x + U_y) * (-U_x + U_y) - 1.5 * (U_x * U_x + U_y * U_y));
				break;
			case 7:
				f_eq[i] = 1.0 / 36.0 * rho[I] * (1.0 + 3.0 * (-U_x - U_y) + 4.5 * (-U_x - U_y) * (-U_x - U_y) - 1.5 * (U_x * U_x + U_y * U_y));
				break;
			case 8:
				f_eq[i] = 1.0 / 36.0 * rho[I] * (1.0 + 3.0 * (U_x - U_y) + 4.5 * (U_x - U_y) * (U_x - U_y) - 1.5 * (U_x * U_x + U_y * U_y));
				break;
		}
	}
}*/