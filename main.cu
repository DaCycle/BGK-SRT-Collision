#include "d2q9.cuh"

__constant__ double w[9];
__constant__ double Ksi[9][2];

int main()
{
	// Initialize Constants
	double h_w[9] = { 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 };
	double h_Ksi[9][2] = { { 0.0, 0.0 }, { 1.0, 0.0 }, { 0.0, 1.0 }, { -1.0, 0.0 }, { 0.0, -1.0 }, { 1.0, 1.0 }, { -1.0, 1.0 }, { -1.0, -1.0 }, { 1.0, -1.0 } };
	cudaMemcpyToSymbol(w, h_w, sizeof(h_w));
	cudaMemcpyToSymbol(Ksi, h_Ksi, sizeof(h_Ksi));

	// Set different grid sizes
	const int sizes = 6;
	const int G[sizes] = { 12, 30, 90, 300, 600, 900 };
	int runs = 100;
	
	for (int i = 0; i < sizes; ++i)
	{
		float moment_ms_vec = 0;
		float eqm_ms_vec = 0;

		// Set size of the grid
		const int N_x = G[i];
		const int N_y = N_x;

		// Allocate CPU memory
		double* h_f = new double[N_x * N_y * 9];
		double* h_f_eq = new double[N_x * N_y * 9];
		double* h_Rho = new double[N_x * N_y];
		double* h_U = new double[N_x * N_y * 2];

		// Initialize the grid
		// f(k, j, i) = h_f(k + (j * 9) + (i * 9 * N_x)) <-- Remember moment_Starts from 0
		double f1[9] = { 1.63, 0.61, 0.41, 0.27, 0.41, 0.15, 0.07, 0.07, 0.16 };
		double f2[9] = { 1.67, 0.42, 0.42, 0.42, 0.42, 0.10, 0.11, 0.10, 0.11 };
		for (int i = 0; i < N_y; ++i) {
			for (int j = 0; j < N_x; ++j) {
				for (int k = 0; k < 9; ++k) {
					if (j % 2 == 0) {
						h_f[i * N_x * 9 + j * 9 + k] = f1[k];
					}
					else {
						h_f[i * N_x * 9 + j * 9 + k] = f2[k];
					}
				}
			}

		}

		// Allocate GPU memory
		double* d_f;
		double* d_f_eq;
		double* d_rho;
		double* d_U;
		cudaMalloc((void**)&d_f, N_x * N_y * 9 * sizeof(double));
		cudaMalloc((void**)&d_f_eq, N_x * N_y * 9 * sizeof(double));
		cudaMalloc((void**)&d_rho, N_x * N_y * sizeof(double));
		cudaMalloc((void**)&d_U, N_x * N_y * 2 * sizeof(double));

		// Copy grid from host to device
		cudaMemcpy(d_f, h_f, N_x * N_y * 9 * sizeof(double), cudaMemcpyHostToDevice);

		for (int j = 0; j < runs; ++j) {

			// start moment kernel timing
			cudaEvent_t moment_Start, moment_Stop;
			cudaEventCreate(&moment_Start);
			cudaEventCreate(&moment_Stop);
			cudaEventRecord(moment_Start);

			// Call moment_rho_u_d2q9 kernel
			int blockSize = 256;
			int numBlocks = (N_x * N_y + blockSize - 1) / blockSize;
			moment_rho_u_d2q9 << <numBlocks, blockSize >> > (d_f, d_rho, d_U, N_x, N_y);

			// Record the moment_Stop event for the kernel execution
			cudaEventRecord(moment_Stop);
			cudaEventSynchronize(moment_Stop);
			if (j > 4) { // Skip first 5 iterations for warm-up
				// Calculate elapsed time
				float moment_ms;
				cudaEventElapsedTime(&moment_ms, moment_Start, moment_Stop);
				moment_ms_vec += moment_ms;
			}
		}

		for (int j = 0; j < runs; ++j) {

			// start eqm kernel timing
			cudaEvent_t eqm_Start, eqm_Stop;
			cudaEventCreate(&eqm_Start);
			cudaEventCreate(&eqm_Stop);
			cudaEventRecord(eqm_Start);

			// Call eqm_d2q9 kernel
			int blockSize = 256;
			int numBlocks = (N_x * N_y * 9 + blockSize - 1) / blockSize;
			eqm_d2q9 << <numBlocks, blockSize >> > (d_rho, d_U, d_f_eq, N_x, N_y);

			// Record the eqm_Stop event for the kernel execution
			cudaEventRecord(eqm_Stop);
			cudaEventSynchronize(eqm_Stop);
			if (j > 4) { // Skip first 5 iterations for warm-up
				// Calculate elapsed time
				float eqm_ms;
				cudaEventElapsedTime(&eqm_ms, eqm_Start, eqm_Stop);
				eqm_ms_vec += eqm_ms;
			}
		}

		// Wait for GPU to finish before accessing on host
		cudaDeviceSynchronize();

		// Copy results back to host
		cudaMemcpy(h_f_eq, d_f_eq, N_x * N_y * 9 * sizeof(double), cudaMemcpyDeviceToHost);

		// Output Sample Results
		//cout << "Rho: " << h_Rho[0] << ", " << h_Rho[1] << ", " << h_Rho[2] << ", " << h_Rho[3] << "\n";
		//cout << "U: (" << h_U[0] << ", " << h_U[0 + N_x * N_y] << "), (" << h_U[1] << ", " << h_U[1 + N_x * N_y] << ")\n";
		//cout << "Time taken for moment kernel execution: " << moment_ms << " ms\n";
		cout << "f_eq(:, 1, 2): " << h_f_eq[9] << ", " << h_f_eq[10] << ", " << h_f_eq[11] << ", " << h_f_eq[12] << ", " 
			 << h_f_eq[13] << ", " << h_f_eq[14] << ", " << h_f_eq[15] << ", " << h_f_eq[16] << ", " 
			<< h_f_eq[17] << "\n";

		// Output average time taken for each grid size
		cout << "Average time taken for grid size " << N_x << ": " << moment_ms_vec / (runs - 5.0f) << " ms\n";
		cout << "Average time taken for eqm kernel execution for grid size " << N_x << ": " << eqm_ms_vec / (runs - 5.0f) << " ms\n";

		// Free allocated memory
		delete[] h_f;
		delete[] h_Rho;
		delete[] h_U;
		cudaFree(d_f);
		cudaFree(d_rho);
		cudaFree(d_U);
	}

	return 0;
}