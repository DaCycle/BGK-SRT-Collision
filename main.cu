#include "d2q9.cuh"

int main()
{
	// Set different grid sizes
	const int G[7] = { 12, 30, 90, 300, 600, 900, 6000 };
	int runs = 100;
	
	for (int i = 0; i < 7; ++i)
	{
		float moment_ms_vec = 0;

		// Set size of the grid
		const int g_x = G[i];
		const int g_y = g_x;

		// Allocate CPU memory
		double* h_f = new double[g_x * g_y * 9];
		double* h_Rho = new double[g_x * g_y];
		double* h_U = new double[g_x * g_y * 2];

		// Initialize the grid
		// f(k, j, i) = h_f(k + (j * 9) + (i * 9 * g_x)) <-- Remember starts from 0
		double f1[9] = { 1.63, 0.61, 0.41, 0.27, 0.41, 0.15, 0.07, 0.07, 0.16 };
		double f2[9] = { 1.67, 0.42, 0.42, 0.42, 0.42, 0.10, 0.11, 0.10, 0.11 };
		for (int i = 0; i < g_y; ++i) {
			for (int j = 0; j < g_x; ++j) {
				for (int k = 0; k < 9; ++k) {
					if (j % 2 == 0) {
						h_f[i * g_x * 9 + j * 9 + k] = f1[k];
					}
					else {
						h_f[i * g_x * 9 + j * 9 + k] = f2[k];
					}
				}
			}

		}

		// Allocate GPU memory
		double* d_f;
		double* d_rho;
		double* d_U;
		cudaMalloc((void**)&d_f, g_x * g_y * 9 * sizeof(double));
		cudaMalloc((void**)&d_rho, g_x * g_y * sizeof(double));
		cudaMalloc((void**)&d_U, g_x * g_y * 2 * sizeof(double));

		// Copy grid from host to device
		cudaMemcpy(d_f, h_f, g_x * g_y * 9 * sizeof(double), cudaMemcpyHostToDevice);

		for (int j = 0; j < runs; ++j) {

			// Start kernel timing
			cudaEvent_t start, stop;
			cudaEventCreate(&start);
			cudaEventCreate(&stop);
			cudaEventRecord(start);

			// Call moment_rho_u_d2q9 kernel
			int blockSize = 256;
			int numBlocks = (g_x * g_y + blockSize - 1) / blockSize;
			moment_rho_u_d2q9 << <numBlocks, blockSize >> > (d_f, d_rho, d_U, g_x, g_y);

			// Record the stop event for the kernel execution
			cudaEventRecord(stop);
			cudaEventSynchronize(stop);
			if (j > 4) { // Skip first 5 iterations for warm-up
				// Calculate elapsed time
				float moment_ms;
				cudaEventElapsedTime(&moment_ms, start, stop);
				moment_ms_vec += moment_ms;
			}

			// Wait for GPU to finish before accessing on host
			cudaDeviceSynchronize();
		}

		// Copy results back to host
		cudaMemcpy(h_Rho, d_rho, g_x * g_y * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(h_U, d_U, g_x * g_y * 2 * sizeof(double), cudaMemcpyDeviceToHost);

		// Output Sample Results
		//cout << "Rho: " << h_Rho[0] << ", " << h_Rho[1] << ", " << h_Rho[2] << ", " << h_Rho[3] << "\n";
		//cout << "U: (" << h_U[0] << ", " << h_U[0 + g_x * g_y] << "), (" << h_U[1] << ", " << h_U[1 + g_x * g_y] << ")\n";
		//cout << "Time taken for moment kernel execution: " << moment_ms << " ms\n";

		// Output average time taken for each grid size
		cout << "Average time taken for grid size " << g_x << ": " << moment_ms_vec / (runs - 5.0f) << " ms\n";

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