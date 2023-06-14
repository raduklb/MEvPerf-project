#include "utils.h"
#include <cblas.h>

double* my_solver(int N, double *A, double *B) {
	double *result, *B_copy;
	
	printf("BLAS SOLVER\n");

	// Copy matrix B.
	B_copy = (double*) calloc(N * N, sizeof(double));
	if (B_copy == NULL) {
		exit(1);
	}
	
	cblas_dcopy(N * N, B, 1, B_copy, 1);

	// Compute A * B. The result is stored in B_copy.
	cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans,
				CblasNonUnit, N, N, 1, A, N, B_copy, N);

	// Compute A * B * B' (B_copy * B'). The result is stored in result.
	result = (double*) calloc(N * N, sizeof(double));
	if (result == NULL) {
		free(B_copy);
		exit(1);
	}

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N,
				1, B_copy, N, B, N, 1, result, N);

	// Compute A' * A. The result is stored in A.
	cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasTrans,
				CblasNonUnit, N, N, 1, A, N, A, N);

	// Compute A' * A + A * B * B' (A + result). The result is stored in result.
	for (int i = 0; i < N * N; i++) {
		*(result + i) = *(A + i) + *(result + i);
	}

	free(B_copy);

	return result;
}
