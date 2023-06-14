#include "utils.h"

double* my_solver(int N, double *A, double* B) {
	double *At, *Bt, *BBt, *first_term, *second_term, *result;
	int i, j, k;
	
	printf("NEOPT SOLVER\n");

	// Allocate memory for At, Bt, BBt, first_term, second_term.
	At = (double*) calloc(N * N, sizeof(double));
	if (At == NULL) {
		exit(1);
	}
	Bt = (double*) calloc(N * N, sizeof(double));
	if (Bt == NULL) {
		free(At);
		exit(1);
	}
	BBt = (double*) calloc(N * N, sizeof(double));
	if (BBt == NULL) {
		free(At);
		free(Bt);
		exit(1);
	}
	first_term = (double*) calloc(N * N, sizeof(double));
	if (first_term == NULL) {
		free(At);
		free(Bt);
		free(BBt);
		exit(1);
	}
	second_term = (double*) calloc(N * N, sizeof(double));
	if (second_term == NULL) {
		free(At);
		free(Bt);
		free(BBt);
		free(first_term);
		exit(1);
	}
	
	// Compute B'. The result is stored in Bt.
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			*(Bt + j * N + i) = *(B + i * N + j);
		}
	}

	// Compute B * B' (B * Bt). The result is stored in BBt.
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			for (k = 0; k < N; k++) {
				*(BBt + i * N + j) += *(B + i * N + k) * *(Bt + k * N + j);
			}
		}
	}

	// Compute A * B * B' (A * BBt). The result is stored in first_term.
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			for (k = N - 1; k >= i; k--) {
				*(first_term + i * N + j) += *(A + i * N + k) * *(BBt + k * N + j);
			}
		}
	}

	// Compute A'. The result is stored in At.
	for (i = 0; i < N; i++) {
		for (j = 0; j <= i; j++) {
			*(At + i * N + j) = *(A + j * N + i);
		}
	}

	// Compute A' * A (At * A). The result is stored in second_term.
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			for (k = 0; k <= j; k++) {
				*(second_term + i * N + j) += *(At + i * N + k) * *(A + k * N + j);
			}
		}
	}

	result = (double*) calloc(N * N, sizeof(double));
	if (result == NULL) {
		exit(1);
	}

	// Compute A' * A + A * B * B' (first_term + second_term).
	// The result is stored in result.
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			*(result + i * N + j) = *(first_term + i * N + j) + *(second_term + i * N + j);
		}
	}

	free(At);
	free(Bt);
	free(BBt);
	free(first_term);
	free(second_term);

	return result;
}
