#include "utils.h"

double* my_solver(int N, double *A, double* B) {
	double *At, *Bt, *BBt, *first_term, *second_term, *result;
	register double *pa, *pat, *orig_pa, *pb, *pc, sum;
	int i, j, k, bi, bj, bk;
	int blockSize = 40;
	
	printf("OPT SOLVER\n");

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
		pa = B + i;
		pat = Bt + i * N;
		for (j = 0; j < N; j++) {
			*pat = *pa;
			pa += N;
			pat++;
		}
	}

	// Compute B * B' (B * Bt). The result is stored in BBt.
	for (bi = 0; bi < N; bi += blockSize) {
        for (bj = 0; bj < N; bj += blockSize) {
            for (bk = 0; bk < N; bk += blockSize) {
                for (i = 0; i < blockSize; i++) {
                    orig_pa = B + (bi + i) * N + bk;
					pc = BBt + (bi + i) * N + bj;
					for (j = 0; j < blockSize; j++) {
                        pa = orig_pa;
						pb = Bt + bk * N + bj + j;
						sum = 0.0;
						for (k = 0; k < blockSize; k++) {
                            sum += *pa * *pb;
							pa++;
							pb += N;
						}
						*pc += sum;
						pc++;
					}
				}
			}
		}
	}

	// Compute A * B * B' (A * BBt). The result is stored in first_term.
	for (bi = 0; bi < N; bi += blockSize) {
        for (bj = 0; bj < N; bj += blockSize) {
            for (bk = bi; bk < N; bk += blockSize) {
                for (i = 0; i < blockSize; i++) {
                    orig_pa = A + (bi + i) * N + bk;
					pc = first_term + (bi + i) * N + bj;
					for (j = 0; j < blockSize; j++) {
                        pa = orig_pa;
						pb = BBt + bk * N + bj + j;
						sum = 0.0;
						for (k = 0; k < blockSize; k++) {
                            sum += *pa * *pb;
							pa++;
							pb += N;
						}
						*pc += sum;
						pc++;
					}
				}
			}
		}
	}

	// Compute A'. The result is stored in At.
	for (i = 0; i < N; i++) {
		pa = A + i;
		pat = At + i * N;
		for (j = 0; j <= i; j++) {
			*pat = *pa;
			pa += N;
			pat++;
		}
	}

	// Compute A' * A (At * A). The result is stored in second_term.
	for (bi = 0; bi < N; bi += blockSize) {
        for (bj = 0; bj < N; bj += blockSize) {
            for (bk = bj; bk >= 0; bk -= blockSize) {
                for (i = 0; i < blockSize; i++) {
                    orig_pa = At + (bi + i) * N + bk;
					pc = second_term + (bi + i) * N + bj;
					for (j = 0; j < blockSize; j++) {
                        pa = orig_pa;
						pb = A + bk * N + bj + j;
						sum = 0.0;
						for (k = 0; k < blockSize; k++) {
                            sum += *pa * *pb;
							pa++;
							pb += N;
						}
						*pc += sum;
						pc++;
					}
				}
			}
		}
	}

	result = (double*) calloc(N * N, sizeof(double));

	// Compute A' * A + A * B * B' (first_term + second_term).
	// The result is stored in result.
	for (i = 0; i < N * N; i++) {
		*(result + i) = *(first_term + i) + *(second_term + i);
	}

	free(At);
	free(Bt);
	free(BBt);
	free(first_term);
	free(second_term);

	return result;
}
