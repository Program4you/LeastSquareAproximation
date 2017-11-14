#include <stdio.h>
#include <math.h>
#include <stdlib.h>

const double eps = 0.01;

double f(double x, double q) {
	return pow(cos(x), 1.0 / q);
}

// Gift :) (linear least square method)
void linear_extrapolation(double *x, double *y, int n, double *a, double *b) {
	double sum_xi = 0;
	double sum_yi = 0;
	double sum_xi_yi = 0;
	double sum_xi_xi = 0;


	for (int i = 0; i < n; i++) {
		sum_xi += x[i];
		sum_yi += y[i];
		sum_xi_yi += x[i] * y[i];
		sum_xi_xi += x[i] * x[i];
	}

	*a = (n * sum_xi_yi - sum_xi * sum_yi) / (n * sum_xi_xi - sum_xi * sum_xi);
	*b = (sum_yi - *a * sum_xi) / n;
}

double** ls_extrapolation(double *x, double *y, int n, int basis) {
	double **matrix = (double **) malloc(basis * sizeof(double *));

	for (int i = 0; i < basis; i++)
		matrix[i] = (double *) malloc((basis + 1) * sizeof(double));
	
	for (int i = 0; i < basis; i++) {
		for (int j = 0; j < basis; j++) {
			double sum_a = 0;

			for (int k = 0; k < n; k++)
				sum_a += pow(x[k], i + j);

			matrix[i][j] = sum_a;
			
		}
	}

	for (int i = 0; i < basis; i++) {
		double sum_b = 0;

		for (int k = 0; k < n; k++)
			sum_b += y[k] * pow(x[k], i);

		matrix[i][basis] = sum_b;
	}

	/*printf("Matrix basis:\n");
	for (int i = 0; i < basis; i++) {
		for (int j = 0; j <= basis; j++)
			printf("%10.5lf", matrix[i][j]);

		printf("\n");
	}*/

	return matrix;
}

void swap_rows(double **matrix, int n, int row1, int row2) {
	for (int i = 0; i < n; i++) {
		double tmp = matrix[row1][i];
		matrix[row1][i] = matrix[row2][i];
		matrix[row2][i] = tmp;
	}
}

double * gauss(double** matrix, int n) {
	double *x = (double *) malloc(n * sizeof(double));

	// direct
	for (int j = 0; j < n; j++) {
		if (matrix[j][j] == 0) {
			int k = j;

			while (k < n && !matrix[k][j])
				k++;

			if (k == n) {
				printf("Eror Gauss\n");
				return NULL;
			}				

			swap_rows(matrix, n + 1, j, k);
		}

		double v = matrix[j][j];

		for (int i = j + 1; i < n; i++) {
			v = -matrix[i][j] / v;

			for (int j1 = 0; j1 <= n; j1++)
				matrix[i][j1] += matrix[j][j1] * v;
		}
	}

	/*printf("\nMatrix basis after Gauss:\n");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j <= n; j++)
			printf("%10.5lf", matrix[i][j]);

		printf("\n");
	}*/

	for (int i = n - 1; i >= 0; i--) {
		double sum = matrix[i][n]; // bi

		for (int j = n - 1; j > i; j--)
			sum -= matrix[i][j] * x[j];

		x[i] = sum / matrix[i][i];
	}

	return x;
}

double aproximate(double x0, double *x, double *y, int n) {
	double s0;
	double *coef = NULL;

	for (int basis = 1; basis < 4; basis++) {
		double **matrix = ls_extrapolation(x, y, n + 1, basis);
		coef = gauss(matrix, basis);

		for (int i = 0; i < basis; i++)
			free(matrix[i]);

		free(matrix);

		double s = 0;

		for (int i = 0; i < n; i++) {
			double sum = 0;

			for (int j = 0; j < basis; j++)
				sum += coef[j] * pow(x[i], j);

			s += pow(y[i] - sum, 2);
		}

		s0 = 0;
		for (int i = 0; i < basis; i++)
			s0 += coef[i] * pow(x0, i);

		free(coef);
		
		if (s < eps)
			break;
	}	

	return s0;
}

int main() {
	const double left = 0;
	const double right = 1;

	int n;
	printf("Enter n: ");
	scanf("%d", &n);

	double q;
	printf("Enter q: ");
	scanf("%lf", &q);

	double h = (right - left) / n;

	double *x = (double *) malloc((n + 1) * sizeof(double));
	double *y = (double *) malloc((n + 1) * sizeof(double));

	for (int i = 0; i <= n; i++) {
		x[i] = left + i * h;
		y[i] = f(x[i], q);
	}

	printf("+============+============+============+=================+\n");
	printf("|     Xn     |    f(x)    |    s(x)    |  |f(x) - s(x)|  |\n");
	printf("+============+============+============+=================+\n");
	for (int i = 0; i <= n; i++) {
		double x0 = left + i * h;
		double f0 = f(x0, q);
		double s0 = aproximate(x0, x, y, n);

		printf("| %10.8lf | %10.8lf | %10.8lf | %15.13lf |\n", x0, f0, s0, fabs(f0 - s0));
	}
	printf("+============+============+============+=================+\n");
}