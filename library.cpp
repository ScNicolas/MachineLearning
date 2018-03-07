
#include <iostream>
#include <cstdlib>     /* srand, rand */
#include <ctime>

extern "C"
{
	int sign(double w[], double x[]) {

		double sum = 0;
		for (int i = 1; i <= 2; i++) {
			sum += w[i] * x[i];
		}
		sum += w[0];

		if (sum > 0)
			return 1;
		else if (sum < 0)
			return -1;
		return 0;
	}

	double *generateWeight(int size, int minRange, int maxRange) {
		double *w = (double *)malloc(sizeof(double) * size);
		for (int i = 0; i < size; i++) {
			srand(time(NULL));
			w[i] = rand() % maxRange + minRange;
		}
		return w;
	}

	double *generateRosenBlatt(int size) {
		return generateWeight(size, -1, 1);
	}

	double *executeRosenBlatt(double w[], double point[], int size) {
		double alpha = 0.1;
		int signPoint = sign(w, point);
		while (signPoint != point[0]) {
			for (int i = 1; i < size; i++) {
				w[i] = w[i] + (alpha * (point[0] - signPoint)) * point[i];
			}
			w[0] = w[0] + (alpha * (point[0] - signPoint));
		}
		return w;
	}

	double *executeLinear2(int size, double points[], int nbrPoints) {
		double *w = generateRosenBlatt(size);
		for (int i = 0; i < nbrPoints; i++) {
			double * point = (double *)malloc(sizeof(double) * size);
			for (int j = 0; j < size; j++) {
				point[j] = points[(i * size) + j];
			}
			executeRosenBlatt(w, point, size);
		}

		return w;
	}
	__declspec(dllexport) double *generateModel(int size) {
		double *w = (double *)malloc(sizeof(double) * size);
		for (int i = 0; i < size; i++) {
			srand(time(NULL));
			w[i] = rand() % 1 + -1;
		}
		return w;
	}

	__declspec(dllexport) double *trainLinear(double  w[], int size, double points[], int nbrPoints) {

		for (int i = 0; i < nbrPoints; i++) {
			double * point = (double *)malloc(sizeof(double) * size);
			for (int j = 0; j < size; j++) {
				point[j] = points[(i * size) + j];
			}
			double alpha = 0.1;
			int signPoint = sign(w, point);
			while (signPoint != point[0]) {
				for (int i = 1; i < size; i++) {
					w[i] = w[i] + (alpha * (point[0] - signPoint) * point[i]);
				}
				w[0] = w[0] + (alpha * (point[0] - signPoint));
				signPoint = sign(w, point);
			}
            free(point);
		}
		return w;
	}

	__declspec(dllexport) double* executeLinear(double w[], int size, double points[], int nbr_points) {
		for (int i = 0; i<nbr_points; i++) {
			double * point = (double *)malloc(sizeof(double) * size);

			for (int j = 0; j < size; j++) {
				point[j] = points[(i * size) + j];
			}
			points[(i * size)] = sign(w, point);
            free(point);
		}

		return points;

	}

	__declspec(dllexport) int hello() {
		return 10;
	}
}
