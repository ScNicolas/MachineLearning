#include <iostream>
#include <cstdlib>     /* srand, rand */
#include <ctime>
#include "Eigen/Dense"
#include <math.h>

extern "C"
{
	struct rbfStruct
	{
		int nbPoints;
		double ** trainingPoints;
		double * weight;
		int alpha;
	};



	__declspec(dllexport)rbfStruct * createModel(double * points[], int nbPoints, int sizePoints , int alpha) 
	{
		rbfStruct * rbf = (rbfStruct*)malloc(sizeof(rbfStruct));
		rbf->nbPoints = nbPoints;
		rbf->trainingPoints = points;
		rbf->weight = (double *)malloc(sizeof(double) * rbf->nbPoints);
		rbf->alpha = alpha;

		return rbf;
	}

	double get_distance(double x [], double y [] , int size) {
		
		double tmpValue = 0.0;
		for (int i = 1; i < size; i++) {
			tmpValue += pow(x[i] - y[i],2);
		}
		return tmpValue;

	}

	__declspec(dllexport) double * train_rbf(double  points[], int  size, int nbPoints, double alpha) {

		Eigen::MatrixXd x(nbPoints, nbPoints);
		Eigen::MatrixXd y(nbPoints, 1);

		for (int i = 0; i < nbPoints; i++)
		{
			for (int j = 0; j < nbPoints; j++)
			{
				double * pointj = (double *)malloc(sizeof(double) * size);
				for (int k = 0; k<size; k++) {
					pointj[k] = points[j*size + k];
				}
				double * pointi = (double *)malloc(sizeof(double) * size);
				for (int k = 0; k<size; k++) {
					pointi[k] = points[i*size + k];
				}

				x(i, j) = exp(alpha * -1 * get_distance(pointi, pointj, size));
				free(pointi);
				free(pointj);
			}
		}
		for (int i = 0; i < nbPoints; i++)
		{
			y(i, 0) = points[i*size];
		}

		Eigen::MatrixXd tmp = x.inverse();
		Eigen::MatrixXd wtmp = tmp*y;

		double * w = (double*)malloc(sizeof(double)*nbPoints);
		for (int i = 0; i < nbPoints; i++) {
			w[i] = wtmp(i, 0);
		}

		return w;
	}

	__declspec(dllexport) double execRbfRegression(double w[], double x[], double * trainingPoints, int nbTrainingsPoints, int size, double alpha) {

		double sum = 0;

		for (int i = 0; i < nbTrainingsPoints; i++)
		{
			double * point = (double *)malloc(sizeof(double) * size);
			for (int j = 0; j<size; j++) {
				point[j] = trainingPoints[i*size + j];
			}
			sum += w[i] * exp(alpha * -1 * get_distance(x, point, size));
			free(point);
		}
		return sum;
	}

	__declspec(dllexport) double execRbfClassif(double w[], double x[], double * trainingPoints, int nbTrainingsPoints, int size, double alpha) {

		double sum = 0;

		for (int i = 0; i < nbTrainingsPoints; i++)
		{
			double * point = (double *)malloc(sizeof(double) * size);
			for (int j = 0; j<size; j++) {
				point[j] = trainingPoints[i*size + j];
			}
			sum += w[i] * exp(alpha * -1 * get_distance(x, point, size));
			free(point);
		}
		if (sum > 0)
			return 1;
		else if (sum < 0)
			return -1;
		else return 0;
	}
	
	int sign(double w[], double x[], int size) {

		double sum = 0;
		for (int i = 1; i < size; i++) {
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
		srand(time(NULL));
		for (int i = 0; i < size; i++) {
			w[i] = rand() % maxRange + minRange;
		}
		return w;
	}

	double *generateRosenBlatt(int size) {
		return generateWeight(size, -1, 2);
	}

	double *executeRosenBlatt(double w[], double point[], int size) {
		double alpha = 0.1;
		int signPoint = sign(w, point, size);
		while (signPoint != point[0]) {
			for (int i = 1; i < size; i++) {
				w[i] = w[i] + (alpha * (point[0] - signPoint)) * point[i];
			}
			w[0] = w[0] + (alpha * (point[0] - signPoint));
			int signPoint = sign(w, point, size);
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
		srand(time(NULL));
		for (int i = 0; i < size; i++) {
			double f = (double)rand() / RAND_MAX;
			w[i] = i; //-1.0 + f * (2.0);
			//w[i] = (double)rand() % 2 + -1.0;
		}
		return w;
	}

	__declspec(dllexport) double *trainLinear(double  w[], int size, double points[], int nbrPoints) {

		int cpt = 0;
		bool find = true;
		while (cpt < 100000 && find == true) {
			find = false;
			for (int i = 0; i < nbrPoints; i++) {
				double * point = (double *)malloc(sizeof(double) * size);
				for (int j = 0; j < size; j++) {
					point[j] = points[(i * size) + j];
				}
				double alpha = 0.1;
				int signPoint = sign(w, point, size);		
				while (signPoint != point[0]) {
					find = true;
					for (int i = 1; i < size; i++) {
						w[i] = w[i] + (alpha * (point[0] - signPoint) * point[i]);
					}
					w[0] = w[0] + (alpha * (point[0] - signPoint));
					signPoint = sign(w, point, size);
				}
				free(point);
			}
			cpt++;
		}
		return w;
	}

	__declspec(dllexport) double* executeLinear(double w[], int size, double points[], int nbr_points) {
		for (int i = 0; i < nbr_points; i++) {
			double * point = (double *)malloc(sizeof(double) * size);

			for (int j = 0; j < size; j++) {
				point[j] = points[(i * size) + j];
			}
			points[(i * size)] = sign(w, point, size);
		}
		return points;

	}

	__declspec(dllexport) double * train_linear_regression(double points[], int nbTrainingSphere, int sizeW) {

		Eigen::MatrixXd x(nbTrainingSphere, (sizeW));
		Eigen::MatrixXd y(nbTrainingSphere, 1);

		for (int i = 0; i < nbTrainingSphere; i++) {
			x(i, 0) = 1;
			y(i, 0) = points[i * (sizeW )];

			for (int j = 1; j<(sizeW); j++) {
				x(i, j) = points[i * (sizeW) + j];
			}
		}

		Eigen::MatrixXd xT = x.transpose();
		Eigen::MatrixXd xTx = xT * x;
		Eigen::MatrixXd xTxInv = xTx.inverse();
		Eigen::MatrixXd xTxInvXT = xTxInv * xT;
		Eigen::MatrixXd matrixW = xTxInvXT * y;

		double * w = (double*)malloc(sizeof(double)*sizeW);
		for (int i = 0; i < sizeW; i++) {
			w[i] = matrixW(i, 0);
		}
		return w;
	}


	__declspec(dllexport) double execRegression(double w[], double x[], int size) {

		double sum = 0;
		for (int i = 1; i < size; i++) {
			sum += w[i] * x[i];
		}
		sum += w[0];
		
		return sum;
	}


	__declspec(dllexport) int evaluateLinear(double w[], int size, double point[]) {

		return sign(w, point, size);
	}


	__declspec(dllexport) int hello() {
		return 10;
	}



}
