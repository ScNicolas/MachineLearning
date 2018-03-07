
#define DllExport   __declspec( dllexport )
#include <iostream>
#include <cstdlib>     /* srand, rand */
#include <ctime>

class DllExport C
{
    void hello() {
        std::cout << "Hello, World!" << std::endl;
    }

    int sign(double w[], double x[]){

        double sum=0;
        for(int i=1;i<=2;i++){
            sum+=w[i]*x[i];
        }
        sum+=w[0];

        if(sum>0)
            return 1;
        else if(sum <0)
            return -1;
        return 0;
    }

    double* generateWeight(int size, int minRange, int maxRange){
        double w[size];
        for(int i=0;i<size;i++){
            srand (time(NULL));
            w[i]= rand() % maxRange + minRange;
        }
        return w;
    }

    double* generateRosenBlatt(int size){
        return generateWeight(size,-1,1);
    }

    double* executeRosenBlatt(double w[], double point[], int size) {
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

    double* executeLinear(int size,double points[], int nbrPoints) {
        double* w=generateRosenBlatt(size);
        for(int i=0;i<nbrPoints;i++) {
            double point[size];
            for(int j=0;j<size;j++){
                point[j]=points[(i*size)+j];
            }
            executeRosenBlatt(w,point,size);
        }

        return w;
    }
};