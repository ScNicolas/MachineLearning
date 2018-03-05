#include "library.h"
#define DllExport   __declspec( dllexport )
#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <time.h>

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

    double generateWeight(int size, int minRange, int maxRange)[]{
        double w[size];
        for(int i=0;i<size;i++){
            srand (time(NULL));
            w[i]= rand() % maxRange + minRange;
        }
        return w;
    }




};