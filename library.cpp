#include "library.h"
#define DllExport   __declspec( dllexport )
#include <iostream>

class DllExport C
{
    void hello() {
        std::cout << "Hello, World!" << std::endl;
    }

    int sign(double w[], double x[]){
        double sum=0;
        for(int i=0;i<2;i++){
            sum+=w[i]*x[i];
        }
        sum+=w[2];

        if(sum>0)
            return 1;
        else if(sum <0)
            return -1;
        return 0;

    }
};