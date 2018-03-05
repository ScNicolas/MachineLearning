#include "library.h"
#define DllExport   __declspec( dllexport )
#include <iostream>

class DllExport C
{
     void hello() {
        std::cout << "Hello, World!" << std::endl;
    }
};