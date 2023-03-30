#include "readData.h"
#define print(x) std::cout<<#x<<" = "<<x<<std::endl;

int main(){

    readData("dummyCoef.txt");
    print(GRAIN[NumberOfGrains]);


    return 0;
}