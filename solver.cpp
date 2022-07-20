#include <iostream>
#include "ChamberPressureSolver.hpp"


int main(){
    

    RungeKutta4(0.25*pi*pow(15.56e-3,2),0.0001,"results4.dat");


    return 0;

}