#include <iostream>
#include "ChamberPressureSolver.hpp"

int main(){
    

    RungeKutta4(0.25*pi*pow(15.56e-3,2),0.0001,"results.dat");
    std::cout<<CharacteristicVelocity(0.25*pi*pow(15.56e-3,2),0.0001)<<"\n";
    return 0;

}