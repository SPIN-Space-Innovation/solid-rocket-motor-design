#include <iostream>
#include "ChamberPressureSolver.hpp"
#include <fstream>

int main(){

    std::vector <double> Pressure;
    std::vector <double> radius;
    std::vector <double> time;

    RungeKutta4(0.25*pi*pow(15.56e-3,2),0.0001,time, radius, Pressure,"results.dat");

    return 0;

}