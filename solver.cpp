#include <iostream>
#include "OdEngineNEW.hpp"
#include <algorithm>
#include <cmath>

// g++ -o solver.exe -I C:\Users\VS0121\eigen-3.4.0 solver.cpp

int main(){
    
    OdEngine dummy(A_t,A_e,0.001,"ChamberPressure.dat", 2.5);
    dummy.ExitConditions();
    dummy.WrightExitPressureToFile("exitPressure.dat");
    dummy.WrightThrustToFile("Thrust.dat");
    dummy.WrightExitVelocityToFile("Velocity.dat");
    dummy.WrightMassFlowRateToFile("MassFlowRate.dat");
    std::cout<<"c* = "<<dummy.CharacteristicVelocity_<<"\n";
    std::cout<<"t_min = "<<dummy.MinThickness_<<"\n";
    std::cout<<"p_max = "<<dummy.p_max_/1e5<<"\n";
    std::cout<<"M_prop = "<<M_prop<<"\n";
    std::cout<<"time interval = "<<dummy.time_(dummy.time_.size()-1)<<"\n";
    std::cout<<"I_t = "<<dummy.TotalImpulse_<<"\n";


    return 0;

}
