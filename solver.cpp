#include <iostream>
#include "CombustionChamber.hpp"
#include "Nozzle.hpp"
#include <algorithm>
#include <cmath>

int main(){
    

    // throat diametre
    double d_t = 15.56e-3;
    // throat cross section area
    double A_t = 0.25*pi*pow(d_t,2);

    CombustionChamber dummyCC(A_t,0.0001,"ChamberPressure.dat",2.5);
    std::cout<<dummyCC.CharacteristicVelocity_<<"\n";
    std::cout<<dummyCC.MinThickness_<<"\n";
    std::cout<<dummyCC.p_max_/1e5<<"\n";

    
    Nozzle dummyNozzle("ChamberPressure.dat",0.25*pi*pow(15.56e-3,2),8, "exit_pressure_sup.dat", "exit_pressure_NSE.dat", "thrust.dat");  
    dummyNozzle.print_massFlowRate("massFlowRate.dat");
    dummyNozzle.print_exitTemperature("exitTemperature.dat");
    std::cout<<dummyNozzle.Impulse_<<"\n";

    return 0;

}
