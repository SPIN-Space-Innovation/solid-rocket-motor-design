#include <iostream>
#include "CombustionChamber.hpp"
#include "Nozzle.hpp"
#include <algorithm>
#include <cmath>

int main(){
    

    CombustionChamber CC(0.25*pi*pow(15.56e-3,2),0.0001,"ChamberPressure.dat",2.5);
    
    //Nozzle dummyNozzle("ChamberPressure.dat",0.25*pi*pow(15.56e-3,2),4, "exit_pressure_sup.dat", "exit_pressure_NSE.dat", "thrust.dat");  
   
    


    
    

    return 0;

}