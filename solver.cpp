#include <iostream>
#include "OdEngine.hpp"
//#include "Nozzle.hpp"
#include <algorithm>
#include <cmath>

int main(){
    
    OdEngine dummy(A_t,A_e,0.005,"ChamberPressure.dat", 2.5);
    std::cout<<"c* = "<<dummy.CharacteristicVelocity_<<"\n";
    std::cout<<"t_min = "<<dummy.MinThickness_<<"\n";
    std::cout<<"p_max = "<<dummy.p_max_/1e5<<"\n";
    std::cout<<"M_prop = "<<M_prop<<"\n";
    


    return 0;

}
