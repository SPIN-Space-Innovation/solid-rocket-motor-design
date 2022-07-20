#pragma once
#include <iostream>
#include "Coef.hpp"
#include <cmath>
#include <vector>
#include <fstream>

double pi = 2*acos(0);

double BurningArea(double r){
    return GRAIN[NumberOfGrains]*(2*pi*r*(GRAIN[Length] - 2*BurningCaseCoef*(r - GRAIN[InnerDiametre]/2)) + 2*BurningCaseCoef*(pi*pow(GRAIN[OutDiametre],2)/4 - pi*r*r));
}
double density_0(double p_0){
    return p_0/(KNDX_PROPELLANT[GasConstant] * KNDX_PROPELLANT[CombustionTemperature]);
}
double EmptyVolume(double r){
    return pi*(GRAIN[Length] - 2*BurningCaseCoef*(r - GRAIN[InnerDiametre]/2))*r*r + 2*BurningCaseCoef*(r - GRAIN[InnerDiametre]/2)*(pi*pow(GRAIN[OutDiametre],2)/4 - pi*r*r) + GRAIN[NumberOfGaps]*pi*GRAIN[Gap]*GRAIN[OutDiametre];
}

double gamma = KNDX_PROPELLANT[SpecificHeatRatio];

double a(double r, double p_0, double A_t){
    return KNDX_PROPELLANT[BurnRateCoef]   *   pow(p_0   ,   KNDX_PROPELLANT[BurnRateExponent]);
}
double b(double r, double p_0, double A_t){
    return (KNDX_PROPELLANT[Density]   -   density_0(p_0));
}
double c(double r, double p_0, double A_t){
    return p_0   *   A_t   *   sqrt(gamma   *   pow(2/(gamma+1),(gamma+1)/(2*(gamma-1)))   /   (KNDX_PROPELLANT[GasConstant]   *   KNDX_PROPELLANT[CombustionTemperature]));
}
double d(double r, double p_0, double A_t){
    return (KNDX_PROPELLANT[GasConstant]   *   KNDX_PROPELLANT[CombustionTemperature])   /   EmptyVolume(r);
}
double f1(double r, double p_0, double A_t){
    return (BurningArea(r) * a(r,p_0,A_t) * b(r,p_0,A_t) - c(r,p_0,A_t)) * d(r,p_0,A_t);
}

double f2(double r, double p_0, double A_t){
    return KNDX_PROPELLANT[BurnRateCoef]   *   pow(p_0   ,   KNDX_PROPELLANT[BurnRateExponent]);
}

// RUNGE KUTTA 4 ORDER K FUNCTIONS
double k1_f1(double r, double p_0, double dt, double A_t){
    return f1(r,p_0,A_t);
}
double k1_f2(double r, double p_0, double dt, double A_t){
    return f2(r,p_0,A_t);
}

double k2_f1(double r, double p_0, double dt, double A_t){
    return f1(r + 0.5*dt*k1_f2(r,p_0,dt,A_t),p_0 + 0.5*dt*k1_f1(r,p_0,dt,A_t),A_t);
}
double k2_f2(double r, double p_0, double dt, double A_t){
    return f2(r + 0.5*dt*k1_f2(r,p_0,dt,A_t),p_0 + 0.5*dt*k1_f1(r,p_0,dt,A_t),A_t);
}

double k3_f1(double r, double p_0, double dt, double A_t){
    return f1(r + 0.5*dt*k2_f2(r,p_0,dt,A_t),p_0 + 0.5*dt*k2_f1(r,p_0,dt,A_t),A_t);
}
double k3_f2(double r, double p_0, double dt, double A_t){
    return f2(r + 0.5*dt*k2_f2(r,p_0,dt,A_t),p_0 + 0.5*dt*k2_f1(r,p_0,dt,A_t),A_t);
}

double k4_f1(double r, double p_0, double dt, double A_t){
    return f1(r + dt*k3_f2(r,p_0,dt,A_t),p_0 + dt*k3_f1(r,p_0,dt,A_t),A_t);
}
double k4_f2(double r, double p_0, double dt, double A_t){
    return f2(r + dt*k3_f2(r,p_0,dt,A_t),p_0 + dt*k3_f1(r,p_0,dt,A_t),A_t);
}

std::vector <double> t;
std::vector <double> R;
std::vector <double> P_0;

void RungeKutta4(double A_throat, double Dt, std::string datfile){
    
    t.push_back(0);
    R.push_back(GRAIN[InnerDiametre]/2);
    P_0.push_back(AmbientPressure);

    while (R.back() < GRAIN[OutDiametre]/2)
    {
        double new_p = P_0.back() + Dt*(k1_f1(R.back(),P_0.back(),Dt,A_throat) + 2*k2_f1(R.back(),P_0.back(),Dt,A_throat) + 2*k3_f1(R.back(),P_0.back(),Dt,A_throat) + k4_f1(R.back(),P_0.back(),Dt,A_throat))/6;
        double new_r = R.back() + Dt*(k1_f2(R.back(),P_0.back(),Dt,A_throat) + 2*k2_f2(R.back(),P_0.back(),Dt,A_throat) + 2*k3_f2(R.back(),P_0.back(),Dt,A_throat) + k4_f2(R.back(),P_0.back(),Dt,A_throat))/6;
        R.push_back(new_r);
        P_0.push_back(new_p);
        t.push_back(t.back() + Dt);

        //std::cout<<t.back()<<"    "<<R.back()<<"    "<<P_0.back()/1e5<<"\n";
    }

    std::ofstream Results;
        Results.open(datfile);
        if (Results.is_open())
        {
            for(int i = 0; i < t.size(); ++i){
                Results<<t.at(i)<<"   "<<P_0.at(i)/1e5<<"\n";
            }
            Results.close();
        }
}
