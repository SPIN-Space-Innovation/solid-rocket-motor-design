#pragma once
#include <iostream>
#include "Coef.hpp"
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>

enum ExitCase{sub, sup}; //flow exit cases: subsonic, supersonic

class OdEngine{
    public:
    std::vector<double> time_; // time array
    std::vector<double> p_c_; // chamber pressure array
    std::vector<double> R_; // grain inner radius array 
    std::vector<double> m_; // mass flow rate array
    std::vector<double> p_e_; // exit pressure  array
    std::vector<double> u_e_; // exit velocity array
    double A_e_; // exit cross section area
    double A_t_; // throat cross section area
    double MinThickness_; // casing minimum thickness
    double CharacteristicVelocity_; // characteristic velocity ( c* )
    double p_max_; // maximum pressure

// constructor 
OdEngine(double A_throat_, double A_exit_, double Dt_, std::string datfile_,double SafetyFactor_)
:A_e_(A_exit_), A_t_(A_throat_) 
{
    //function to solve the system of the 2 ode, to find the functions p_c = p_c(t), R = R(t)
    RungeKutta4(Dt_,datfile_);
    // function to calculate the minimum thickness of the casing taking as input the safety factor
    MinThickness_ = MinThickness(SafetyFactor_);
    //function to calculate the characteristic velocity (c*) 
    CharacteristicVelocity_ = CharacteristicVelocity(Dt_);
    //maximum chamber pressure 
    p_max_ = *max_element(p_c_.begin(), p_c_.end());


}

OdEngine() = default;
OdEngine(const OdEngine &other) = default;
~OdEngine() = default;


// function to calculate the chamber temperature taking as input the time moment
// the temperature is assumed to increase exponentially reaching the theoretical combustion chamber at time moment, t = tc (tc is defined by the user in the Coef.hpp header file)
// for t >= tc, the temperature remains constant
double T(double time){
    if (time<tc)
    {
        return (273.15+20)*exp(log(KNDX_PROPELLANT[CombustionTemperature]/(273.15+20))*time/tc);
    }
    if(time>=tc){
        return KNDX_PROPELLANT[CombustionTemperature];
    }
    
}

// function to calculate the total burning area taking as input  the grain inner radius, r
double BurningArea(double r){
    return GRAIN[NumberOfGrains]*(2*pi*r*(GRAIN[Length] - 2*BurningCaseCoef*(r - GRAIN[InnerDiametre]/2)) + 2*BurningCaseCoef*(pi*pow(GRAIN[OutDiametre],2)/4 - pi*r*r));
}
// function to calculate the density of the gas inside the chamber, using the perfect gas law and taking as inputs the time moment, time and the pressure
double density_0(double p_0, double time){
    return p_0/(KNDX_PROPELLANT[GasConstant] * T(time));
}
// function to calculate the empty volume inside the chamber taking as input the grain inner radius, r
double EmptyVolume(double r){
    return GRAIN[NumberOfGrains]*(  pi*r*r*(GRAIN[Length] - 2*(r - GRAIN[InnerDiametre]/2)) + 2*(r - GRAIN[InnerDiametre]/2)*pi*(pow(GRAIN[OutDiametre],2)/4 - r*r)  ) + GRAIN[NumberOfGaps]*pi*GRAIN[Gap]*pow(GRAIN[OutDiametre],2)/4;
}
// gamma constant of the exhaust gasses
double gamma = KNDX_PROPELLANT[SpecificHeatRatio];

// function to calculate the value of the mathematical function F(M,ratio) =  AMR(M,ratio) - ratio, as it is defined in the relation 4.9
double F(double M, double RATIO){
    return (1/M)*pow((2/(gamma+1))*(1+((gamma-1)/2)*M*M),(gamma+1)/(2*(gamma-1)))-RATIO;
}
// derivative of function F with respect of the Mach number, M
double dF(double M){
    double a1 = pow(2,1+(gamma+1)/(2*(gamma-1)));
    double a2 = pow((0.5*(gamma-1)*M*M+1)/(gamma+1),(gamma+1)/(2*(gamma-1)));
    double a3 = M*M*((gamma-1)*M*M+2);
    return a1*(M*M-1)*a2/a3;
}
    
    // area mach number relation, solved for the mach number using newton-raphson method, M_new = M_old - F(M_old)/[dF/dM(M_old)]
    // CASE = sup or sub
double AMR(int CASE){
        double ratio = A_e_/A_t_;
        double gamma = KNDX_PROPELLANT[SpecificHeatRatio];
        double err = 1;
        double M_0;
        switch (CASE)
        {
            case sub:
                M_0 = 0.1; // initialize M for the subsonic case
                while(fabs(err)>1e-3){
                    err = -F(M_0,ratio)/dF(M_0);
                    M_0 -= F(M_0,ratio)/dF(M_0);;
                }
                return M_0;
                break;
            case sup:
                M_0 = 2; // initialize M for the supersonic case
                while(fabs(err)>1e-3){
                    err = -F(M_0,ratio)/dF(M_0);
                    M_0 -= F(M_0,ratio)/dF(M_0);;
                }
                return M_0;
                break;
            default:
                break;
        }
        
    }


// a, b, c, m_n, d function used to calculate f1 function
// the system of the 2 odes is written in the form dp/dt = f1(r, p, t), dr/dt = f2(r, p, t)
double a(double r, double p_0){
    return KNDX_PROPELLANT[BurnRateCoef]   *   pow(p_0   ,   KNDX_PROPELLANT[BurnRateExponent]);
}
double b(double r, double p_0, double time){
    return (KNDX_PROPELLANT[Density]   -   density_0(p_0,time));
}
double m_n(double r, double p_0, double time){
    //double p_a = 1e5;
    //double M_e = sqrt((2/(gamma-1))*(pow(p_0/p_a, (gamma-1)/gamma) - 1));
    //if (M_e < AMR(sub))
    //{
    //   return p_0 * A_e_ * sqrt((2*gamma/(gamma-1)) * (1/(KNDX_PROPELLANT[GasConstant]*KNDX_PROPELLANT[CombustionTemperature])) * pow(p_0/p_a, 2/gamma) * pow(1-p_a/p_0, (gamma-1)/gamma));
    //
    //}
    //else{
        return p_0   *   A_t_   *   sqrt(gamma   *   pow(2/(gamma+1),(gamma+1)/(2*(gamma-1)))   /   (KNDX_PROPELLANT[GasConstant]   *   T(time)));
    //}
}
double d(double r, double p_0, double time){
    return (KNDX_PROPELLANT[GasConstant]   *   T(time))   /   EmptyVolume(r);
}

// function to calculate the value of f1 
double f1(double r, double p_0, double time){
    return (BurningArea(r) * a(r,p_0) * b(r,p_0,time) - m_n(r,p_0,time)) * d(r,p_0,time);
}

// function to calculate the value of f2
double f2(double r, double p_0, double time){
    return KNDX_PROPELLANT[BurnRateCoef]   *   pow(p_0   ,   KNDX_PROPELLANT[BurnRateExponent]);
}


// RUNGE KUTTA 4 ORDER K FUNCTIONS
double k1_f1(double r, double p_0, double dt, double time){
    return dt*f1(r,p_0,time);
}
double k1_f2(double r, double p_0, double dt, double time){
    return dt*f2(r,p_0,time);
}

double k2_f1(double r, double p_0, double dt, double time){
    return dt*f1(r + k1_f2(r,p_0,dt,time)/2, p_0 + k1_f1(r,p_0,dt,time)/2, time + dt/2);
}
double k2_f2(double r, double p_0, double dt, double time){
    return dt*f2(r + k1_f2(r,p_0,dt,time)/2, p_0 + k1_f1(r,p_0,dt,time)/2, time + dt/2);
}

double k3_f1(double r, double p_0, double dt, double time){
    return dt*f1(r + k2_f2(r,p_0,dt,time)/2, p_0 + k2_f1(r,p_0,dt,time)/2, time + dt/2);
}
double k3_f2(double r, double p_0, double dt, double time){
    return dt*f2(r + k2_f2(r,p_0,dt,time)/2, p_0 + k2_f1(r,p_0,dt,time)/2, time + dt/2);
}

double k4_f1(double r, double p_0, double dt, double time){
    return dt*f1(r + k3_f2(r,p_0,dt,time), p_0 + k3_f1(r,p_0,dt,time), time + dt);
}
double k4_f2(double r, double p_0, double dt, double time){
    return dt*f2(r + k3_f2(r,p_0,dt,time), p_0 + k3_f1(r,p_0,dt,time), time + dt);
}

void RungeKutta4(double Dt, std::string datfile){
    
    time_.push_back(0);
    R_.push_back(GRAIN[InnerDiametre]/2);
    p_c_.push_back(AmbientPressure);

    while (R_.back() < GRAIN[OutDiametre]/2)
    {
        double new_p = p_c_.back() + (k1_f1(R_.back(),p_c_.back(),Dt,time_.back()) + 2*k2_f1(R_.back(),p_c_.back(),Dt,time_.back()) + 2*k3_f1(R_.back(),p_c_.back(),Dt,time_.back()) + k4_f1(R_.back(),p_c_.back(),Dt,time_.back()))/6;
        double new_r = R_.back() + (k1_f2(R_.back(),p_c_.back(),Dt,time_.back()) + 2*k2_f2(R_.back(),p_c_.back(),Dt,time_.back()) + 2*k3_f2(R_.back(),p_c_.back(),Dt,time_.back()) + k4_f2(R_.back(),p_c_.back(),Dt,time_.back()))/6;
        R_.push_back(new_r);
        p_c_.push_back(new_p);
        time_.push_back(time_.back() + Dt);

        //std::cout<<t.back()<<"    "<<R.back()<<"    "<<P_0.back()/1e5<<"\n";
    }

    std::ofstream Results;
        Results.open(datfile);
        if (Results.is_open())
        {
            Results<<"# time    pressure    mass flow rate"<<"\n";
            for(int i = 0; i < time_.size(); ++i){
                Results<<time_.at(i)<<"         "<<p_c_.at(i)/1e5<<"        "<<m_n(0,p_c_.at(i),time_.at(i))<<"\n";
            }
            Results.close();
        }
}

double MinThickness(double SafetyFactor){
    double Pmax = *max_element(p_c_.begin(), p_c_.end());
    return 1e3*SafetyFactor*Pmax*GRAIN[OutDiametre]/(2*Sy); //mm
}

double CharacteristicVelocity(double Dt){
    double Mp = GRAIN[NumberOfGrains]*(pi*GRAIN[OutDiametre]*GRAIN[OutDiametre]/4 - pi*GRAIN[InnerDiametre]*GRAIN[InnerDiametre]/4)*GRAIN[Length]*KNDX_PROPELLANT[Density];
    double I = 0;
    for (int i = 0; i < time_.size()-1; i++)
    {
        I += Dt*(p_c_.at(i+1) + p_c_.at(i))/2;
    }
    return I*A_t_/Mp;
}


};