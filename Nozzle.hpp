#pragma once
#include <fstream>
#include <string>
#include <cmath>
#include "Coef.hpp"


enum ExitCase{sub, sup};

class Nozzle{
    public:

    // All these vectors/variables are functions of time
    std::vector <double> t_; // time moment vector
    std::vector <double> p_c_; // chamber pressure vector
    std::vector <double> p_e_sup_; // ambient pressure for supersonic flow inside the divergent part of the nozzle
    std::vector <double> p_e_NSE_; // ambient pressure, for which normal shock at the exit of the nozzle occurs
    std::vector <double> m_dot_; // mass flow rate
    std::vector <double> exit_temperature_; //exit temperature
    std::vector <double> exit_velocity_; // exit flow velocity
    std::vector <double> thrust_; //thrust

    double A_t_; // throat cross area
    double ratio_; // expansion ratio 
    double M_e_sup_; // supersonic mach number at the exit of the nozzle ( M_sup(ratio, M>1) )
    double Impulse_; // total impulse

    // datfile stores the data for the combustion chamber
    // datfile_p_e_sup will store the data of the vector p_e_sup_
    // datfile_p_e_NSE will store the data of the vector p_e_NSE_
    // datfile_thrust will store the data of the vector thrust_
    Nozzle(std::string datfile, double A_t, double ratio, std::string datfile_p_e_sup, std::string datfile_p_e_NSE, std::string datfile_thrust)
    :A_t_(A_t), ratio_(ratio)
    {
        std::fstream readDat;
            readDat.open(datfile);
            if(readDat.is_open()){
                double t, p;
                while(readDat >> t >> p){
                    t_.push_back(t);
                    p_c_.push_back(p*1e5);
                }
            }
            readDat.close();
            M_e_sup_ = AMR(ratio_,sup);
            P_e_sup(M_e_sup_);
            P_e_NSE(M_e_sup_);
            M_dot();
            Exit_temperatute();
            Exit_velocity();
            Thrust();
            print_p_e_sup(datfile_p_e_sup);
            print_p_e_NSE(datfile_p_e_NSE);
            print_thrust(datfile_thrust);
            Impulse();
    }

    double F(double M, double RATIO){
        double gamma = KNDX_PROPELLANT[SpecificHeatRatio];
            return (1/M)*pow((2/(gamma+1))*(1+((gamma-1)/2)*M*M),(gamma+1)/(2*(gamma-1)))-RATIO;
    }
    double a1(double M){
        double gamma = KNDX_PROPELLANT[SpecificHeatRatio];
        return pow(2,1+(gamma+1)/(2*(gamma-1)));
    }
    double a2(double M){
        double gamma = KNDX_PROPELLANT[SpecificHeatRatio];
        return pow((0.5*(gamma-1)*M*M+1)/(gamma+1),(gamma+1)/(2*(gamma-1)));
    }
    double a3(double M){
        double gamma = KNDX_PROPELLANT[SpecificHeatRatio];
        return M*M*((gamma-1)*M*M+2);
    }
    double dF(double M){
        return a1(M)*(M*M-1)*a2(M)/a3(M);
    }
    
    // area mach number relation, solved for the mach number using newton-raphson method 
    // CASE = sup or sub
    double AMR(double ratio, int CASE){
        double gamma = KNDX_PROPELLANT[SpecificHeatRatio];
        double err = 1;
        double M_0;
        switch (CASE)
        {
            case sub:
                M_0 = 0.1;
                while(fabs(err)>1e-3){
                err = -F(M_0,ratio)/dF(M_0);
                //std::cout<<"ERROR = "<<err<<"   ";
                M_0 -= F(M_0,ratio)/dF(M_0);;
                //std::cout<<"Me = "<<M_0<<"\n";
                }
                return M_0;
                break;
            case sup:
                M_0 = 2;
                while(fabs(err)>1e-3){
                err = -F(M_0,ratio)/dF(M_0);
                //std::cout<<"ERROR = "<<err<<"   ";
                M_0 -= F(M_0,ratio)/dF(M_0);;
                //std::cout<<"Me = "<<M_0<<"\n";
                }
                return M_0;
                break;
            default:
                break;
        }
        
    }

    // temperature
    double T(double time){
    if (time<tc)
    {
        //return 273.15+20 + KNDX_PROPELLANT[CombustionTemperature]*time/tc;
        return (273.15+20)*exp(log(KNDX_PROPELLANT[CombustionTemperature]/(273.15+20))*time/tc);
    }
    if(time>=tc){
        return KNDX_PROPELLANT[CombustionTemperature];
    }
    
}

    void P_e_sup(double Me){
        double gamma = KNDX_PROPELLANT[SpecificHeatRatio];
        for (int i = 0; i < p_c_.size(); i++)
        {
            p_e_sup_.push_back((pow(1+((gamma-1)/2)*pow(Me,2),-gamma/(gamma-1)))*p_c_.at(i));
        }
        
    }
    
    void P_e_NSE(double Me){
        double gamma = KNDX_PROPELLANT[SpecificHeatRatio];
        double pb_pe = 1 + (2*gamma/(gamma+1))*(Me*Me-1);
        for (int i = 0; i < p_e_sup_.size(); i++)
        {
            p_e_NSE_.push_back(pb_pe*p_e_sup_.at(i));
        }
    }

    void M_dot(){
        double gamma = KNDX_PROPELLANT[SpecificHeatRatio];
        double G = sqrt(gamma)*pow(2/(gamma+1),(gamma+1)/(2*(gamma-1)));
        for (int i = 0; i < p_c_.size(); i++)
        {
            m_dot_.push_back(G*A_t_*p_c_.at(i)/sqrt(KNDX_PROPELLANT[GasConstant]*T(t_.at(i))));
        }
        
    }
    
    void Exit_temperatute(){
        double gamma = KNDX_PROPELLANT[SpecificHeatRatio];
        for (int i = 0; i < p_c_.size(); i++)
        {
            exit_temperature_.push_back(T(t_.at(i))*pow(1 + ((gamma-1)/2)*M_e_sup_*M_e_sup_,-1));
        }
        
    }

    void Exit_velocity(){
        double gamma = KNDX_PROPELLANT[SpecificHeatRatio];
        for (int i = 0; i < exit_temperature_.size(); i++)
        {
            exit_velocity_.push_back(M_e_sup_*sqrt(gamma*KNDX_PROPELLANT[GasConstant]*exit_temperature_.at(i)));
        }
        
    }

    void Thrust(){
        for (int i = 0; i < p_c_.size(); i++)
        {
            thrust_.push_back(m_dot_.at(i)*exit_velocity_.at(i) + (p_e_sup_.at(i) - 1e5)*ratio_*A_t_);
        }
        
    }

    void Impulse(){
        double I = 0;
        for(int i = 1; i < t_.size(); ++i){
            if (thrust_.at(i-1)>0)
            {
                I += (thrust_.at(i)+thrust_.at(i-1))*(t_.at(i)-t_.at(i-1))/2;
            }
        }
        Impulse_ = I;
    }

    //print functions

    void print_p_e_sup(std::string datfile){
        std::ofstream Results;
        Results.open(datfile);
        if (Results.is_open())
        {
            for(int i = 0; i < t_.size(); ++i){
                Results<<t_.at(i)<<"   "<<p_e_sup_.at(i)/1e5<<"\n";
            }
            Results.close();
        }
    }

    void print_p_e_NSE(std::string datfile){
        std::ofstream Results;
        Results.open(datfile);
        if (Results.is_open())
        {
            for(int i = 0; i < t_.size(); ++i){
                Results<<t_.at(i)<<"   "<<p_e_NSE_.at(i)/1e5<<"\n";
            }
            Results.close();
        }
    }

    void print_thrust(std::string datfile){
        std::ofstream Results;
        Results.open(datfile);
        if (Results.is_open())
        {
            for(int i = 0; i < t_.size(); ++i){
                if (thrust_.at(i)>0)
                {
                    Results<<t_.at(i)<<"   "<<thrust_.at(i)<<"\n";
                }
            }
            Results.close();
        }
    }

    void print_massFlowRate(std::string datfile){
        std::ofstream Results;
        Results.open(datfile);
        if (Results.is_open())
        {
            for(int i = 0; i < t_.size(); ++i){
                Results<<t_.at(i)<<"   "<<m_dot_.at(i)<<"\n";
            }
            Results.close();
        }
    }

    void print_exitTemperature(std::string datfile){
        std::ofstream Results;
        Results.open(datfile);
        if (Results.is_open())
        {
            for(int i = 0; i < t_.size(); ++i){
                Results<<t_.at(i)<<"   "<<exit_temperature_.at(i)<<"\n";
            }
            Results.close();
        }
    }

    
    Nozzle() = default;
    Nozzle(const Nozzle &other) = default;
    ~Nozzle() = default;

};

