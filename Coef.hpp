#pragma once
#include <iostream>
#include <cmath>


enum prop_coef {BurnRateCoef, BurnRateExponent, Density, CombustionTemperature, SpecificHeatRatio, ExchaustMolarMass, GasConstant};
enum other_coef {NumberOfGrains, NumberOfGaps, Gap, InnerDiametre, OutDiametre, Length};

double pi = 2*acos(0);

// throat diametre
double d_t = 15.56e-3;
// throat cross section area
double A_t = 0.25*pi*pow(d_t,2);
// exit diametre
double d_e = 44.26e-3;
// throat cross section area
double A_e = 0.25*pi*pow(d_e,2);

double AmbientPressure = 1e5;
double KNDX_PROPELLANT[GasConstant+1] = {
    11e-6, //BurnRateCoef
    0.4437, //BurnRateExponent
    1783.22563, //Density
    1625, //CombustionTemperature
    1.1308, //SpecificHeatRatio
    42.39e-3, //ExchaustMolarMass
    8.314463/KNDX_PROPELLANT[ExchaustMolarMass] //GasConstant
};
double LengthOfCasing = 375e-3;

double GRAIN[Length+1] = {
    4, //NumberOfGrains
    5, //NumberOfGaps
    5e-3, //Gap
    20e-3, //InnerDiametre
    56e-3, //OutDiametre
    (LengthOfCasing - GRAIN[NumberOfGaps]*GRAIN[Gap])/4 //Length
};



double M_prop = GRAIN[NumberOfGrains]*3.14*GRAIN[Length]*(0.25*(pow(GRAIN[OutDiametre],2) - pow(GRAIN[InnerDiametre],2)))*KNDX_PROPELLANT[Density];

double Sy = 200e6; // Yield strength of casing material

double BurningCaseCoef = 1; //1 -> great combustion
                               //0 -> bad combustion
double tc = 0.2;
