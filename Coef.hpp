#pragma once
#include <iostream>


enum prop_coef {BurnRateCoef, BurnRateExponent, Density, CombustionTemperature, SpecificHeatRatio, ExchaustMolarMass, GasConstant};
enum other_coef {NumberOfGrains, NumberOfGaps, Gap, InnerDiametre, OutDiametre, Length};

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
double GRAIN[Length+1] = {
    4, //NumberOfGrains
    0, //NumberOfGaps
    5e-3, //Gap
    20e-3, //InnerDiametre
    56e-3, //OutDiametre
    100e-3 //Length
};

double Sy = 200e6; // Yield strength of casing material

double BurningCaseCoef = 1; //1 -> great combustion
                               //0 -> bad combustion
double tc = 0.3;
