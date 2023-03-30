#pragma once
#include "Coef.hpp"
#include <fstream>
#include <sstream>
#include <string>


void readData(std::string nameFile) {
    std::ifstream inFile(nameFile);

    if (inFile.is_open()) {
        std::string line;

        while (getline(inFile, line)) {
            std::stringstream ss(line);
            std::string valueString, nameString;
            double value;

            ss >> valueString >> nameString;

            try {
                value = std::stod(valueString);
            }
            catch (const std::invalid_argument& e) {
                std::cerr << "Invalid argument: " << e.what() << '\n';
                continue;
            }

            if (nameString == "<BurnRateCoef>")
                KNDX_PROPELLANT[BurnRateCoef] = value;
            else if (nameString == "<BurnRateExponent>")
                KNDX_PROPELLANT[BurnRateExponent] = value;
            else if (nameString == "<Density>")
                KNDX_PROPELLANT[Density] = value;
            else if (nameString == "<CombustionTemperature>")
                KNDX_PROPELLANT[CombustionTemperature] = value;
            else if (nameString == "<SpecificHeatRatio>")
                KNDX_PROPELLANT[SpecificHeatRatio] = value;
            else if (nameString == "<ExchaustMolarMass>")
                KNDX_PROPELLANT[ExchaustMolarMass] = value;
            else if (nameString == "<GasConstant>")
                KNDX_PROPELLANT[GasConstant] = value;
            else if (nameString == "<NumberOfGrains>")
                GRAIN[NumberOfGrains] = value;
            else if (nameString == "<NumberOfGaps>")
                GRAIN[NumberOfGaps] = value;
            else if (nameString == "<Gap>")
                GRAIN[Gap] = value;
            else if (nameString == "<InnerDiametre>")
                GRAIN[InnerDiametre] = value;
            else if (nameString == "<OutDiametre>")
                GRAIN[OutDiametre] = value;
            else if (nameString == "<Length>")
                GRAIN[Length] = value;
            else {
                std::cerr << "Unknown coefficient name: " << nameString << '\n';
                continue;
            }
        }

        inFile.close();
    }
    else {
        std::cerr << "Unable to open file: " << nameFile << '\n';
    }
}
