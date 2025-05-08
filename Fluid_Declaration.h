#ifndef __FLUID_DECLARATION__
#define __FLUID_DECLARATION__

/*
unit of measurements:
density				kg/m3
viscosity			Pa.s or 0.001xcp
Temperature			degC
Salinity			ppm or mg/L
*/

//#include "pch.h"

#include <string>

class Fluid
{
public:
	void Set_Density(double Density);					//used for direct user input
	void Set_Viscosity(double Viscosity);				//used for direct user input

	void Set_Density(double* Density_Array, double Temperature);			//used for direct data entery and interpolation
	void Set_Viscosity(double* Viscosity_Array, double Temperature);		//used for direct data entery and interpolation

	void Set_Density(std::string Formula_str, std::string Temprature_str, const std::string* Parameter_str, double* Parameter_Val, unsigned int Parameter_number, double Temperature);			//used for parsing user specified formula with parameters
	void Set_Viscosity(std::string Formula_str, std::string Temprature_str, const std::string* Parameter_str, double* Parameter_Val, unsigned int Parameter_number, double Temperature);		//used for parsing user specified formula with parameters

	void Set_Density_DeadOil(std::string Density_Model_String, double Temperature);			//used for predefined each crude oil models (only crude oil)
	void Set_Viscosity_DeadOil(std::string Viscosity_Model_String, double Temperature);		//used for predefined each crude oil models (only crude oil)

	void Set_Density_DeadOil(std::string Density_Model_String, double STD_Density, double Temperature);			//used for predefined general crude oil models (only crude oil)
	void Set_Viscosity_DeadOil(std::string Viscosity_Model_String, double STD_Density, double Temperature);		//used for predefined general crude oil models (only crude oil)

	void Set_Density_LiveOil(std::string Density_SolubilityModel_String, double Gas_SG, double Pressure);			//used for predefined general crude oil models (only for correction of dissolved gas on crude oil)
	void Set_Viscosity_LiveOil(std::string Viscosity_SolubilityModel_String, double Gas_SG, double Pressure);		//used for predefined general crude oil models (only for correction of dissolved gas on crude oil)

	void Set_Density_Water(std::string Density_Model_String, double Temperature, double Salinity);			//used for predefined models (only saline water)
	void Set_Viscosity_Water(std::string Viscosity_Model_String, double Temperature, double Salinity);		//Temperature shall be in degC, used for predefined models (only saline water)

	//pvt files reading
	//compositional EOS

public:
	bool Aqueous_Flag;      //0 for oil and 1 for water
	double Temperature;       
	double Density;
	double Viscosity;

};

#endif