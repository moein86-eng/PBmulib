#ifndef __MP_FLUID_DECLARATION__
#define __MP_FLUID_DECLARATION__

#include "Fluid_Declaration.h"


class MP_Fluid
{
public:
	MP_Fluid(const Fluid* Fluid1, const Fluid* Fluid2, double Vol_Frac, std::string Phase_Inversion_Model_Str = "Model1");    //volume fraction is defined as fraction of 1st fluid volume to the total volume

	void Set_Tur_Ener_Diss_Rate_Av(double Tur_Ener_Diss_Rate_Av);
	void Set_Tur_Ener_Diss_Rate_Av(double Power, double Vol); //power in (W), ro_c is continuous phase density in (kg/m3),  Vol is volume in m3

	void Set_Surface_Tension(double Surface_Tension);						//sets value for oil-water surface tension (N/m)
	void Set_Surface_Tension(double* Surface_Tension_Array, double Temperature);		//used for direct data entery and interpolation
	void Set_Surface_Tension(std::string Formula_str, std::string Temprature_str, const std::string* Parameter_str, double* Parameter_Val, unsigned int Parameter_number);		//used for parsing user specified formula with parameters
	void Set_Surface_Tension(std::string Surface_Tension_Model_String, double Temperature);		//used for predefined each crude oil models

private:
	void Phase_Inversion_Model1();    //determines the dispersed and continuous phases
	void Set_Density_Mix();          //calculates the mixture density according to weighted average
	void Set_Viscosity_Mix();        //calculates the mixture viscosity according to weighted average (there are two relations one for oil in water and one for water in oil)

public:
	double* Turbulence;              //Turbulence filed (m2/s2)
	double* Velocity;                //velocity field (m/s)
	double* Shear_Rate;              //shear rate field (s-1)

	double Reynolds;				 //Reynold number (only for pipe)
	bool Turbulent_Flag;			 //0 for laminar and 1 for turbulent (only for pipe)

	double Tur_Ener_Diss_Rate_Av;    //average turbulent energy dissipation rate (m2/s3)

	const Fluid* Fluid_C;            //pointer to the continuous fluid
	const Fluid* Fluid_D;            //pointer to the dispersed fluid
	double Vol_Frac;
	unsigned int Continuous_Flag;    //1 if Fluid1 is Continuous 2 if the Fluid2 is continuous

	double Density_Mix;			     //bulk density for the mixture
	double Viscosity_Mix;			 //bulk density for the mixture

	double Surface_Tension;          //surface tension (N/m)
};

#endif