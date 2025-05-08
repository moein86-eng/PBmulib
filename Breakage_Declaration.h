#ifndef __BREAKAGE_DECLARATION__
#define __BREAKAGE_DECLARATION__

//#include "pch.h"

#include "exprtk.hpp"
#include <string>
#include <cmath>
#include <regex>

#include "MP_Fluid_Declaration.h"


template<typename Grid_type>
class Breakage
{
public:
	Breakage(const Grid_type* Grid, const MP_Fluid* Fluid_Val, std::string Formula_str, std::string Variable_str, std::string VarType_str, const std::string* Parameter_str, double* Parameter_Val, unsigned int Parameter_number, const std::string* Fluid_str, const std::string* Fluid_ID, unsigned int Fluid_number);   //constructor used for user-specified relation using parser
	Breakage(const Grid_type* Grid, const MP_Fluid* Fluid_Val, std::string ModelName_str, double* Parameter_Val);   //constructor used for included models from literature
	Breakage(const Grid_type* Grid, const MP_Fluid* Fluid_Val, double* B_f, double** B_f_Shifted);   //constructor used for importing array of the calculated values
	~Breakage();
	std::string Get_Formula() const;      //returns formula string
	std::string Get_Variable() const;     //returns variable string
	unsigned int Get_Parameter_number() const;     //return number of parameters
	std::string Get_Parameter_str();     //returns all parameters' strings
	double* Get_Parameter_Val() const;    //returns a pointer to parameter value array
	unsigned int Get_Fluid_number() const;     //return number of properties
	std::string Get_Fluid_str();       //returns all the Fluid strings
	double* Get_Fluid_Val() const;     //returns all Fluid value strings
	std::string Get_VarType() const;      //returns variable type string ("volume_DL" or "radius_DL")

private:
	void clean_up(double* arr);		      //clean up the dynamic meamory as required (two dimension array)
	void clean_up(double** arr);			//clean up the dynamic meamory as required (two dimension array)
	void Parser_fun();                 //parse user-specified string relations 
	void Parser_fun(std::string For_str, std::string r); //parse user-specified string relations (used for converted string from "volume_DL" to "radius_DL")
	void Model1();					  //model1 from literature
	void Model2();					  //model2 from literature
	std::string vol_to_rad();             //returns converted string from "volume_DL" to "radius_DL"

private:
	double* Parameter_Val;                //array quadrature weight
	unsigned int Parameter_number;        //number of parameters in the relation
	const std::string* Parameter_str;           //array to quadrature weight
	std::string Parameter_str_sum;        //string for all parameters
	unsigned int Fluid_number;			  //number of physio-chemical fluid properties in the relation
	const std::string* Fluid_str;				  //array of physio-chemical fluid properties strings
	std::string Fluid_str_sum;			  //string for all Properties
	const std::string* Fluid_ID;		     	  //string for all fluid Properties as identifiers

	std::string Formula_str;              //string for breakage formula
	std::string Variable_str;             //string for grid variable
	std::string VarType_str;              //type of variable in distribution formula ("volume_DL" or "radius_DL"). Only dimensionless volume or radius can be used.

public:
	const Grid_type* Grid;                //pointer to the grid
	const MP_Fluid* Fluid_Val;	      //pointer to the Fluid type (physio-chemical properties values)
	double* B_f;						//array for breakage frequency (s-1). 
	double** B_f_Shifted;				//array for shifted breakage frequency (s-1). 
};

#endif