#ifndef __OW_INTERFACE_DECLARATION__
#define __OW_INTERFACE_DECLARATION__

//#include "pch.h"

#include<string>
class OW_Interface
{
public:
	OW_Interface(double* Grid_Points_, unsigned int Grid_Number, std::string Formula_str, std::string Variable_str, std::string VarType_str, const std::string* Parameter_str, double* Parameter_val, unsigned int Parameter_number);   //constructor used for user-specified relation using parser
	OW_Interface(double* Grid_Points_, unsigned int Grid_Number, std::string ModelName_str, double* Parameter_val);   //constructor used for included models from literature
	~OW_Interface();
	std::string Get_Formula() const;      //returns formula string
	std::string Get_Variable() const;     //returns variable string
	std::string* Get_Parameter_str() const;     //returns a pointer to parameter string array
	double* Get_Parameter_val() const;    //returns a pointer to parameter value array
	std::string Get_VarType() const;      //returns variable type string ("volume_DL" or "radius_DL")

private:
	void clean_up(double** arr);		  //clean up the dynamic meamory as required (two dimension array)
	double** Parser_fun();                //parse user-specified string relations 
	double** Parser_fun(std::string For_str, std::string r); //parse user-specified string relations (used for converted string from "volume_DL" to "radius_DL")
	double** Model1();					  //model1 from literature
	double** Model2();					  //model2 from literature
	std::string vol_to_rad();             //returns converted string from "volume_DL" to "radius_DL"

private:
	unsigned int Grid_Number;             //number of grid points
	const double* Grid_Points_;           //array of grid points
	double* Parameter_val;                //array quadrature weight
	unsigned int Parameter_number;        //number of parameters in the relation
	std::string* Parameter_str;           //array to quadrature weight
	std::string Formula_str;              //string for coalescence formula
	std::string Variable_str;             //string for grid variable
	std::string VarType_str;              //type of variable in distribution formula ("volume_DL" or "radius_DL"). Only dimensionless volume or radius can be used.

public:
	double** IC_r;						  //array for interfacial Coalescence rate. 
};

#endif