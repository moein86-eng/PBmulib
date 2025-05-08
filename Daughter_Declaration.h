#ifndef __DAUGHTER_DECLARATION__
#define __DAUGHTER_DECLARATION__

//#include "pch.h"

#include "exprtk.hpp"
#include <string>
#include <cmath>
#include <regex>

template<typename Grid_type>
class Daughter
{
public:
	Daughter(const Grid_type* Grid, std::string Formula_str, std::string Variable_str_j, std::string Variable_str_i, std::string VarType_str = "radius_DL");   //constructor used for user-specified relation using parser
	Daughter(Grid_type* Grid, std::string ModelName_str = "Model1");   //constructor used for included models from literature
	Daughter(const Grid_type* Grid, double** beta, double** beta_Shifted);   //constructor used for importing array of the calculated values
	~Daughter();
	std::string Get_Formula() const;      //returns formula string
	std::string Get_Variable() const;     //returns variable string
	std::string Get_VarType() const;      //returns variable type string ("volume_DL" or "radius_DL")
	double Get_Int();                     //returns daughter distribution integration over domain

private:
	void clean_up(double** arr);		  //clean up the dynamic meamory as required (two dimension array)
	void InputChecking();                 //checks the relations and returns the error or warning logs (for detail refer to the main function body)
	void Parser_fun();                //parse user-specified string relations 
	void Parser_fun(std::string For_str, std::string rj, std::string ri); //parse user-specified string relations (used for converted string from "volume_DL" to "radius_DL")
	void Model1();					  //model1 from literature
	void Model2();					  //model2 from literature
	std::string vol_to_rad();             //returns converted string from "volume_DL" to "radius_DL"

private:
	std::string Formula_str;              //string for daughter distribution formula
	std::string Variable_str_j;           //string for mother grid variable
	std::string Variable_str_i;           //string for daughter grid variable
	std::string VarType_str;              //type of variable in distribution formula ("volume_DL" or "radius_DL"). Only dimensionless volume or radius can be used.

public:
	const Grid_type* Grid;                      //pointer to the grid
	double** beta;						  //array for daughter distribution. (beta_b(j,i)   j>=i j is breaked down to i) / beta is dimensionless (beta_DL=beta_NDL/R_max) 
	double** beta_Shifted;						//array for shifted daughter distribution. (beta_b(j,i)   j>=i j is breaked down to i) / beta is dimensionless (beta_DL=beta_NDL/R_max) 
};

#endif