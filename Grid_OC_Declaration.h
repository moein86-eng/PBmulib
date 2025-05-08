#ifndef __GRID_OC_DECLARATION__ 
#define __GRID_OC_DECLARATION__

#include <cmath>
//#include "pch.h"

class Grid_OC
{
public:
	Grid_OC(unsigned int colloc_points, unsigned int left = 1, unsigned int right = 1, double R_max=0, double t_c=1);
	~Grid_OC();
	double** _1st_derivatic_weight();     //returns an 2D array for the 1st derivative weights
	double** _2nd_derivatic_weight();     //returns an 2D array for the 2nd derivative weights

private:
	double* jcobi_roots(double alpha = 0, double beta = 0);       //return an array for the jacobian polynomial roots  
	double* Quadrature_weight();          //returns an array for the quadrature weights
	double* dif();                        //internal function for the method
	double* dfopr(const unsigned int &i, const unsigned int &id, double dif[]);  //internal function for the method
	void clean_up(double dif[]);          //clean up the dynamic meamory as required (one dimension array)
	void clean_up(double** arr);		  //clean up the dynamic meamory as required (two dimension array)

private:
	unsigned int colloc_points;           //number of middle collocation points
	unsigned int left;					  //lower boundary point (1 if included, 0 if not)
	unsigned int right;                   //upper boundary point (1 if included, 0 if not)
	bool _1st_deriv_bool;                 //boolian only used for deconstructor 
	bool _2nd_deriv_bool;                 //boolian only used for deconstructor
	double Tol;                           //Newton method tolerance for root determination

public:
	unsigned int TotalPointsNumber;            //number of total orthogonal collocation points
	double* Grid_Points_;                 // array of jacobian roots
	double* Quadrature_weight_;           //array quadrature weight
	double** _1st_derivatic_weight_;      //will be calculated only if respective function is called
	double** _2nd_derivatic_weight_;      //will be calculated only if respective function is called
	double R_max;                         //maximum droplet radius used for generating the inner dimension vector
	double t_c;                                 //characteristic time for dimensionalizing (default value is 1 sec)
};

#endif