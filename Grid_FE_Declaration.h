#ifndef __GRID_FE_DECLARATION__ 
#define __GRID_FE_DECLARATION__

//#include "pch.h"

#include <cmath>
#include <iostream>

class Grid_FE
{
public:
	Grid_FE(const double* ElementBoundaryNodes, const unsigned int* ElementPointsNumber, unsigned int ElementNumber, bool _1st_derivatic_bool = false, bool _2nd_derivatic_bool = false, double R_max=0 , double t_c = 1);
	Grid_FE(double Int_01, double Mean, double Int_99, const unsigned int* ElementPointsNumber, unsigned int ElementNumber, bool _1st_derivatic_bool = false, bool _2nd_derivatic_bool = false, double R_max = 0, double t_c = 1); //automatic grid based on initial condition

	~Grid_FE();

private:
	void clean_up(double arr_1D[]);				//clean up the dynamic meamory as required (one dimension array)
	void clean_up(double** arr_2D);				//clean up the dynamic meamory as required (two dimension array)
	void weights();                             //generates gid and weights

private:
	const double* ElementBoundaryNodes;         //an array of left and right nodes for all elements
	const unsigned int* ElementPointsNumber;    //an array of number of total points in each element 
	unsigned int ElementNumber;                 //number of elements
	bool _1st_derivatic_bool;					//boolian only used for deconstructor 
	bool _2nd_derivatic_bool;					//boolian only used for deconstructor 
	bool Boundary_Delete;						//boolian only used for deconstructor 

public:
	unsigned int TotalPointsNumber;				//number of total points in the global coordinate
	double* Grid_Points_;						//generated grid points for all elements according to global coordinate 
	double* Quadrature_weight_;					//generated quadrature weight array for all elements according to global coordinate
	double** _1st_derivatic_weight_;			//generated 1st derivative weight matrix for all elements according to global coordinate
	double** _2nd_derivatic_weight_;			//generated 1st derivative weight matrix for all elements according to global coordinate
	double R_max;								//maximum droplet radius used for generating the inner dimension vector
	double t_c;                                 //characteristic time for dimensionalizing (default value is 1 sec)
};

#endif