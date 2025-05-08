#ifndef __DISTRIBUTION_DECLARATION__ 
#define __DISTRIBUTION_DECLARATION__

//#include "pch.h"

#include <cmath>

/*constructor main input:
dependent variable: dimensionless radius
independent variable: dimensionless distribution (radius-based)*/


//notes
//1: add formula parser functions
//2: add more prespecified distribution functions
//3: MP_Fluid should be used in it

template<typename Grid_type>
class Distribution
{
public:
	Distribution(const Grid_type* Grid, double* dist_sai_RB);  //(for post processing)when the LogNormal distribution is used, and array for grid points are defined.
	Distribution(const Grid_type* Grid, double mu, double fi, double VolFrac); //(for initial condition)when the LogNormal distribution is used, and array for grid points are defined.
	Distribution(const Grid_type* Grid, const double* User_dist, const double* User_grid, unsigned int user_N); //(for initial condition)when the both grid points and distribution points are defined.
	Distribution(double* User_dist, double* User_grid, unsigned int user_N); //(for auto grids)
	Distribution(double mu, double fi, double VolFrac, unsigned int user_N); //(for auto grids)

	void Set_Grid(const Grid_type* Grid);//(for auto grids)

	~Distribution();

	double* Droplet_Grid_Radius();			//returns inner variable based on doplet radius (non dimensionless, [m])
	double* Droplet_Grid_Volume();		    //returns inner variable based on doplet volume (non dimensionless, [m3])
	double* Droplet_Grid_Diameter();		//returns inner variable based on doplet diameter (non dimensionless, [m])
	double* Droplet_Grid_Radius_DL();		//returns inner variable based on doplet radius (dimensionless)
	double* Droplet_Grid_Volume_DL();		//returns inner variable based on doplet volume (dimensionless)

	double* Density_RBased_DL();			//returns dimensionless density based on doplet radius (dimensionless)
	double* Density_VBased_DL();			//returns dimensionless density based on doplet volume (dimensionless)
	double* Volume_Density_RBased();		//returns volume density based on doplet radius (non dimensionless, [m-1])
	double* Volume_Density_VBased();		//returns volume density based on doplet volume (non dimensionless, [m-3])
	double* Number_Density_RBased();		//returns number density based on doplet radius (non dimensionless, [m-3m-1])
	double* Number_Density_VBased();		//returns number density based on doplet volume (non dimensionless, [m-3m-3])
	
	double moment__RBased(short int degree);		//returns the moment of order "degree" based on doplet radius (dimensionless)
	double moment__VBased(short int degree);		//returns the moment of order "degree" based on doplet volume (dimensionless)

	double Total_Number();					//retuns total droplet number (non dimensionless, [m-3])
	double Total_Number_DL();				//retuns total droplet number (dimensionless)

	double Volume_Fraction();               //retuns volume fraction (dimensionless)

	double Average_Radius_RBased_DL();		//retuns average droplet radius based on doplet radius (dimensionless)
	double Average_Volume_RBased_DL();		//retuns average droplet volume based on doplet radius (dimensionless)
	double Average_Radius_VBased_DL();		//retuns average droplet radius based on doplet volume (dimensionless)
	double Average_Volume_VBased_DL();		//retuns average droplet volume based on doplet volume (dimensionless)
	double Average_Radius_RBased();			//retuns average droplet radius based on doplet radius (non dimensionless, [m])
	double Average_Volume_RBased();			//retuns average droplet volume based on doplet radius (non dimensionless, [m3])
	double Average_Radius_VBased();			//retuns average droplet radius based on doplet volume (non dimensionless, [m])
	double Average_Volume_VBased();			//retuns average droplet volume based on doplet volume (non dimensionless, [m3])

	double Standard_Deviation_RBased_DL() ;	//retuns standard deviation based on doplet radius (dimensionless)
	double Standard_Deviation_VBased_DL();	//retuns standard deviation based on doplet volume (dimensionless)
	double Standard_Deviation_RBased();	    //retuns standard deviation based on doplet radius (non dimensionless, [m])
	double Standard_Deviation_VBased();		//retuns standard deviation based on doplet volume (non dimensionless, [m3])

	double Sauter_Mean_Radius_DL();			//retuns sauter mean radius (dimensionless)
	double Sauter_Mean_Radius();			//retuns sauter mean radius (non dimensionless, [m])

	double head_x(double criteria = 0.1);		//returns distribution head (criteria is the percentage of peak to detect the head)
	double tail_x(double criteria = 0.1);	    //returns distribution tail (criteria is the percentage of peak to detect the tail)
	double peak_x();							//returns distribution peak 

	double Exp_Trapz(const double* x, const double* y, unsigned int N);    //calculates the integral of the user specified distribution using trapezoidal method
	double Exp_Int_99(const double* x, const double* y, unsigned int N);   //calculates the head point on user specified grid with more than 99% portion of the distribution integral
	double Exp_Int_01(const double* x, const double* y, unsigned int N);   //calculates the tail point on user specified grid with less than 01% portion of the distribution integral
	double Exp_Ave(const double* x, const double* y, unsigned int N);      //calculates the average on user specified grid based on the distribution

	double zerotail(const double* User_grid, const double* User_dist, const unsigned int &user_size);    //auxilary function for Lagrange interpolation
	double zerohead(const double* User_grid, const double* User_dist, const unsigned int &user_size);    //auxilary function for Lagrange interpolation
	double* Lagrange_Interpolation(const double* User_grid, const double* User_dist, const double* R_grid, const unsigned int &user_N, const unsigned int &Grid_N); //Lagrange interpolation specially adapted for probability distribution

	unsigned int Binary_Search(const double* x, const double xhat, const unsigned int n);       //Binary search to find index i such that x(i) <= xhat <= x(i + 1)
	double* PCHIP_Mono_Interpolation(const double* User_grid, const double* User_dist, const double* R_grid, const unsigned int& user_N, const unsigned int& Grid_N, double Inter_Para = 1.5); //piecewise cubic Hermite interpolating polynomial specially adapted for monotonicity (Inter_Para shall be between 1.5 to 3, 1.5 provides similar results to MATLAB)

	void InputChecking(const double* User_grid, const double* User_dist, const unsigned int &user_N);      //checks the input arrays and returns the error or warning logs (for detail refer to the main function body)

	void clean_up(double arr[]);            //clean up the dynamic meamory as required
	
	void Normalize(double VF0);						//normalize according to specified volume fraction only the main distribution array (dist_sai_RB), other values should be updated if required after


public:
	//double* R_grid;						    //array for dimensionless droplet radius (radius-based)
	double* dist_sai_RB;              //array for dimensionless distribution (radius-based)
	double* Droplet_Grid_Radius_;     //will be calculated only if respective function is called 
	double* Droplet_Grid_Volume_;		//will be calculated only if respective function is called
	double* Droplet_Grid_Diameter_;	//will be calculated only if respective function is called
	double* Droplet_Grid_Volume_DL_;	//will be calculated only if respective function is called
	double* Density_VBased_DL_;		//will be calculated only if respective function is called
	double* Volume_Density_RBased_;	//will be calculated only if respective function is called
	double* Volume_Density_VBased_;	//will be calculated only if respective function is called
	double* Number_Density_RBased_;	//will be calculated only if respective function is called
	double* Number_Density_VBased_;	//will be calculated only if respective function is called
	double* User_grid;
	double* User_dist;
	unsigned int user_N;

	const Grid_type* Grid;                  //pointer to the grid
	bool dist_sai_RB_Bool;					//boolian only used for different deconstructor
	bool Delete_User;                       //boolian only used for different deconstructor
	double Int_99;                          //distribution head used for auto grid generation
	double Int_01;                          //distribution tail used for auto grid generation
	double Mean;                            //distribution mean used for auto grid generation
};

#endif