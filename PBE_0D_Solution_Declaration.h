#ifndef __PBE_0D_SOLUTION_DECLARATION__
#define __PBE_0D_SOLUTION_DECLARATION__
#pragma warning(disable: 4996)

#include "Distribution_declaration.h"
#include "Daughter_Declaration.h"
#include "Breakage_Declaration.h"
#include "Coalescence_Declaration.h"

//#include "pch.h"

#include <cmath>
#include <regex>
#include <iostream>
#include <chrono>
#include <vector>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/container/vector.hpp>
#include <boost/numeric/ublas/vector.hpp>

//#include <boost/numeric/odeint/external/vexcl/vexcl.hpp>

typedef std::vector<double*> result_dist_type;
typedef std::vector<double> result_time_type;

//typedef std::vector<double> state_type;                       //3.5 sec
//typedef boost::array< double, 238 > state_type;               //1.9 sec
//typedef std::array<double, 238> state_type;                   //2.2 sec
//typedef  boost::container::vector<double> state_type;         //did not run
typedef boost::numeric::ublas::vector< double > state_type;     //more than 2 times slower than std::vector

template<typename Grid_Type>
class PBE_0D_Solution
{
public:
	PBE_0D_Solution(Distribution<Grid_Type>* Int_Dist, Daughter<Grid_Type>* Daughter_, Breakage<Grid_Type>* Breakage_, Coalescence<Grid_Type>* Coalescence_);   //constructor for 0D-PBE module without inlet and outlet (both coalescence and breakage)
	PBE_0D_Solution(Distribution<Grid_Type>* Int_Dist, Distribution<Grid_Type>* Inlet_Dist, Daughter<Grid_Type>* Daughter_, Breakage<Grid_Type>* Breakage_, Coalescence<Grid_Type>* Coalescence_, double Residence_Time);   //constructor for 0D-PBE module with inlet and outlet (both coalescence and breakage)
	PBE_0D_Solution(Distribution<Grid_Type>* Int_Dist, Daughter<Grid_Type>* Daughter_, Breakage<Grid_Type>* Breakage_, Coalescence<Grid_Type>* Coalescence_, double Velocity);   //constructor for 0D-PBE PipeFlow module (both coalescence and breakage)
	
	PBE_0D_Solution(Distribution<Grid_Type>* Int_Dist, Coalescence<Grid_Type>* Coalescence_);   //constructor for 0D-PBE module without inlet and outlet (coalescence only)
	PBE_0D_Solution(Distribution<Grid_Type>* Int_Dist, Distribution<Grid_Type>* Inlet_Dist, Coalescence<Grid_Type>* Coalescence_, double Residence_Time);   //constructor for 0D-PBE module with inlet and outlet (coalescence only)
	PBE_0D_Solution(Distribution<Grid_Type>* Int_Dist, Coalescence<Grid_Type>* Coalescence_, double Velocity);   //constructor for 0D-PBE PipeFlow module (coalescence only)
	
	PBE_0D_Solution(Distribution<Grid_Type>* Int_Dist, Daughter<Grid_Type>* Daughter_, Breakage<Grid_Type>* Breakage_);   //constructor for 0D-PBE module without inlet and outlet (Breakage only)
	PBE_0D_Solution(Distribution<Grid_Type>* Int_Dist, Distribution<Grid_Type>* Inlet_Dist, Daughter<Grid_Type>* Daughter_, Breakage<Grid_Type>* Breakage_, double Residence_Time);   //constructor for 0D-PBE module with inlet and outlet (Breakage only)
	PBE_0D_Solution(Distribution<Grid_Type>* Int_Dist, Daughter<Grid_Type>* Daughter_, Breakage<Grid_Type>* Breakage_, double Velocity);   //constructor for 0D-PBE PipeFlow module (Breakage only)

	PBE_0D_Solution(Distribution<Grid_Type>* Int_Dist, Distribution<Grid_Type>* Inlet_Dist, Grid_Type* Grid_NBC, double Residence_Time);   //constructor for 0D-PBE module with inlet and outlet (no coalescence or breakage)

	~PBE_0D_Solution();
	
	void Set_Solver(std::string Integrator_lib="ODEINT", std::string Integrator_Met="runge_kutta_cash_karp54", double Integrator_RelTol=1e-6, double Integrator_AbsTol=1e-6, unsigned int Save_Option =1, double dt_Initial=0.1);      //sets time integrator options-inputs
	void Simulate(double Final_Time_Length);      //solve the set ode equations, input should be time if 0D-PBE module with or without input/output to be used, for pipeflow module length to be used.

private:
	void clean_up(double** arr);	  	  //clean up the dynamic meamory as required (two dimension array)
	void clean_up(double arr[]);          //clean up the dynamic meamory as required
	void Coalescence_Death_Birth_Coef();  //generates required constant death / birth matrices for coalescence terms (3 matrices)
	void Breakage_Death_Birth_Coef();     //generates required constant death / birth matrices for breakage terms (2 matrices)
	double* PCHIP_Mono_Interpolation(const double* User_grid, const state_type& User_dist, const double* R_grid, const unsigned int& user_N, const unsigned int& Grid_N, double* yhat, double Inter_Para = 1.5); //piecewise cubic Hermite interpolating polynomial specially adapted for monotonicity (Inter_Para shall be between 1.5 to 3, 1.5 provides similar results to MATLAB)
	unsigned int Binary_Search(const double* x, const double xhat, const unsigned int n);    //Binary search to find index i such that x(i) <= xhat <= x(i + 1)
	void ODE_Set_B(const state_type &sai, state_type& dsai_dt, const double& t);		 //returns time differntial vector
	void ODE_Set_C(const state_type& sai, state_type& dsai_dt, const double& t);        //returns time differntial vector
	void ODE_Set_BC(const state_type& sai, state_type& dsai_dt, const double& t);       //returns time differntial vector
	void ODE_Set_B_IO(const state_type& sai, state_type& dsai_dt, const double& t);     //returns time differntial vector
	void ODE_Set_C_IO(const state_type& sai, state_type& dsai_dt, const double& t);     //returns time differntial vector
	void ODE_Set_BC_IO(const state_type& sai, state_type& dsai_dt, const double& t);    //returns time differntial vector
	void ODE_Set_NBC_IO(const state_type& sai, state_type& dsai_dt, const double& t);   //returns time differntial vector
	void Save_Results(const state_type& sai, const double& t);                         //auxillary function to save the results 

private:
	double** Br_B_Coef;					  //constant birth matrix coefficient for breakage
	double* Br_D_Coef;				      //constant death matrix coefficient for breakage
	double** Co_B_Coef;					  //constant birth matrix coefficient for coalescence
	double** Co_D_Coef;					  //constant death matrix coefficient for coalescence
	double** X_Mat_1;					  //constant matrix for local to global remapping used if coalescence is involved
	double** X_Mat_2;					  //constant matrix for local to global remapping used if coalescence is involved
	double** X_Mat_3;					  //constant matrix for local to global remapping used if Breakage is involved

	double* sai1;						  //temporary array
	double* sai2;						  //temporary array
	double* sai3;						  //temporary array
	double* dy_dx;						  //temporary array
	double* m;							  //temporary array
	double* c_coef;						  //temporary array
	double* d_coef;						  //temporary array
	unsigned int Index;                   //temporary index
	unsigned int Grid_Size;               //temporary index
	double Time_DL;						  //dimensionless time (characteristic time / residence time)
	void (PBE_0D_Solution::* ODE_Set_fun)(const state_type& sai, state_type& dsai_dt, const double& t);     //pointer to the ODE set function 

public:
	double Velocity;					  //velocity in the pipe, only used for pipe flow module
	double Residence_Time;				  //Residence time, only used for CSTR module
	
	std::string Integrator_Lib;			  //time integrator library "odeint" / "sundail"
	std::string Integrator_Method;		  //time integrator method "RK_fehlberg78" / "RK_cash_karp54" / "RK_dopri5"
	double RelTol;						  //time integrator relative tolerance
	double AbsTol;						  //time integrator absolute tolerance
	unsigned int Save_Option;             //Save_Option specifies the number of steps when the data should be saved
	double dt_Initial;                    //time integrator starting time step

	bool Flag_Solution;                   //0: for not solved, 1: solved 
	unsigned int Flag_Module;			  //1: without inlet and outlet module, 2: with inlet and outlet module, 3: PipeFlow module
	unsigned int Flag_BC;				  //0: No breakage or coalescence, 1: Both Breakage and Coalescence, 2: only Coalescence, 3: only Breakage
	double Final_Time;					  //final silulation time
	double Final_Length;				  //final silulation length
				  
	result_dist_type Dist_Sol;            //generated solution matrix by time integrator based on dimensionless distribution (radius-based)
	result_time_type Time_Sol;            //generated time vector from the time integrator
	unsigned int Step_Number;			  //number of solution points with respect to time
	
	double Execution_Time;                //solver execution time in micro-sec
	double PrePros_Time;                  //solver pre-processing time in micro-sec

	Distribution<Grid_Type>* Int_Dist;	  //initial distribution vector
	Distribution<Grid_Type>* Inlet_Dist;  //inlet distribution vector
	Daughter<Grid_Type>* Daughter_;		  //daughter distribution vector
	Breakage<Grid_Type>* Breakage_;		  //breakage frequency vector
	Coalescence<Grid_Type>* Coalescence_; //coalescence rate vector
	Grid_Type* Grid_NBC;                  //Grid only used when no breakage or coalescence is used
};

#endif