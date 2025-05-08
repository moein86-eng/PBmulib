#ifndef __PBE_0D_POSTPROCESS_DECLARATION__ 
#define __PBE_0D_POSTPROCESS_DECLARATION__

#include "Distribution_declaration.h"

class PBE_0D_PostProcess
{
public:
	PBE_0D_PostProcess(PBE_0D_Solution results);  
	~PBE_0D_PostProcess();
	void PostProcess();                           //Performs post-processing

private:
	void clean_up(Distribution arr[]);            //clean up the dynamic meamory as required

public:
	bool Droplet_Grid_Volume_bool;		    //boolean for generating inner variable based on doplet volume (non dimensionless, [m3])
	bool Droplet_Grid_Diameter_bool;		//boolean for generating inner variable based on doplet diameter (non dimensionless, [m])
	bool Droplet_Grid_Radius_DL_bool;		//boolean for generating inner variable based on doplet radius (dimensionless)
	bool Droplet_Grid_Volume_DL_bool;		//boolean for generating inner variable based on doplet volume (dimensionless)

	bool Density_VBased_DL_bool;			//boolean for generating dimensionless density based on doplet volume (dimensionless)
	bool Volume_Density_RBased_bool;		//boolean for generating volume density based on doplet radius (non dimensionless, [m-1])
	bool Volume_Density_VBased_bool;		//boolean for generating volume density based on doplet volume (non dimensionless, [m-3])
	bool Number_Density_RBased_bool;		//boolean for generating number density based on doplet radius (non dimensionless, [m-3m-1])
	bool Number_Density_VBased_bool;		//boolean for generating number density based on doplet volume (non dimensionless, [m-3m-3])

	bool moment__RBased_bool;		    //boolean for generating moment of order "degree" based on doplet radius (dimensionless)
	bool moment__VBased_bool;	    	//boolean for generating moment of order "degree" based on doplet volume (dimensionless)
	short int moment__RBased_degree;    //degree for generating moment based on doplet volume
	short int moment__VBased_degree;

	bool Total_Number_bool;					//boolean for generating total droplet number (non dimensionless, [m-3])
	bool Total_Number_DL_bool;				//boolean for generating total droplet number (dimensionless)

	bool Volume_Fraction_bool;               //boolean for generating volume fraction (dimensionless)

	bool Average_Radius_RBased_DL_bool;		//boolean for generating average droplet radius based on doplet radius (dimensionless)
	bool Average_Volume_RBased_DL_bool;		//boolean for generating average droplet volume based on doplet radius (dimensionless)
	bool Average_Radius_VBased_DL_bool;		//boolean for generating average droplet radius based on doplet volume (dimensionless)
	bool Average_Volume_VBased_DL_bool;		//boolean for generating average droplet volume based on doplet volume (dimensionless)
	bool Average_Radius_RBased_bool;			//boolean for generating average droplet radius based on doplet radius (non dimensionless, [m])
	bool Average_Volume_RBased_bool;			//boolean for generating average droplet volume based on doplet radius (non dimensionless, [m3])
	bool Average_Radius_VBased_bool;			//boolean for generating average droplet radius based on doplet volume (non dimensionless, [m])
	bool Average_Volume_VBased_bool;			//boolean for generating average droplet volume based on doplet volume (non dimensionless, [m3])

	bool Standard_Deviation_RBased_DL_bool;	//retuns standard deviation based on doplet radius (dimensionless)
	bool Standard_Deviation_VBased_DL_bool;	//retuns standard deviation based on doplet volume (dimensionless)
	bool Standard_Deviation_RBased_bool;	    //retuns standard deviation based on doplet radius (non dimensionless, [m])
	bool Standard_Deviation_VBased_bool;		//retuns standard deviation based on doplet volume (non dimensionless, [m3])

	bool Sauter_Mean_Radius_DL_bool;			//retuns sauter mean radius (dimensionless)
	bool Sauter_Mean_Radius_bool;			//retuns sauter mean radius (non dimensionless, [m])

	PBE_0D_Solution results;               //raw simulation results
	Distribution* results_pp;              //post processed result
};
#endif