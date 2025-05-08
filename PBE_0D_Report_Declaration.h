#ifndef __PBE_0D_REPORT_DECLARATION__
#define __PBE_0D_REPORT_DECLARATION__

//#include "pch.h"
#include<string>
#include "PBE_0D_Solution_Declaration.h"
#include "PBE_0D_PostProcess_Declaration.h"


class PBE_0D_Report
{
public:
	PBE_0D_Report(PBE_0D_Solution results, PBE_0D_PostProcess results_pp);   //constructor for 0D-PBE
	
	std::string Get_path() const;             //returns the directory path for saving files
	void Set_path(std::string path);          //sets the directory path for saving files

private:
	void generate_text();		          //generates the main simulation report

	void generate_binary(double* data, std::string name, unsigned int size);		      //generates a binary file from a 1D vector
	void generate_binary(double** data, std::string name, unsigned int size);		      //generates a binary file from a 2D matrix

private:
	std::string path;         //string for path of the directory to save the files
	PBE_0D_Solution results;
};

#endif