#ifndef __LOGGER_DECLARATION__ 
#define __LOGGER_DECLARATION__

//#include "pch.h"

#include <iostream>
#include <string>

class Logger
{          
public:
	Logger(const Logger&) = delete;

	static Logger& Instance();

	void Add_Error(std::string _str);
	void Add_Warning(std::string _str);
	void Add_Info(std::string _str);

	void Print_Error();
	void Print_Warning();
	void Print_Info();

	const std::string& Get_Error_str() const;
	const std::string& Get_Warning_str() const;
	const std::string& Get_Info_str() const;
	   	
	const bool& Get_Error_flag() const;
	const bool& Get_Warning_flag() const;
	const bool& Get_Info_flag() const;

private:
	Logger(); 
	std::string Error_str;  
	std::string Warning_str;
	std::string Info_str;
	bool Error_flag;
	bool Warning_flag;
	bool Info_flag;
};

#endif