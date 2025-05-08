#include "Logger_Declaration.h"

Logger::Logger() 
	: Error_str("(Error - Num)                           (Error - Description)\n"),
	Warning_str("(Warning-Num)                           (Warning-Description)\n"),
	Info_str("(Info-Num)                              (Info-Description)\n"),
	Error_flag(false), Warning_flag(false), Info_flag(false) {}

Logger& Logger::Instance() {
	static Logger theInstance;
	return theInstance;
}

void Logger::Add_Error(std::string _str) { 
	Error_str += _str;
	Error_flag = true;
}
void Logger::Add_Warning(std::string _str) {
	Warning_str += _str;
	Warning_flag = true;
}
void Logger::Add_Info(std::string _str) {
	Info_str += _str;
	Info_flag = true;
}

void Logger::Print_Error() { std::cout << Error_str << "\n\n"; }
void Logger::Print_Warning() { std::cout << Warning_str << "\n\n"; }
void Logger::Print_Info() { std::cout << Info_str << "\n\n"; }

const std::string& Logger::Get_Error_str() const { return Error_str; }
const std::string& Logger::Get_Warning_str() const { return Warning_str; }
const std::string& Logger::Get_Info_str() const { return Info_str; }

const bool& Logger::Get_Error_flag() const { return Error_flag; }
const bool& Logger::Get_Warning_flag() const { return Warning_flag; }
const bool& Logger::Get_Info_flag() const { return Info_flag; }

/*
string a1 = "(Error-Num)                             (Error-Description)\n";
string a2 = "Distribution.E1				User-specified grid is outside 0 - 1.\n";
string a3 = "Distribution.E2				User-specified distribution is negative.\n";

string b1 = "(Warning-Num)                           (Warning-Description)\n";
string b2 = "Distribution.W1				User-specified grid is not sorted or there are duplicated points.\n";
string b3 = "Distribution.W2				User-specified distribution has more than 1 peak.\n";
*/