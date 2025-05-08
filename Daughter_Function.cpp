#include "Daughter_Declaration.h"
#include "Logger_Declaration.h"


template<typename Grid_type>
Daughter<Grid_type>::Daughter(const Grid_type* Grid, std::string Formula_str, std::string Variable_str_j, std::string Variable_str_i, std::string VarType_str) {

	this->Grid = Grid;

	this->Formula_str = Formula_str;
	this->Variable_str_j = Variable_str_j;
	this->Variable_str_i = Variable_str_i;
	this->VarType_str = VarType_str;

	if (VarType_str == "radius_DL") {
		Daughter::Parser_fun();
		Daughter::InputChecking();
	}

	else if (VarType_str == "volume_DL") {
		std::string For_str = Daughter::vol_to_rad();
		Daughter::Parser_fun(For_str, "rj", "ri");
		Daughter::InputChecking();
	}

	else {
		Logger::Instance().Add_Error("Daughter.E2	     			Specified variable type is unknown.\n");
		//exit(-1);
	}
}

template<typename Grid_type>
Daughter<Grid_type>::Daughter(Grid_type* Grid, std::string ModelName_str) {

	this->Grid = Grid;

	const char* ModelName_str_temp = ModelName_str.c_str();
	if (ModelName_str == "Model1") {
		Daughter::Model1();
	}
	else if (ModelName_str == "Model2") {
		Daughter::Model2();
	}
	else {
		Logger::Instance().Add_Error("Daughter.E1	     			Specified daughter distribution model is unknown.\n");
		//exit(-1);
	}
	Daughter::InputChecking();
}

template<typename Grid_type>
Daughter<Grid_type>::Daughter(const Grid_type* Grid, double** beta, double** beta_Shifted) {
	this->Grid = Grid;
	this->beta = beta;
	this->beta_Shifted = beta_Shifted;
}

template<typename Grid_type>
Daughter<Grid_type>::~Daughter() {
	Daughter::clean_up(beta);
	Daughter::clean_up(beta_Shifted);
}
template<typename Grid_type>
std::string Daughter<Grid_type>::Get_VarType() const { return VarType_str; }

template<typename Grid_type>
std::string Daughter<Grid_type>::Get_Variable() const { return Variable_str_j + " & " + Variable_str_i; }

template<typename Grid_type>
std::string Daughter<Grid_type>::Get_Formula() const { return Formula_str; }

template<typename Grid_type>
void Daughter<Grid_type>::clean_up(double** arr) {
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		delete[] arr[i];
	}
	delete[] arr;
}

template<typename Grid_type>
void Daughter<Grid_type>::InputChecking() {
	double sum = Daughter::Get_Int();
	if (abs(sum - 1) > 0.01) {
		auto w1a = std::to_string(sum);
		std::string w1 = "Daughter.W1	     			 Daughter distribution integration over domain is " + w1a + ".\n";
		Logger::Instance().Add_Warning(w1);
	}
}

template<typename Grid_type>
void Daughter<Grid_type>::Parser_fun() {
	typedef exprtk::symbol_table<double> symbol_table_t;
	typedef exprtk::expression<double>     expression_t;
	typedef exprtk::parser<double>             parser_t;

	const std::string expression_string = Formula_str;

	double xj, xi;

	symbol_table_t symbol_table;
	symbol_table.add_variable(Variable_str_j, xj);
	symbol_table.add_variable(Variable_str_i, xi);
	symbol_table.add_constants();

	expression_t expression;
	expression.register_symbol_table(symbol_table);

	parser_t parser;
	parser.compile(expression_string, expression);
	double** beta_temp;
	beta_temp = new double* [Grid->TotalPointsNumber];
	for (unsigned int j = 0; j < Grid->TotalPointsNumber; j++) {
		beta_temp[j] = new double[Grid->TotalPointsNumber];
		std::fill(beta_temp[j] + j, beta_temp[j] + Grid->TotalPointsNumber, 0);
		for (unsigned int i = 0; i < j; i++) {
			xj = Grid->Grid_Points_[j];
			xi = Grid->Grid_Points_[i];
			beta_temp[j][i] = expression.value();
		}
	}
	this->beta = beta_temp;

	double** beta_Shifted_temp;
	const double* xx = Grid->Grid_Points_;
	beta_Shifted_temp = new double* [Grid->TotalPointsNumber];
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		beta_Shifted_temp[i] = new double[Grid->TotalPointsNumber];
		for (unsigned int j = 0; j < Grid->TotalPointsNumber; j++) {
			xj = xx[j] * (1 - xx[i]) + xx[i];
			xi = xx[i];
			beta_Shifted_temp[i][j] = expression.value();
		}
	}
	beta_Shifted_temp[0][0] = 0;
	this->beta_Shifted = beta_Shifted_temp;
}

template<typename Grid_type>
void Daughter<Grid_type>::Parser_fun(std::string For_str, std::string rj, std::string ri) {
	typedef exprtk::symbol_table<double> symbol_table_t;
	typedef exprtk::expression<double>     expression_t;
	typedef exprtk::parser<double>             parser_t;

	const std::string expression_string = For_str;

	double xj, xi;

	symbol_table_t symbol_table;
	symbol_table.add_variable(rj, xj);
	symbol_table.add_variable(ri, xi);
	symbol_table.add_constants();

	expression_t expression;
	expression.register_symbol_table(symbol_table);

	parser_t parser;
	parser.compile(expression_string, expression);
	double** beta_temp;
	beta_temp = new double* [Grid->TotalPointsNumber];
	for (unsigned int j = 0; j < Grid->TotalPointsNumber; j++) {
		beta_temp[j] = new double[Grid->TotalPointsNumber];
		std::fill(beta_temp[j] + j, beta_temp[j] + Grid->TotalPointsNumber, 0);
		for (unsigned int i = 0; i < j; i++) {
			xj = Grid->Grid_Points_[j];
			xi = Grid->Grid_Points_[i];
			beta_temp[j][i] = expression.value();
		}
	}
	this->beta = beta_temp;

	double** beta_Shifted_temp;
	const double* xx = Grid->Grid_Points_;
	beta_Shifted_temp = new double* [Grid->TotalPointsNumber];
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		beta_Shifted_temp[i] = new double[Grid->TotalPointsNumber];
		for (unsigned int j = 0; j < Grid->TotalPointsNumber; j++) {
			xj = xx[j] * (1 - xx[i]) + xx[i];
			xi = xx[i];
			beta_Shifted_temp[i][j] = expression.value();
		}
	}
	beta_Shifted_temp[0][0] = 0;
	this->beta_Shifted = beta_Shifted_temp;
}

template<typename Grid_type>
void Daughter<Grid_type>::Model1() {
	this->Formula_str = "7.2/rj^3*exp(-4.5*(2*ri^3-rj^3)^2/rj^6)*ri^2";
	this->Variable_str_j = "rj";
	this->Variable_str_i = "ri";
	this->VarType_str = "radius_DL";

	double** beta_temp;
	beta_temp = new double* [Grid->TotalPointsNumber];
	for (unsigned int j = 0; j < Grid->TotalPointsNumber; j++) {
		beta_temp[j] = new double[Grid->TotalPointsNumber];
		std::fill(beta_temp[j] + j, beta_temp[j] + Grid->TotalPointsNumber, 0);
		for (unsigned int i = 0; i < j; i++) {
			beta_temp[j][i] = 7.2 * exp(-4.5 * pow((2 * pow(Grid->Grid_Points_[i], 3) - pow(Grid->Grid_Points_[j], 3)), 2) / pow(Grid->Grid_Points_[j], 6)) * pow(Grid->Grid_Points_[i], 2) / pow(Grid->Grid_Points_[j], 3);
		}
	}
	this->beta = beta_temp;

	double xj, xi;
	double** beta_Shifted_temp;
	const double* xx = Grid->Grid_Points_;
	beta_Shifted_temp = new double* [Grid->TotalPointsNumber];
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		beta_Shifted_temp[i] = new double[Grid->TotalPointsNumber];
		for (unsigned int j = 0; j < Grid->TotalPointsNumber; j++) {
			xj = xx[j] * (1 - xx[i]) + xx[i];
			xi = xx[i];
			beta_Shifted_temp[i][j] = 7.2 * exp(-4.5 * pow((2 * pow(xi, 3) - pow(xj, 3)), 2) / pow(xj, 6)) * pow(xi, 2) / pow(xj, 3);
		}
	}
	beta_Shifted_temp[0][0] = 0;
	this->beta_Shifted = beta_Shifted_temp;
}

template<typename Grid_type>
void Daughter<Grid_type>::Model2() {
	this->Formula_str = "3*ri^2/rj^3";
	this->Variable_str_j = "rj";
	this->Variable_str_i = "ri";
	this->VarType_str = "radius_DL";

	double** beta_temp;
	beta_temp = new double* [Grid->TotalPointsNumber];
	for (unsigned int j = 0; j < Grid->TotalPointsNumber; j++) {
		beta_temp[j] = new double[Grid->TotalPointsNumber];
		std::fill(beta_temp[j] + j, beta_temp[j] + Grid->TotalPointsNumber, 0);
		for (unsigned int i = 0; i < j; i++) {
			beta_temp[j][i] = 3 * pow(Grid->Grid_Points_[i], 2) / pow(Grid->Grid_Points_[j], 3);
		}
	}
	this->beta = beta_temp;

	double xj, xi;
	double** beta_Shifted_temp;
	const double* xx = Grid->Grid_Points_;
	beta_Shifted_temp = new double* [Grid->TotalPointsNumber];
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		beta_Shifted_temp[i] = new double[Grid->TotalPointsNumber];
		for (unsigned int j = 0; j < Grid->TotalPointsNumber; j++) {
			xj = xx[j] * (1 - xx[i]) + xx[i];
			xi = xx[i];
			beta_Shifted_temp[i][j] = 3 * pow(xi, 2) / pow(xj, 3);
		}
	}
	beta_Shifted_temp[0][0] = 0;
	this->beta_Shifted = beta_Shifted_temp;
}

template<typename Grid_type>
double Daughter<Grid_type>::Get_Int() {
	double sum = 0;
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++)
		sum += beta[Grid->TotalPointsNumber - 1][i] * Grid->Quadrature_weight_[i];
	return sum;
}

template<typename Grid_type>
std::string Daughter<Grid_type>::vol_to_rad() {
	std::string string1(Formula_str);
	string1 = std::regex_replace(string1, std::regex(Variable_str_j), "(rj^3)");
	string1 = std::regex_replace(string1, std::regex(Variable_str_i), "(ri^3)");
	string1 = "(" + string1 + ")" + "*(3*ri^2)";
	return string1;
}