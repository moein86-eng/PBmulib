#include "Coalescence_Declaration.h"
#include "Logger_Declaration.h"
#include "MP_Fluid_Declaration.h"


template<typename Grid_type>
Coalescence<Grid_type>::Coalescence(const Grid_type* Grid, const MP_Fluid* Fluid_Val, std::string Formula_str, std::string Variable_str_i, std::string Variable_str_j, std::string VarType_str, const std::string* Parameter_str, const double* Parameter_Val, unsigned int Parameter_number, const std::string* Fluid_str, const std::string* Fluid_ID, unsigned int Fluid_number){
	this->Grid = Grid;
	this->Fluid_Val = Fluid_Val;

	this->Parameter_str = Parameter_str;
	this->Parameter_Val = Parameter_Val;
	this->Parameter_number = Parameter_number;

	this->Fluid_str = Fluid_str;
	this->Fluid_ID = Fluid_ID;
	this->Fluid_number = Fluid_number;

	this->Formula_str = Formula_str;
	this->Variable_str_i = Variable_str_i;
	this->Variable_str_j = Variable_str_j;
	this->VarType_str = VarType_str;

	if ((VarType_str == "radius_DL") || (VarType_str == "radius")) {
		Coalescence::Parser_fun();
	}

	else if ((VarType_str == "volume_DL") || (VarType_str == "volume")) {
		std::string For_str = Coalescence::vol_to_rad();
		Coalescence::Parser_fun(For_str, "ri", "rj");
	}

	else {
		Logger::Instance().Add_Error("Coalescence.E2	     			Specified variable type is unknown.\n");
		//exit(-1);
	}
}

template<typename Grid_type>
Coalescence<Grid_type>::Coalescence(const Grid_type* Grid, const MP_Fluid* Fluid_Val, std::string ModelName_str, const double* Parameter_Val) {
	this->Grid = Grid;
	this->Parameter_Val = Parameter_Val;
	this->Fluid_Val = Fluid_Val;


	const char* ModelName_str_temp = ModelName_str.c_str();
	if (ModelName_str == "Model1") {
		Coalescence::Model1();
	}
	else if (ModelName_str == "Model2") {
		Coalescence::Model2();
	}
	else {
		Logger::Instance().Add_Error("Coalescence.E1	     			Specified Coalescence distribution model is unknown.\n");
		//exit(-1);
	}
}

template<typename Grid_type>
Coalescence<Grid_type>::Coalescence(const Grid_type* Grid, const MP_Fluid* Fluid_Val, double** C_r, double** C_r_Shifted) {
	this->Grid = Grid;
	this->Fluid_Val = Fluid_Val;
	this->C_r = C_r;
	this->C_r_Shifted = C_r_Shifted;
}

template<typename Grid_type>
Coalescence<Grid_type>::~Coalescence() {
	Coalescence::clean_up(C_r);
	Coalescence::clean_up(C_r_Shifted);
}

template<typename Grid_type>
std::string Coalescence<Grid_type>::Get_VarType() const { return VarType_str;}

template<typename Grid_type>
std::string Coalescence<Grid_type>::Get_Variable() const { return Variable_str_i + " & " + Variable_str_j; }

template<typename Grid_type>
std::string Coalescence<Grid_type>::Get_Formula() const { return Formula_str;}

template<typename Grid_type>
unsigned int Coalescence<Grid_type>::Get_Parameter_number() const { return Parameter_number; }

template<typename Grid_type>
std::string Coalescence<Grid_type>::Get_Parameter_str() {
	if (Parameter_str && (Parameter_number > 0)) {
		Parameter_str_sum = " ";
		for (unsigned int i = 0; i < Parameter_number; i++)
			Parameter_str_sum = Parameter_str_sum + " , " + Parameter_str[i];
	}
	return Parameter_str_sum; }

template<typename Grid_type>
const double* Coalescence<Grid_type>::Get_Parameter_Val() const { return Parameter_Val; }

template<typename Grid_type>
unsigned int Coalescence<Grid_type>::Get_Fluid_number() const { return Fluid_number; }

template<typename Grid_type>
std::string Coalescence<Grid_type>::Get_Fluid_str() {
	if(Fluid_str && (Fluid_number > 0)){
		Fluid_str_sum = " ";
		for (unsigned int i = 0; i < Fluid_number;i++)
			Fluid_str_sum = Fluid_str_sum + " , " + Fluid_str[i];
		}
	return Fluid_str_sum; }

template<typename Grid_type>
double* Coalescence<Grid_type>::Get_Fluid_Val() const { return Fluid_Val; }

template<typename Grid_type>
void Coalescence<Grid_type>::clean_up(double** arr) {
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		delete[] arr[i];
	}
	delete[] arr;
}

template<typename Grid_type>
void Coalescence<Grid_type>::Parser_fun() {
	typedef exprtk::symbol_table<double> symbol_table_t;
	typedef exprtk::expression<double>     expression_t;
	typedef exprtk::parser<double>             parser_t;

	const std::string expression_string = Formula_str;

	double xj, xi;

	symbol_table_t symbol_table;
	symbol_table.add_variable(Variable_str_j, xj);
	symbol_table.add_variable(Variable_str_i, xi);

	for (unsigned int i = 0; i < Parameter_number;i++)
		symbol_table.add_constant(Parameter_str[i], Parameter_Val[i]);

	for (unsigned int i = 0; i < Fluid_number; i++) {
		if (Fluid_ID[i] == "Tur_Ener_Diss_Rate_Av") {
			symbol_table.add_constant(Fluid_str[i], Fluid_Val->Tur_Ener_Diss_Rate_Av);
			continue;
		}
		if (Fluid_ID[i] == "Density_Dispersed") {
			symbol_table.add_constant(Fluid_str[i], Fluid_Val->Fluid_D->Density);
			continue;
		}
		if (Fluid_ID[i] == "Density_Continuous") {
			symbol_table.add_constant(Fluid_str[i], Fluid_Val->Fluid_C->Density);
			continue;
		}
		if (Fluid_ID[i] == "Viscosity_Dispersed") {
			symbol_table.add_constant(Fluid_str[i], Fluid_Val->Fluid_D->Viscosity);
			continue;
		}
		if (Fluid_ID[i] == "Viscosity_Continuous") {
			symbol_table.add_constant(Fluid_str[i], Fluid_Val->Fluid_C->Viscosity);
			continue;
		}
		if (Fluid_ID[i] == "Surface_Tension") {
			symbol_table.add_constant(Fluid_str[i], Fluid_Val->Surface_Tension);
			continue;
		}
		Logger::Instance().Add_Error("Coalescence.E3	     			Specified fluid ID (" + Fluid_ID[i] + ") is unknown.\n");
	}

	symbol_table.add_constants();
	expression_t expression;
	expression.register_symbol_table(symbol_table);

	parser_t parser;

	parser.compile(expression_string, expression);
	double** C_r_temp;
	C_r_temp = new double* [Grid->TotalPointsNumber];
	for (unsigned int j = 0; j < Grid->TotalPointsNumber; j++) {
		C_r_temp[j] = new double[Grid->TotalPointsNumber];
		for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
			if ((VarType_str == "volume") || (VarType_str == "radius")){
				xi = Grid->Grid_Points_[i] * Grid->R_max;
				xj = Grid->Grid_Points_[j] * Grid->R_max;
			}
			else {
				xi = Grid->Grid_Points_[i];
				xj = Grid->Grid_Points_[j];
			}
			C_r_temp[j][i] = expression.value();
		}
	}
	this-> C_r = C_r_temp;

	double** C_r_Shifted_temp;
	double Rm = Grid->R_max;
	const double* x = Grid->Grid_Points_;
	C_r_Shifted_temp = new double* [Grid->TotalPointsNumber];
	for (unsigned int j = 0; j < Grid->TotalPointsNumber; j++) {
		C_r_Shifted_temp[j] = new double[Grid->TotalPointsNumber];
		for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
			if ((VarType_str == "volume") || (VarType_str == "radius")) {
				xi = x[j] * x[i] * Rm / pow(2.0, (1.0 / 3.0));
				xj = pow((pow((x[i] * Rm), 3) - pow(xi, 3)), (1.0 / 3.0));
			}
			else {
				xi = x[j] * x[i] / pow(2.0, (1.0 / 3.0));
				xj = pow((pow(x[i], 3) - pow(xi, 3)), (1.0 / 3.0));
			}
			C_r_Shifted_temp[j][i] = expression.value();
		}
	}
	this->C_r_Shifted = C_r_Shifted_temp;
}

template<typename Grid_type>
void Coalescence<Grid_type>::Parser_fun(std::string For_str, std::string rj, std::string ri) {
	typedef exprtk::symbol_table<double> symbol_table_t;
	typedef exprtk::expression<double>     expression_t;
	typedef exprtk::parser<double>             parser_t;

	const std::string expression_string = For_str;

	double xj, xi;

	symbol_table_t symbol_table;
	symbol_table.add_variable(rj, xj);
	symbol_table.add_variable(ri, xi);
	
	for (unsigned int i = 0; i < Parameter_number;i++)
		symbol_table.add_constant(Parameter_str[i], Parameter_Val[i]);

	for (unsigned int i = 0; i < Fluid_number; i++) {
		if (Fluid_ID[i] == "Tur_Ener_Diss_Rate_Av") {
			symbol_table.add_constant(Fluid_str[i], Fluid_Val->Tur_Ener_Diss_Rate_Av);
			continue;
		}
		if (Fluid_ID[i] == "Density_Dispersed") {
			symbol_table.add_constant(Fluid_str[i], Fluid_Val->Fluid_D->Density);
			continue;
		}
		if (Fluid_ID[i] == "Density_Continuous") {
			symbol_table.add_constant(Fluid_str[i], Fluid_Val->Fluid_C->Density);
			continue;
		}
		if (Fluid_ID[i] == "Viscosity_Dispersed") {
			symbol_table.add_constant(Fluid_str[i], Fluid_Val->Fluid_D->Viscosity);
			continue;
		}
		if (Fluid_ID[i] == "Viscosity_Continuous") {
			symbol_table.add_constant(Fluid_str[i], Fluid_Val->Fluid_C->Viscosity);
			continue;
		}
		if (Fluid_ID[i] == "Surface_Tension") {
			symbol_table.add_constant(Fluid_str[i], Fluid_Val->Surface_Tension);
			continue;
		}
		Logger::Instance().Add_Error("Coalescence.E3	     			Specified fluid ID (" + Fluid_ID[i] + ") is unknown.\n");
	}

	symbol_table.add_constants();

	expression_t expression;
	expression.register_symbol_table(symbol_table);

	parser_t parser;
	parser.compile(expression_string, expression);

	double** C_r_temp;
	C_r_temp = new double* [Grid->TotalPointsNumber];
	for (unsigned int j = 0; j < Grid->TotalPointsNumber; j++) {
		C_r_temp[j] = new double[Grid->TotalPointsNumber];
		for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
			if ((VarType_str == "volume") || (VarType_str == "radius")) {
				xi = Grid->Grid_Points_[i] * Grid->R_max;
				xj = Grid->Grid_Points_[j] * Grid->R_max;
			}
			else {
				xi = Grid->Grid_Points_[i];
				xj = Grid->Grid_Points_[j];
			}
			C_r_temp[j][i] = expression.value();
		}
	}
	this->C_r = C_r_temp;

	double** C_r_Shifted_temp;
	double Rm = Grid->R_max;
	const double* x = Grid->Grid_Points_;
	C_r_Shifted_temp = new double* [Grid->TotalPointsNumber];
	for (unsigned int j = 0; j < Grid->TotalPointsNumber; j++) {
		C_r_Shifted_temp[j] = new double[Grid->TotalPointsNumber];
		for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
			if ((VarType_str == "volume") || (VarType_str == "radius")) {
				xi = x[j] * x[i] * Rm / pow(2.0, (1.0 / 3.0));
				xj = pow((pow((x[i] * Rm), 3) - pow(xi, 3)), (1.0 / 3.0));
			}
			else {
				xi = x[j] * x[i] / pow(2.0, (1.0 / 3.0));
				xj = pow((pow(x[i], 3) - pow(xi, 3)), (1.0 / 3.0));
			}
			C_r_Shifted_temp[j][i] = expression.value();
		}
	}
	this->C_r_Shifted = C_r_Shifted_temp;
}

template<typename Grid_type>
void Coalescence<Grid_type>::Model1() {
	this->Formula_str = "4*2^(1/3)*kc1*e^(1/3)*(ri+rj)^2*(ri^(2/3)+rj^(2/3))^0.5*exp(-kc2*ro_c^0.5*e^(1/3)/2^(1/6)/to_wo^0.5*(1/2*(1/ri+1/rj)^-1)^(5/6))";
	this->Variable_str_j = "rj";
	this->Variable_str_i = "ri";
	this->VarType_str = "radius_DL";
	this->Parameter_number = 2;
	this->Parameter_str_sum = "kc1 , kc2";
	this->Fluid_number = 2;
	this->Fluid_str_sum = "ro_c , to_wo, e";

	double ro_c = Fluid_Val->Fluid_C->Density;
	double to_wo = Fluid_Val->Surface_Tension;
	double e = Fluid_Val->Tur_Ener_Diss_Rate_Av;
	double kc1 = Parameter_Val[0];
	double kc2 = Parameter_Val[1];
	double w, r_eq, sai;
	double Rm = Grid->R_max;
	const double* x = Grid->Grid_Points_;

	double** C_r_temp;
	C_r_temp = new double* [Grid->TotalPointsNumber];
	for (unsigned int j = 0; j < Grid->TotalPointsNumber; j++) {
		C_r_temp[j] = new double[Grid->TotalPointsNumber];
		for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
			w = 4 * pow(2, (1.0 / 3.0)) * kc1 * pow(e, (1.0 / 3.0)) * (x[i]* Rm + x[j] * Rm) * (x[i] * Rm + x[j] * Rm) * pow((pow(x[i] * Rm, (2.0 / 3.0)) + pow(x[j] * Rm, (2.0 / 3.0))), 0.5);
			r_eq= 0.5 * pow((1.0 / (x[i] * Rm) + 1.0 / (x[j] * Rm)) , -1);
			sai = exp(-kc2 * pow(ro_c, 0.5) * pow(e, (1.0 / 3.0)) / pow(2, (1.0 / 6.0)) / pow(to_wo, 0.5) * pow(r_eq, (5.0 / 6.0)));
			C_r_temp[j][i] = w * sai;
		}
	}
	this->C_r = C_r_temp;

	double** C_r_Shifted_temp;
	double rr1;
	double rr2;
	C_r_Shifted_temp = new double* [Grid->TotalPointsNumber];
	for (unsigned int j = 0; j < Grid->TotalPointsNumber; j++) {
		C_r_Shifted_temp[j] = new double[Grid->TotalPointsNumber];
		for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
			rr1 = x[j]*x[i]*Rm/pow(2.0,(1.0/3.0));
			rr2 = pow((pow((x[i]*Rm),3) - pow(rr1,3)), (1.0/3.0));
			w = 4 * pow(2, (1.0 / 3.0)) * kc1 * pow(e, (1.0 / 3.0)) * (rr1 + rr2) * (rr1 + rr2) * pow((pow(rr1, (2.0 / 3.0)) + pow(rr2, (2.0 / 3.0))), 0.5);
			r_eq = 0.5 * pow((1.0 / rr1 + 1.0 / rr2), -1);
			sai = exp(-kc2 * pow(ro_c, 0.5) * pow(e, (1.0 / 3.0)) / pow(2, (1.0 / 6.0)) / pow(to_wo, 0.5) * pow(r_eq, (5.0 / 6.0)));
			C_r_Shifted_temp[j][i] = w * sai;
		}
	}
	this->C_r_Shifted = C_r_Shifted_temp;
}

template<typename Grid_type>
void Coalescence<Grid_type>::Model2() {
	this->Formula_str = "const";
	this->Variable_str_j = "rj";
	this->Variable_str_i = "ri";
	this->VarType_str = "radius_DL";
	this->Parameter_number = 1;
	this->Parameter_str_sum = "const";
	this->Fluid_number = 0;
	this->Fluid_str = nullptr;
	this->Fluid_str_sum = "-";

	double** C_r_temp;
	C_r_temp = new double* [Grid->TotalPointsNumber];
	for (unsigned int j = 0; j < Grid->TotalPointsNumber; j++) {
		C_r_temp[j] = new double[Grid->TotalPointsNumber];
		for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
			C_r_temp[j][i] = Parameter_Val[0];
		}
	}
	this->C_r = C_r_temp;

	double** C_r_Shifted_temp;
	C_r_Shifted_temp = new double* [Grid->TotalPointsNumber];
	for (unsigned int j = 0; j < Grid->TotalPointsNumber; j++) {
		C_r_Shifted_temp[j] = new double[Grid->TotalPointsNumber];
		for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
			C_r_Shifted_temp[j][i] = Parameter_Val[0];
		}
	}
	this->C_r_Shifted = C_r_Shifted_temp;
}

template<typename Grid_type>
std::string Coalescence<Grid_type>::vol_to_rad() {
	std::string string1(Formula_str);
	string1 = std::regex_replace(string1, std::regex(Variable_str_j), "(rj^3)");
	string1 = std::regex_replace(string1, std::regex(Variable_str_i), "(ri^3)");
	return string1;
}