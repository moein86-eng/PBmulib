#include "Breakage_Declaration.h"
#include "Logger_Declaration.h"
#include "MP_Fluid_Declaration.h"


template<typename Grid_type>
Breakage<Grid_type>::Breakage(const Grid_type* Grid, const MP_Fluid* Fluid_Val, std::string Formula_str, std::string Variable_str, std::string VarType_str, const std::string* Parameter_str, double* Parameter_Val, unsigned int Parameter_number, const std::string* Fluid_str, const std::string* Fluid_ID, unsigned int Fluid_number) {
	this->Grid = Grid;
	this->Fluid_Val = Fluid_Val;

	this->Formula_str = Formula_str;
	this->Variable_str = Variable_str;
	this->VarType_str = VarType_str;

	this->Parameter_str = Parameter_str;
	this->Parameter_Val = Parameter_Val;
	this->Parameter_number = Parameter_number;

	this->Fluid_str = Fluid_str;
	this->Fluid_ID = Fluid_ID;
	this->Fluid_number = Fluid_number;

	if ((VarType_str == "radius_DL") || (VarType_str == "radius")) {
		Breakage::Parser_fun();
	}

	else if ((VarType_str == "volume_DL") || (VarType_str == "volume")) {
		std::string For_str = Breakage::vol_to_rad();
		Breakage::Parser_fun(For_str, "r");
	}

	else {
		Logger::Instance().Add_Error("Breakage.E2	     			Specified variable type is unknown.\n");
		//exit(-1);
	}
}

template<typename Grid_type>
Breakage<Grid_type>::Breakage(const Grid_type* Grid, const MP_Fluid* Fluid_Val, std::string ModelName_str, double* Parameter_Val) {
	this->Grid = Grid;
	this->Parameter_Val = Parameter_Val;
	this->Fluid_Val = Fluid_Val;

	const char* ModelName_str_temp = ModelName_str.c_str();
	if (ModelName_str == "Model1") {
		Breakage::Model1();
	}
	else if (ModelName_str == "Model2") {
		Breakage::Model2();
	}
	else {
		Logger::Instance().Add_Error("Breakage.E1	     			Specified Breakage distribution model is unknown.\n");
		//exit(-1);
	}
}


template<typename Grid_type>
Breakage<Grid_type>::Breakage(const Grid_type* Grid, const MP_Fluid* Fluid_Val, double* B_f, double** B_f_Shifted) {
	this->Grid = Grid;
	this->Fluid_Val = Fluid_Val;
	this->B_f = B_f;
	this->B_f_Shifted = B_f_Shifted;
}

template<typename Grid_type>
Breakage<Grid_type>::~Breakage() {
	Breakage::clean_up(B_f);
	Breakage::clean_up(B_f_Shifted);
}

template<typename Grid_type>
std::string Breakage<Grid_type>::Get_VarType() const { return VarType_str; }

template<typename Grid_type>
std::string Breakage<Grid_type>::Get_Variable() const { return Variable_str; }

template<typename Grid_type>
std::string Breakage<Grid_type>::Get_Formula() const { return Formula_str; }

template<typename Grid_type>
unsigned int Breakage<Grid_type>::Get_Parameter_number() const { return Parameter_number; }

template<typename Grid_type>
std::string Breakage<Grid_type>::Get_Parameter_str() {
	if (Parameter_str && (Parameter_number>0)) {
		Parameter_str_sum = Parameter_str[0];
		for (unsigned int i = 1; i < Parameter_number; i++)
			Parameter_str_sum = Parameter_str_sum + " , " + Parameter_str[i];
	}
	return Parameter_str_sum;
}

template<typename Grid_type>
double* Breakage<Grid_type>::Get_Parameter_Val() const { return Parameter_Val; }

template<typename Grid_type>
unsigned int Breakage<Grid_type>::Get_Fluid_number() const { return Fluid_number; }

template<typename Grid_type>
std::string Breakage<Grid_type>::Get_Fluid_str() {
	if (Fluid_str && (Fluid_number>0)) {
		Fluid_str_sum = Fluid_str[0];
		for (unsigned int i = 1; i < Fluid_number; i++)
			Fluid_str_sum = Fluid_str_sum + " , " + Fluid_str[i];
	}
	return Fluid_str_sum;
}

template<typename Grid_type>
double* Breakage<Grid_type>::Get_Fluid_Val() const { return Fluid_Val; }

template<typename Grid_type>
void Breakage<Grid_type>::clean_up(double* arr) {
	delete[] arr;
}

template<typename Grid_type>
void Breakage<Grid_type>::clean_up(double** arr) {
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		delete[] arr[i];
	}
	delete[] arr;
}

template<typename Grid_type>
void Breakage<Grid_type>::Parser_fun() {
	typedef exprtk::symbol_table<double> symbol_table_t;
	typedef exprtk::expression<double>     expression_t;
	typedef exprtk::parser<double>             parser_t;

	const std::string expression_string = Formula_str;

	double x;

	symbol_table_t symbol_table;
	symbol_table.add_variable(Variable_str, x);

	for(unsigned int i=0; i< Parameter_number;i++)
		symbol_table.add_constant(Parameter_str[i], Parameter_Val[i]);

	for (unsigned int i = 0; i < Fluid_number; i++){
			if(Fluid_ID[i] == "Tur_Ener_Diss_Rate_Av"){
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
			Logger::Instance().Add_Error("Breakage.E3	     			Specified fluid ID (" + Fluid_ID[i]  + ") is unknown.\n");
	}

	symbol_table.add_constants();
	expression_t expression;
	expression.register_symbol_table(symbol_table);

	parser_t parser;
	parser.compile(expression_string, expression);
	
	double* B_f_temp = new double [Grid->TotalPointsNumber];
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		if((VarType_str=="volume") || (VarType_str == "radius"))
			x = Grid->Grid_Points_[i] * Grid->R_max;
		else
			x = Grid->Grid_Points_[i];
    
		B_f_temp[i] = expression.value();
		}
	B_f_temp[0] = 0;
	this->B_f = B_f_temp;

	double Rm = Grid->R_max;
	const double* xx = Grid->Grid_Points_;
	double** B_f_Shifted_temp = new double* [Grid->TotalPointsNumber];
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		B_f_Shifted_temp[i] = new double[Grid->TotalPointsNumber];
		for (unsigned int j = 0; j < Grid->TotalPointsNumber; j++) {
			if ((VarType_str == "volume") || (VarType_str == "radius"))
				x = xx[j] * (Rm-xx[i] * Rm) + xx[i]*Rm;
			else
				x = xx[j] * (1 - xx[i]) + xx[i];

			B_f_Shifted_temp[i][j] = expression.value();
		}
	}
	B_f_Shifted_temp[0][0] = 0;
	this->B_f_Shifted = B_f_Shifted_temp;

}

template<typename Grid_type>
void Breakage<Grid_type>::Parser_fun(std::string For_str, std::string r) {
	typedef exprtk::symbol_table<double> symbol_table_t;
	typedef exprtk::expression<double>     expression_t;
	typedef exprtk::parser<double>             parser_t;

	const std::string expression_string = For_str;

	double x;

	symbol_table_t symbol_table;
	symbol_table.add_variable(r, x);
	
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
		Logger::Instance().Add_Error("Breakage.E3	     			Specified fluid ID (" + Fluid_ID[i] + ") is unknown.\n");
	}

	symbol_table.add_constants();

	expression_t expression;
	expression.register_symbol_table(symbol_table);

	parser_t parser;
	parser.compile(expression_string, expression);
	double* B_f_temp = new double [Grid->TotalPointsNumber];
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		if ((VarType_str == "volume") || (VarType_str == "radius"))
			x = Grid->Grid_Points_[i] * Grid->R_max;
		else
			x = Grid->Grid_Points_[i];

		B_f_temp[i] = expression.value();
	}
	B_f_temp[0] = 0;
	this->B_f = B_f_temp;

	double Rm = Grid->R_max;
	const double* xx = Grid->Grid_Points_;
	double** B_f_Shifted_temp = new double* [Grid->TotalPointsNumber];
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		B_f_Shifted_temp[i] = new double[Grid->TotalPointsNumber];
		for (unsigned int j = 0; j < Grid->TotalPointsNumber; j++) {
			if ((VarType_str == "volume") || (VarType_str == "radius"))
				x = xx[j] * (Rm - xx[i] * Rm) + xx[i] * Rm;
			else
				x = xx[j] * (1 - xx[i]) + xx[i];

			B_f_Shifted_temp[i][j] = expression.value();
		}
	}
	B_f_Shifted_temp[0][0] = 0;
	this->B_f_Shifted = B_f_Shifted_temp;
}

template<typename Grid_type>
void Breakage<Grid_type>::Model1() {
	this->Formula_str = "kb1*e^(1/3)/2^(2/3)/r^(2/3)*(ro_d/ro_c)^0.5*exp(-kb2*to_wo/ro_d/2^(5/3)/e^(2/3)/r^(5/3))";
	this->Variable_str = "r";
	this->VarType_str = "radius_DL";
	this->Parameter_number = 2;
	this->Parameter_str_sum = "kb1 , kb2";
	this->Fluid_number = 4;
	this->Fluid_str_sum = "ro_d , ro_c, to_wo, e";

	double ro_c = Fluid_Val->Fluid_C->Density;
	double ro_d = Fluid_Val->Fluid_D->Density;
	double to_wo = Fluid_Val->Surface_Tension;
	double e = Fluid_Val->Tur_Ener_Diss_Rate_Av;
	double kb1 = Parameter_Val[0];
	double kb2 = Parameter_Val[1];
	double Rm = Grid->R_max;
	const double* rr = Grid->Grid_Points_;
	double rrj;

	double* B_f_temp = new double [Grid->TotalPointsNumber];
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		B_f_temp[i] = kb1 * pow(e , (1.0 / 3.0)) / pow(2.0 , (2.0 / 3.0)) / pow(rr[i]*Rm , (2.0 / 3.0)) * pow((ro_d / ro_c) , 0.5) * exp(-kb2 * to_wo / ro_d / pow(2.0 , (5.0 / 3.0)) / pow(e , (2.0 / 3.0)) / pow(rr[i] *Rm , (5.0 / 3.0)));
	}
	B_f_temp[0] = 0;
	this->B_f = B_f_temp;

	double** B_f_Shifted_temp = new double* [Grid->TotalPointsNumber];
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		B_f_Shifted_temp[i] = new double[Grid->TotalPointsNumber];
		for (unsigned int j = 0; j < Grid->TotalPointsNumber; j++) {
			rrj = rr[j] * (Rm - rr[i] * Rm) + rr[i] * Rm;
			B_f_Shifted_temp[i][j] = kb1 * pow(e, (1.0 / 3.0)) / pow(2.0, (2.0 / 3.0)) / pow(rrj, (2.0 / 3.0)) * pow((ro_d / ro_c), 0.5) * exp(-kb2 * to_wo / ro_d / pow(2.0, (5.0 / 3.0)) / pow(e, (2.0 / 3.0)) / pow(rrj, (5.0 / 3.0)));;
		}
	}
	B_f_Shifted_temp[0][0] = 0;
	this->B_f_Shifted = B_f_Shifted_temp;
}

template<typename Grid_type>
void Breakage<Grid_type>::Model2() {
	this->Formula_str = "const";
	this->Variable_str = "r";
	this->VarType_str = "radius_DL";
	this->Parameter_number = 1;
	this->Parameter_str_sum = "const";
	this->Fluid_number = 0;
	this->Fluid_str_sum = " ";

	double* B_f_temp = new double [Grid->TotalPointsNumber];
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		B_f_temp[i] = Parameter_Val[0];
	}
	this->B_f = B_f_temp;

	double** B_f_Shifted_temp = new double* [Grid->TotalPointsNumber];
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		B_f_Shifted_temp[i] = new double[Grid->TotalPointsNumber];
		for (unsigned int j = 0; j < Grid->TotalPointsNumber; j++) {
			B_f_Shifted_temp[i][j] = Parameter_Val[0];
		}
	}
	this->B_f_Shifted = B_f_Shifted_temp;
}

template<typename Grid_type>
std::string Breakage<Grid_type>::vol_to_rad() {
	std::string string1(Formula_str);
	string1 = std::regex_replace(string1, std::regex(Variable_str), "(rj^3)");
	return string1;
}