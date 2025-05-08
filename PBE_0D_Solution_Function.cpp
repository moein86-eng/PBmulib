#pragma warning(disable: 4996)
#include "PBE_0D_Solution_Declaration.h"
#include "Logger_Declaration.h"

using namespace std;
using namespace boost::numeric::odeint;

#define pi_const 3.141592653589793238462643383279502

template<typename Grid_Type>
PBE_0D_Solution<Grid_Type>::PBE_0D_Solution(Distribution<Grid_Type>* Int_Dist, Daughter<Grid_Type>* Daughter_, Breakage<Grid_Type>* Breakage_, Coalescence<Grid_Type>* Coalescence_) {
	auto t1 = std::chrono::high_resolution_clock::now();
	this->Int_Dist = Int_Dist;
	this->Daughter_ = Daughter_;
	this->Breakage_ = Breakage_;
	this->Coalescence_ = Coalescence_;
	this->Flag_Module = 1;
	this->Flag_BC = 1;
	this->Index = 0;
	this->Grid_Size = Coalescence_->Grid->TotalPointsNumber;
	Coalescence_Death_Birth_Coef();
	Breakage_Death_Birth_Coef();
	//ODE_Set_BC(Int_Dist->dist_sai_RB);
	this->ODE_Set_fun = (&PBE_0D_Solution<Grid_Type>::ODE_Set_BC);
	//(this->*ODE_Set_fun)(Int_Dist->dist_sai_RB);
	auto t2 = std::chrono::high_resolution_clock::now();
	this->PrePros_Time = (double)(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()/1000);
}

template<typename Grid_Type>
PBE_0D_Solution<Grid_Type>::PBE_0D_Solution(Distribution<Grid_Type>* Int_Dist, Distribution<Grid_Type>* Inlet_Dist, Daughter<Grid_Type>* Daughter_, Breakage<Grid_Type>* Breakage_, Coalescence<Grid_Type>* Coalescence_, double Residence_Time) {
	auto t1 = std::chrono::high_resolution_clock::now();
	this->Int_Dist = Int_Dist;
	this->Inlet_Dist = Inlet_Dist;
	this->Daughter_ = Daughter_;
	this->Breakage_ = Breakage_;
	this->Coalescence_ = Coalescence_;
	this->Residence_Time = Residence_Time;
	this->Time_DL = Breakage_->Grid->t_c / Residence_Time;
	this->Flag_Module = 2;
	this->Flag_BC = 1;
	this->Index = 0;
	this->Grid_Size = Coalescence_->Grid->TotalPointsNumber;
	Coalescence_Death_Birth_Coef();
	Breakage_Death_Birth_Coef();
	//ODE_Set_BC_IO(Int_Dist->dist_sai_RB);
	this->ODE_Set_fun = (&PBE_0D_Solution<Grid_Type>::ODE_Set_BC_IO);
	//(this->*ODE_Set_fun)(Int_Dist->dist_sai_RB);
	auto t2 = std::chrono::high_resolution_clock::now();
	this->PrePros_Time = (double)(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000);
}

template<typename Grid_Type>
PBE_0D_Solution<Grid_Type>::PBE_0D_Solution(Distribution<Grid_Type>* Int_Dist, Daughter<Grid_Type>* Daughter_, Breakage<Grid_Type>* Breakage_, Coalescence<Grid_Type>* Coalescence_, double Velocity) {
	auto t1 = std::chrono::high_resolution_clock::now();
	this->Int_Dist = Int_Dist;
	this->Daughter_ = Daughter_;
	this->Breakage_ = Breakage_;
	this->Coalescence_ = Coalescence_;
	this->Velocity = Velocity;
	this->Flag_Module = 3;
	this->Flag_BC = 1;
	this->Index = 0;
	this->Grid_Size = Coalescence_->Grid->TotalPointsNumber;
	Coalescence_Death_Birth_Coef();
	Breakage_Death_Birth_Coef();
	//ODE_Set_BC(Int_Dist->dist_sai_RB);
	this->ODE_Set_fun = (&PBE_0D_Solution<Grid_Type>::ODE_Set_BC);
	//(this->*ODE_Set_fun)(Int_Dist->dist_sai_RB);
	auto t2 = std::chrono::high_resolution_clock::now();
	this->PrePros_Time = (double)(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000);
}

template<typename Grid_Type>
PBE_0D_Solution<Grid_Type>::PBE_0D_Solution(Distribution<Grid_Type>* Int_Dist, Coalescence<Grid_Type>* Coalescence_) {
	auto t1 = std::chrono::high_resolution_clock::now();
	this->Int_Dist = Int_Dist;
	this->Coalescence_ = Coalescence_;
	this->Flag_Module = 1;
	this->Flag_BC = 2;
	this->Index = 0;
	this->Grid_Size = Coalescence_->Grid->TotalPointsNumber;
	Coalescence_Death_Birth_Coef();
	//ODE_Set_C(Int_Dist->dist_sai_RB);
	this->ODE_Set_fun = (&PBE_0D_Solution<Grid_Type>::ODE_Set_C);
	//(this->*ODE_Set_fun)(Int_Dist->dist_sai_RB);
	auto t2 = std::chrono::high_resolution_clock::now();
	this->PrePros_Time = (double)(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000);
}

template<typename Grid_Type>
PBE_0D_Solution<Grid_Type>::PBE_0D_Solution(Distribution<Grid_Type>* Int_Dist, Distribution<Grid_Type>* Inlet_Dist, Coalescence<Grid_Type>* Coalescence_, double Residence_Time) {
	auto t1 = std::chrono::high_resolution_clock::now();
	this->Int_Dist = Int_Dist;
	this->Inlet_Dist = Inlet_Dist;
	this->Coalescence_ = Coalescence_;
	this->Residence_Time = Residence_Time;
	this->Time_DL = Coalescence_->Grid->t_c / Residence_Time;
	this->Flag_Module = 2;
	this->Flag_BC = 2;
	this->Index = 0;
	this->Grid_Size = Coalescence_->Grid->TotalPointsNumber;
	Coalescence_Death_Birth_Coef();
	//ODE_Set_C_IO(Int_Dist->dist_sai_RB);
	this->ODE_Set_fun = (&PBE_0D_Solution<Grid_Type>::ODE_Set_C_IO);
	//(this->*ODE_Set_fun)(Int_Dist->dist_sai_RB);
	auto t2 = std::chrono::high_resolution_clock::now();
	this->PrePros_Time = (double)(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000);
}

template<typename Grid_Type>
PBE_0D_Solution<Grid_Type>::PBE_0D_Solution(Distribution<Grid_Type>* Int_Dist, Coalescence<Grid_Type>* Coalescence_, double Velocity) {
	auto t1 = std::chrono::high_resolution_clock::now();
	this->Int_Dist = Int_Dist;
	this->Coalescence_ = Coalescence_;
	this->Velocity = Velocity;
	this->Flag_Module = 3;
	this->Flag_BC = 2;
	this->Index = 0;
	this->Grid_Size = Coalescence_->Grid->TotalPointsNumber;
	Coalescence_Death_Birth_Coef();
	//ODE_Set_C(Int_Dist->dist_sai_RB);
	this->ODE_Set_fun = (&PBE_0D_Solution<Grid_Type>::ODE_Set_C);
	//(this->*ODE_Set_fun)(Int_Dist->dist_sai_RB);
	auto t2 = std::chrono::high_resolution_clock::now();
	this->PrePros_Time = (double)(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000);
}

template<typename Grid_Type>
PBE_0D_Solution<Grid_Type>::PBE_0D_Solution(Distribution<Grid_Type>* Int_Dist, Daughter<Grid_Type>* Daughter_, Breakage<Grid_Type>* Breakage_) {
	auto t1 = std::chrono::high_resolution_clock::now();
	this->Int_Dist = Int_Dist;
	this->Daughter_ = Daughter_;
	this->Breakage_ = Breakage_;
	this->Flag_Module = 1;
	this->Flag_BC = 3;
	this->Index = 0;
	this->Grid_Size = Breakage_->Grid->TotalPointsNumber;
	Breakage_Death_Birth_Coef();
	//ODE_Set_B(Int_Dist->dist_sai_RB);
	this->ODE_Set_fun = (&PBE_0D_Solution<Grid_Type>::ODE_Set_B);
	//(this->*ODE_Set_fun)(Int_Dist->dist_sai_RB);
	auto t2 = std::chrono::high_resolution_clock::now();
	this->PrePros_Time = (double)(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000);
}

template<typename Grid_Type>
PBE_0D_Solution<Grid_Type>::PBE_0D_Solution(Distribution<Grid_Type>* Int_Dist, Distribution<Grid_Type>* Inlet_Dist, Daughter<Grid_Type>* Daughter_, Breakage<Grid_Type>* Breakage_, double Residence_Time) {
	auto t1 = std::chrono::high_resolution_clock::now();
	this->Int_Dist = Int_Dist;
	this->Inlet_Dist = Inlet_Dist;
	this->Daughter_ = Daughter_;
	this->Breakage_ = Breakage_;
	this->Residence_Time = Residence_Time;
	this->Time_DL = Breakage_->Grid->t_c / Residence_Time;
	this->Flag_Module = 2;
	this->Flag_BC = 3;
	this->Index = 0;
	this->Grid_Size = Breakage_->Grid->TotalPointsNumber;
	Breakage_Death_Birth_Coef();
	//ODE_Set_B_IO(Int_Dist->dist_sai_RB);
	this->ODE_Set_fun = (&PBE_0D_Solution<Grid_Type>::ODE_Set_B_IO);
	//(this->*ODE_Set_fun)(Int_Dist->dist_sai_RB);
	auto t2 = std::chrono::high_resolution_clock::now();
	this->PrePros_Time = (double)(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000);
}

template<typename Grid_Type>
PBE_0D_Solution<Grid_Type>::PBE_0D_Solution(Distribution<Grid_Type>* Int_Dist, Daughter<Grid_Type>* Daughter_, Breakage<Grid_Type>* Breakage_, double Velocity) {
	auto t1 = std::chrono::high_resolution_clock::now();
	this->Int_Dist = Int_Dist;
	this->Daughter_ = Daughter_;
	this->Breakage_ = Breakage_;
	this->Velocity = Velocity;
	this->Flag_Module = 3;
	this->Flag_BC = 3;
	this->Index = 0;
	this->Grid_Size = Breakage_->Grid->TotalPointsNumber;
	Breakage_Death_Birth_Coef();
	//ODE_Set_B(Int_Dist->dist_sai_RB);
	this->ODE_Set_fun = (&PBE_0D_Solution<Grid_Type>::ODE_Set_B);
	//(this->*ODE_Set_fun)(Int_Dist->dist_sai_RB);
	auto t2 = std::chrono::high_resolution_clock::now();
	this->PrePros_Time = (double)(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000);
}

template<typename Grid_Type>
PBE_0D_Solution<Grid_Type>::PBE_0D_Solution(Distribution<Grid_Type>* Int_Dist, Distribution<Grid_Type>* Inlet_Dist, Grid_Type* Grid_NBC, double Residence_Time) {
	auto t1 = std::chrono::high_resolution_clock::now();
	this->Int_Dist = Int_Dist;
	this->Inlet_Dist = Inlet_Dist;
	this->Residence_Time = Residence_Time;
	this->Grid_NBC = Grid_NBC;
	this->Time_DL = Grid_NBC->t_c / Residence_Time;
	this->Flag_Module = 2;
	this->Flag_BC = 0;
	this->Index = 0;
	this->Grid_Size = Grid_NBC->TotalPointsNumber;
	//ODE_Set_NBC_IO(Int_Dist->dist_sai_RB);
	this->ODE_Set_fun = (&PBE_0D_Solution<Grid_Type>::ODE_Set_NBC_IO);
	//(this->*ODE_Set_fun)(Int_Dist->dist_sai_RB);
	auto t2 = std::chrono::high_resolution_clock::now();
	this->PrePros_Time = (double)(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000);
}

template<typename Grid_Type>
PBE_0D_Solution<Grid_Type>::~PBE_0D_Solution() {

	if (Flag_BC == 1) {   //both breakage and coalescence
		PBE_0D_Solution::clean_up(Br_B_Coef);
		PBE_0D_Solution::clean_up(Br_D_Coef);
		PBE_0D_Solution::clean_up(X_Mat_3);
		PBE_0D_Solution::clean_up(sai3);
		PBE_0D_Solution::clean_up(Co_B_Coef);
		PBE_0D_Solution::clean_up(Co_D_Coef);
		PBE_0D_Solution::clean_up(X_Mat_1);
		PBE_0D_Solution::clean_up(X_Mat_2);
		PBE_0D_Solution::clean_up(sai1);
		PBE_0D_Solution::clean_up(sai2);
		PBE_0D_Solution::clean_up(dy_dx);
		PBE_0D_Solution::clean_up(m);
		PBE_0D_Solution::clean_up(c_coef);
		PBE_0D_Solution::clean_up(d_coef);
	}
	else if (Flag_BC == 2) { //only coalescence
		PBE_0D_Solution::clean_up(Co_B_Coef);
		PBE_0D_Solution::clean_up(Co_D_Coef);
		PBE_0D_Solution::clean_up(X_Mat_1);
		PBE_0D_Solution::clean_up(X_Mat_2);
		PBE_0D_Solution::clean_up(sai1);
		PBE_0D_Solution::clean_up(sai2);
		PBE_0D_Solution::clean_up(dy_dx);
		PBE_0D_Solution::clean_up(m);
		PBE_0D_Solution::clean_up(c_coef);
		PBE_0D_Solution::clean_up(d_coef);
	}
	else if (Flag_BC == 3)   //only breakage
	{
		PBE_0D_Solution::clean_up(Br_B_Coef);
		PBE_0D_Solution::clean_up(Br_D_Coef);
		PBE_0D_Solution::clean_up(X_Mat_3);
		PBE_0D_Solution::clean_up(sai3);
		PBE_0D_Solution::clean_up(dy_dx);
		PBE_0D_Solution::clean_up(m);
		PBE_0D_Solution::clean_up(c_coef);
		PBE_0D_Solution::clean_up(d_coef);
	}
}

template<typename Grid_Type>
void PBE_0D_Solution<Grid_Type>::clean_up(double** arr) {
	for (unsigned int i = 0; i < Grid_Size; i++) {
		delete[] arr[i];
	}
	delete[] arr;
}

template<typename Grid_Type>
void PBE_0D_Solution<Grid_Type>::clean_up(double arr[]) {
	delete[] arr;
}

template<typename Grid_Type>
void PBE_0D_Solution<Grid_Type>::Coalescence_Death_Birth_Coef() {
	double** C_r = Coalescence_->C_r;
	double** C_r_s = Coalescence_->C_r_Shifted;
	double Rm = Coalescence_->Grid->R_max;
	double Vm = 4.0 / 3.0 * pi_const * pow(Rm, 3);
	double t_c = Coalescence_->Grid->t_c;
	double* x = Coalescence_->Grid->Grid_Points_;
	double* W = Coalescence_->Grid->Quadrature_weight_;
	unsigned int I = Coalescence_->Grid->TotalPointsNumber;

	double** X_Mat_1_temp = new double* [I];
	for (unsigned int i = 0; i < I; i++) {
		X_Mat_1_temp[i] = new double[I];
		for (unsigned int j = 0; j < I; j++) {
			X_Mat_1_temp[i][j] = x[i] * x[j] / pow(2.0, (1.0 / 3.0));
		}
	}
	this->X_Mat_1 = X_Mat_1_temp;

	double** X_Mat_2_temp = new double* [I];
	for (unsigned int i = 0; i < I; i++) {
		X_Mat_2_temp[i] = new double[I];
		for (unsigned int j = 0; j < I; j++) {
			X_Mat_2_temp[i][j] = pow((pow(x[i], 3) - pow(X_Mat_1[i][j], 3)), (1.0 / 3.0));
		}
	}
	this->X_Mat_2 = X_Mat_2_temp;

	double** Co_B_Coef_temp = new double* [I];
	for (unsigned int i = 0; i < I; i++) {
		Co_B_Coef_temp[i] = new double[I];
		for (unsigned int j = 0; j < I - 1; j++) {
			if ((i == 0) || (j == 0))
				Co_B_Coef_temp[i][j] = 0;
			else
			Co_B_Coef_temp[i][j] = t_c / Vm * 2.0 / pow(2.0, (1.0 / 3.0)) * C_r_s[j][i] * W[j] / pow(x[i], 2) / pow(x[j], 3) / pow((1.0 - pow(x[j], 3.0) / 2.0), (5.0 / 3.0));
		}
		Co_B_Coef_temp[i][I - 1] = 0;
	}
	this->Co_B_Coef = Co_B_Coef_temp;

	double** Co_D_Coef_temp = new double* [I];
	for (unsigned int i = 0; i < I; i++) {
		Co_D_Coef_temp[i] = new double[I];
		for (unsigned int j = 0; j < I; j++) {
			if (j==0)
				Co_D_Coef_temp[i][j] = 0;
			else
			Co_D_Coef_temp[i][j] = (t_c / Vm) * W[j] * C_r[j][i] / pow(x[j], 3);
		}
	}
	this->Co_D_Coef = Co_D_Coef_temp;

	double* dy_dx_temp = new double[I];
	double* m_temp = new double[I];
	double* c_coef_temp = new double[I - 1];
	double* d_coef_temp = new double[I - 1];
	double* sai1_temp = new double[I];
	double* sai2_temp = new double[I];

	this->dy_dx = dy_dx_temp;
	this->m = m_temp;
	this->c_coef = c_coef_temp;
	this->d_coef = d_coef_temp;
	this->sai1 = sai1_temp;
	this->sai2 = sai2_temp;
}

template<typename Grid_Type>
void PBE_0D_Solution<Grid_Type>::Breakage_Death_Birth_Coef() {
	double** beta = Daughter_->beta;
	double** beta_s = Daughter_->beta_Shifted;
	double* B_f = Breakage_->B_f;
	double** B_f_s = Breakage_->B_f_Shifted;
	double* x = Breakage_->Grid->Grid_Points_;
	double* W = Breakage_->Grid->Quadrature_weight_;
	double t_c = Breakage_->Grid->t_c;
	unsigned int I = Breakage_->Grid->TotalPointsNumber;

	double** X_Mat_3_temp = new double* [I];
	for (unsigned int i = 0; i < I; i++) {
		X_Mat_3_temp[i] = new double[I];
		for (unsigned int j = 0; j < I; j++) {
			X_Mat_3_temp[i][j] = x[i] * (1-x[j]) + x[j];
		}
	}
	this->X_Mat_3 = X_Mat_3_temp;

	double** Br_B_Coef_temp = new double* [I];
	for (unsigned int i = 0; i < I; i++) {
		Br_B_Coef_temp[i] = new double[I];
		for (unsigned int j = 0; j < I; j++) {
			Br_B_Coef_temp[i][j] = 2.0 * t_c * pow(x[i], 3) * (1- x[i]) / pow((x[j] * (1- x[i])+ x[i]), 3) * B_f_s[i][j] * beta_s[i][j] * W[j];
		}
	}
	Br_B_Coef_temp[0][0] = 0;
	this->Br_B_Coef = Br_B_Coef_temp;

	double* Br_D_Coef_temp = new double[I];
	for (unsigned int i = 0; i < I; i++)
		Br_D_Coef_temp[i] = t_c * B_f[i];
	this->Br_D_Coef = Br_D_Coef_temp;

	double* sai3_temp = new double[I];
	this->sai3 = sai3_temp;

	if (Flag_BC==3){
		double* dy_dx_temp = new double[I];
		double* m_temp = new double[I];
		double* c_coef_temp = new double[I - 1];
		double* d_coef_temp = new double[I - 1];

		this->dy_dx = dy_dx_temp;
		this->m = m_temp;
		this->c_coef = c_coef_temp;
		this->d_coef = d_coef_temp;
	}
}

template<typename Grid_type>
double* PBE_0D_Solution<Grid_type>::PCHIP_Mono_Interpolation(const double* x, const state_type& y, const double* a, const unsigned int& n, const unsigned int& n1, double* yhat, double Inter_Para) {
	unsigned int i, index;
	double alpha, beta, to, xx;
	/*
	double* dy_dx = new double[n];
	double* m = new double[n];
	double* c_coef = new double[n - 1];
	double* d_coef = new double[n - 1];
	double* yhat = new double[n1];
	*/

	dy_dx[0] = (y[1] - y[0]) / (x[1] - x[0]);
	m[0] = dy_dx[0];

	for (i = 1; i < n - 1; i++)
	{
		dy_dx[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
		if ((dy_dx[i - 1] > 0) == (dy_dx[i] < 0))
			//if ((dy_dx[i - 1] * dy_dx[i])<=0)  //slower
			m[i] = 0;
		else
			m[i] = (dy_dx[i - 1] + dy_dx[i]) / 2;

		if (dy_dx[i - 1] == dy_dx[i]) {
			m[i] = 0;
			m[i - 1] = 0;
		}
		else {
			alpha = m[i - 1] / dy_dx[i - 1];
			beta = m[i] / dy_dx[i - 1];
			if (alpha < 0) { m[i - 1] = 0; }
			if (beta < 0) { m[i] = 0; }
			to = Inter_Para / pow(alpha * alpha + beta * beta, 0.5);
			if (to < 1) {
				m[i - 1] = to * alpha * dy_dx[i - 1];
				m[i] = to * beta * dy_dx[i - 1];
			}
		}
	}

	i = n - 1;
	dy_dx[i] = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
	m[i] = dy_dx[i];
	//std::cout << "Hello I am Here" << std::endl;

	if ((dy_dx[i - 1] > 0) == (dy_dx[i] < 0))
		//if ((dy_dx[i - 1] * dy_dx[i])<=0)  //slower
		m[i] = 0;
	else
		m[i] = (dy_dx[i - 1] + dy_dx[i]) / 2;

	if (dy_dx[i - 1] == dy_dx[i]) {
		m[i] = 0;
		m[i - 1] = 0;
	}
	else {
		alpha = m[i - 1] / dy_dx[i - 1];
		beta = m[i] / dy_dx[i - 1];
		if (alpha < 0) { m[i - 1] = 0; }
		if (beta < 0) { m[i] = 0; }
		to = Inter_Para / pow(alpha * alpha + beta * beta, 0.5);
		if (to < 1) {
			m[i - 1] = to * alpha * dy_dx[i - 1];
			m[i] = to * beta * dy_dx[i - 1];
		}
	}

	for (i = 0; i < n - 1; i++)
	{
		c_coef[i] = (3 * dy_dx[i] - 2 * m[i] - m[i + 1]) / (x[i + 1] - x[i]);
		d_coef[i] = (m[i] - 2 * dy_dx[i] + m[i + 1]) / pow((x[i + 1] - x[i]), 2);
	}

	for (i = 0; i < n1; i++) {
		index = Binary_Search(x, a[i], n1);
		xx = a[i] - x[index];
		yhat[i] = y[index] + xx * (m[index] + xx * (c_coef[index] + xx * d_coef[index]));
	}

	/*
	clean_up(dy_dx);
	clean_up(m);
	clean_up(c_coef);
	clean_up(d_coef);
	*/
	return yhat;
}

template<typename Grid_type>
unsigned int PBE_0D_Solution<Grid_type>::Binary_Search(const double* x, const double xhat, const unsigned int n) {
	unsigned int ia = 0, im;
	unsigned int ib = n - 1;
	while (ib - ia > 1) {
		im = (unsigned int)((ia + ib) / 2);
		if (x[im] < xhat)
			ia = im;
		else
			ib = im;
	}
	return ia;
}

template<typename Grid_type>
void PBE_0D_Solution<Grid_type>::ODE_Set_B(const state_type& sai, state_type& dsai_dt, const double& t) {
	double Br_B_Rate;
	double Br_D_Rate;
	double Br_Rate;
	for (unsigned int i = 0; i < Breakage_->Grid->TotalPointsNumber; i++) {
		Br_D_Rate = Br_D_Coef[i] * sai[i];
		Br_B_Rate = 0;
		sai3 = PCHIP_Mono_Interpolation(Breakage_->Grid->Grid_Points_, sai, X_Mat_3[i], Breakage_->Grid->TotalPointsNumber, Breakage_->Grid->TotalPointsNumber, sai3);
		for (unsigned int j = 0; j < Breakage_->Grid->TotalPointsNumber; j++) {
			Br_B_Rate += Br_B_Coef[i][j] * sai3[j];
		}
		dsai_dt[i] = Br_B_Rate - Br_D_Rate;
		//std::cout << Br_B_Rate << "      " << Br_D_Rate << std::endl;
		//std::cout << Br_Rate << std::endl;
	}
}

template<typename Grid_type>
void PBE_0D_Solution<Grid_type>::ODE_Set_C(const state_type& sai, state_type& dsai_dt, const double& t) {
	double Co_B_Rate;
	double Co_D_Rate;
	double Co_Rate;
	for (unsigned int i = 0; i < Coalescence_->Grid->TotalPointsNumber; i++) {
		Co_B_Rate = 0;
		Co_D_Rate = 0;
		sai1 = PCHIP_Mono_Interpolation(Coalescence_->Grid->Grid_Points_, sai, X_Mat_1[i], Coalescence_->Grid->TotalPointsNumber, Coalescence_->Grid->TotalPointsNumber, sai1);
		sai2 = PCHIP_Mono_Interpolation(Coalescence_->Grid->Grid_Points_, sai, X_Mat_2[i], Coalescence_->Grid->TotalPointsNumber, Coalescence_->Grid->TotalPointsNumber, sai2);
		for (unsigned int j = 0; j < Coalescence_->Grid->TotalPointsNumber; j++) {
			Co_B_Rate += Co_B_Coef[i][j] * sai1[j] * sai2[j];
			Co_D_Rate += Co_D_Coef[i][j] * sai[j];
		}
		Co_D_Rate *= sai[i];
		dsai_dt[i] = Co_B_Rate - Co_D_Rate;
		//std::cout << Co_B_Rate << "      " << Co_D_Rate << std::endl;
		//std::cout << Co_Rate << std::endl;
	}
}

template<typename Grid_type>
void PBE_0D_Solution<Grid_type>::ODE_Set_BC(const state_type& sai, state_type& dsai_dt, const double& t) {
	double Br_B_Rate;
	double Br_D_Rate;
	double Br_Rate;
	double Co_B_Rate;
	double Co_D_Rate;
	double Co_Rate;
	//double Rate;

	for (unsigned int i = 0; i < Coalescence_->Grid->TotalPointsNumber; i++) {
		Br_D_Rate = Br_D_Coef[i] * sai[i];
		Br_B_Rate = 0;
		Co_B_Rate = 0;
		Co_D_Rate = 0;
		sai1 = PCHIP_Mono_Interpolation(Coalescence_->Grid->Grid_Points_, sai, X_Mat_1[i], Coalescence_->Grid->TotalPointsNumber, Coalescence_->Grid->TotalPointsNumber, sai1);
		sai2 = PCHIP_Mono_Interpolation(Coalescence_->Grid->Grid_Points_, sai, X_Mat_2[i], Coalescence_->Grid->TotalPointsNumber, Coalescence_->Grid->TotalPointsNumber, sai2);
		sai3 = PCHIP_Mono_Interpolation(Breakage_->Grid->Grid_Points_, sai, X_Mat_3[i], Breakage_->Grid->TotalPointsNumber, Breakage_->Grid->TotalPointsNumber, sai3);
		for (unsigned int j = 0; j < Coalescence_->Grid->TotalPointsNumber; j++) {
			Br_B_Rate += Br_B_Coef[i][j] * sai3[j];
			Co_B_Rate += Co_B_Coef[i][j] * sai1[j] * sai2[j];
			Co_D_Rate += Co_D_Coef[i][j] * sai[j];
		}
		Br_Rate = Br_B_Rate - Br_D_Rate;
		Co_D_Rate *= sai[i];
		Co_Rate = Co_B_Rate - Co_D_Rate;
		dsai_dt[i] = Br_Rate + Co_Rate;
		//std::cout << Co_B_Rate << "      " << Co_D_Rate << std::endl;
		//std::cout << Rate << std::endl;
	}
}

template<typename Grid_type>
void PBE_0D_Solution<Grid_type>::ODE_Set_B_IO(const state_type& sai, state_type& dsai_dt, const double& t) {
	double Br_B_Rate;
	double Br_D_Rate;
	double Br_Rate;
	//double Rate;
	for (unsigned int i = 0; i < Breakage_->Grid->TotalPointsNumber; i++) {
		Br_D_Rate = Br_D_Coef[i] * sai[i];
		Br_B_Rate = 0;
		sai3 = PCHIP_Mono_Interpolation(Breakage_->Grid->Grid_Points_, sai, X_Mat_3[i], Breakage_->Grid->TotalPointsNumber, Breakage_->Grid->TotalPointsNumber, sai3);
		for (unsigned int j = 0; j < Breakage_->Grid->TotalPointsNumber; j++) {
			Br_B_Rate += Br_B_Coef[i][j] * sai3[j];
		}
		Br_Rate = Br_B_Rate - Br_D_Rate;
		dsai_dt[i] = Br_Rate + Time_DL * (Inlet_Dist->dist_sai_RB[i] - sai[i]);
		//std::cout << Br_B_Rate << "      " << Br_D_Rate << std::endl;
		//std::cout << Rate << std::endl;
	}
}

template<typename Grid_type>
void PBE_0D_Solution<Grid_type>::ODE_Set_C_IO(const state_type& sai, state_type& dsai_dt, const double& t) {
	double Co_B_Rate;
	double Co_D_Rate;
	double Co_Rate;
	//double Rate;
	for (unsigned int i = 0; i < Coalescence_->Grid->TotalPointsNumber; i++) {
		Co_B_Rate = 0;
		Co_D_Rate = 0;
		sai1 = PCHIP_Mono_Interpolation(Coalescence_->Grid->Grid_Points_, sai, X_Mat_1[i], Coalescence_->Grid->TotalPointsNumber, Coalescence_->Grid->TotalPointsNumber, sai1);
		sai2 = PCHIP_Mono_Interpolation(Coalescence_->Grid->Grid_Points_, sai, X_Mat_2[i], Coalescence_->Grid->TotalPointsNumber, Coalescence_->Grid->TotalPointsNumber, sai2);
		for (unsigned int j = 0; j < Coalescence_->Grid->TotalPointsNumber; j++) {
			Co_B_Rate += Co_B_Coef[i][j] * sai1[j] * sai2[j];
			Co_D_Rate += Co_D_Coef[i][j] * sai[j];
		}
		Co_D_Rate *= sai[i];
		Co_Rate = Co_B_Rate - Co_D_Rate;
		dsai_dt[i] = Co_Rate + Time_DL * (Inlet_Dist->dist_sai_RB[i] - sai[i]);
		//std::cout << Co_B_Rate << "      " << Co_D_Rate << std::endl;
		//std::cout << Rate << std::endl;
	}
}

template<typename Grid_type>
void PBE_0D_Solution<Grid_type>::ODE_Set_BC_IO(const state_type& sai, state_type& dsai_dt, const double& t) {
	double Br_B_Rate;
	double Br_D_Rate;
	double Br_Rate;
	double Co_B_Rate;
	double Co_D_Rate;
	double Co_Rate;
	//double Rate;
	for (unsigned int i = 0; i < Coalescence_->Grid->TotalPointsNumber; i++) {
		Br_D_Rate = Br_D_Coef[i] * sai[i];
		Br_B_Rate = 0;
		Co_B_Rate = 0;
		Co_D_Rate = 0;
		sai1 = PCHIP_Mono_Interpolation(Coalescence_->Grid->Grid_Points_, sai, X_Mat_1[i], Coalescence_->Grid->TotalPointsNumber, Coalescence_->Grid->TotalPointsNumber, sai1);
		sai2 = PCHIP_Mono_Interpolation(Coalescence_->Grid->Grid_Points_, sai, X_Mat_2[i], Coalescence_->Grid->TotalPointsNumber, Coalescence_->Grid->TotalPointsNumber, sai2);
		sai3 = PCHIP_Mono_Interpolation(Breakage_->Grid->Grid_Points_, sai, X_Mat_3[i], Breakage_->Grid->TotalPointsNumber, Breakage_->Grid->TotalPointsNumber, sai3);
		for (unsigned int j = 0; j < Coalescence_->Grid->TotalPointsNumber; j++) {
			Br_B_Rate += Br_B_Coef[i][j] * sai3[j];
			Co_B_Rate += Co_B_Coef[i][j] * sai1[j] * sai2[j];
			Co_D_Rate += Co_D_Coef[i][j] * sai[j];
		}
		Br_Rate = Br_B_Rate - Br_D_Rate;
		Co_D_Rate *= sai[i];
		Co_Rate = Co_B_Rate - Co_D_Rate;
		dsai_dt[i] = Br_Rate + Co_Rate + Time_DL * (Inlet_Dist->dist_sai_RB[i] - sai[i]);
		//std::cout << Co_B_Rate << "      " << Co_D_Rate << std::endl;
		//std::cout << Rate << std::endl;
	}
}

template<typename Grid_type>
void PBE_0D_Solution<Grid_type>::ODE_Set_NBC_IO(const state_type& sai, state_type& dsai_dt, const double& t) {
	//double Rate;
	for (unsigned int i = 0; i < Grid_NBC->TotalPointsNumber; i++) {
		dsai_dt[i] = Time_DL * (Inlet_Dist->dist_sai_RB[i] - sai[i]);
		//std::cout << Rate << std::endl;
	}
}

template<typename Grid_type>
void PBE_0D_Solution<Grid_type>::Save_Results(const state_type& sai, const double& t) {
	if (((Index % Save_Option) == 0) || (t == Final_Time)) {
		double* sai_temp = new double[Grid_Size];

		for (unsigned int i = 0; i < Grid_Size; i++)
			sai_temp[i] = sai[i];

		Dist_Sol.push_back(sai_temp);
		Time_Sol.push_back(t);
		//cout << t << '\t' << x[0] << '\t' << x[1] << '\t' << x[2] << endl;
	}
	Index += 1;
}

template<typename Grid_type>
void PBE_0D_Solution<Grid_type>::Set_Solver(std::string Integrator_Lib, std::string Integrator_Method, double RelTol, double AbsTol, unsigned int Save_Option, double dt_Initial) {
	this->Integrator_Lib = Integrator_Lib;
	this->Integrator_Method = Integrator_Method;
	this->RelTol = RelTol;
	this->AbsTol = AbsTol;
	this->Save_Option = Save_Option;
	this->dt_Initial = dt_Initial;
}

template<typename Grid_type>
void PBE_0D_Solution<Grid_type>::Simulate(double Final_Time_Length) {
	auto t1 = std::chrono::high_resolution_clock::now();
	if (Flag_Module == 3) {
		this->Final_Time = Final_Time_Length / Velocity;
	}
	else {
		this->Final_Time = Final_Time_Length;
	}

	//state_type sai0;
	state_type sai0(Grid_Size,0.0);   //for ublas vector and arrays

	for (unsigned int i = 0; i < Grid_Size; i++) {
		//sai0.push_back(Int_Dist->dist_sai_RB[i]);
		sai0[i] = Int_Dist->dist_sai_RB[i];  //for ublas vector and arrays
	}
	if (Integrator_Lib == "ODEINT") {
		if (Integrator_Method == "runge_kutta_cash_karp54") {
			typedef runge_kutta_cash_karp54< state_type > stepper_type;
			integrate_adaptive(make_controlled(AbsTol, RelTol, stepper_type()),
				[&](const state_type& sai, state_type& dsai_dt, const double& t) {
					(this->*ODE_Set_fun)(sai, dsai_dt, t); },
				sai0, 0.0, Final_Time, dt_Initial,
						[&](const state_type& sai, const double& t) {
						Save_Results(sai, t); });
			this->Step_Number = Time_Sol.size();
			this->Flag_Solution = true;
		}
		else if (Integrator_Method == "runge_kutta_dopri5") {
			typedef runge_kutta_dopri5< state_type > stepper_type;
			integrate_adaptive(make_controlled(AbsTol, RelTol, stepper_type()),
				[&](const state_type& sai, state_type& dsai_dt, const double& t) {
					(this->*ODE_Set_fun)(sai, dsai_dt, t); },
				sai0, 0.0, Final_Time, dt_Initial,
						[&](const state_type& sai, const double& t) {
						Save_Results(sai, t); });
			this->Step_Number = Time_Sol.size();
			this->Flag_Solution = true;
		}
		else if (Integrator_Method == "runge_kutta_fehlberg78") {
			typedef runge_kutta_fehlberg78< state_type > stepper_type;
			integrate_adaptive(make_controlled(AbsTol, RelTol, stepper_type()),
				[&](const state_type& sai, state_type& dsai_dt, const double& t) {
					(this->*ODE_Set_fun)(sai, dsai_dt, t); },
				sai0, 0.0, Final_Time, dt_Initial,
						[&](const state_type& sai, const double& t) {
						Save_Results(sai, t); });
			this->Step_Number = Time_Sol.size();
			this->Flag_Solution = true;
		}

		else if (Integrator_Method == "STI") {
			typedef runge_kutta_dopri5< state_type > stepper_type;
			//typedef rosenbrock4_dense_output< state_type > stepper_type;

			//integrate_const(make_controlled(AbsTol, RelTol, stepper_type()),
			//integrate_adaptive(make_controlled(AbsTol, RelTol, stepper_type()),
			//integrate_const(make_dense_output(AbsTol, RelTol, stepper_type()),
			integrate_adaptive(make_dense_output(AbsTol, RelTol, stepper_type()),
				[&](const state_type& sai, state_type& dsai_dt, const double& t) {
					(this->*ODE_Set_fun)(sai, dsai_dt, t); },
				sai0, 0.0, Final_Time, dt_Initial,
						[&](const state_type& sai, const double& t) {
						Save_Results(sai, t); });
			this->Step_Number = Time_Sol.size();
			this->Flag_Solution = true;
		}
		else {
			Logger::Instance().Add_Error("PBE_0D_Solution_Function.E2	         Specified Integrator Method is unknown.\n");
			//exit(-1);
			this->Flag_Solution = false;
		}
	}
	else {
		Logger::Instance().Add_Error("PBE_0D_Solution_Function.E1	         Specified Integrator Library is unknown.\n");
		//exit(-1);
		this->Flag_Solution = false;
	}
	if (Flag_Module == 3) {
		for (unsigned int i = 0; i < Step_Number; i++) {
			Time_Sol[i] = Time_Sol[i] * Velocity;
		}
	}
	auto t2 = std::chrono::high_resolution_clock::now();
	this->Execution_Time = (int)(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()/1000);
}