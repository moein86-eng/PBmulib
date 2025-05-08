#include "Logger_Declaration.h"
#include "Distribution_declaration.h"


#define pi_const 3.141592653589793238462643383279502


template<typename Grid_type>
Distribution<Grid_type>::Distribution(const Grid_type* Grid, double* dist_sai_RB)
	: dist_sai_RB_Bool(0),
	Droplet_Grid_Radius_(nullptr), 
	Droplet_Grid_Volume_(nullptr),
	Droplet_Grid_Volume_DL_(nullptr),
	Density_VBased_DL_(nullptr),
	Droplet_Grid_Diameter_(nullptr),
	Volume_Density_RBased_(nullptr),
	Volume_Density_VBased_(nullptr),
	Number_Density_RBased_(nullptr),
	Number_Density_VBased_(nullptr) {
	this->Grid = Grid;
	this->dist_sai_RB = dist_sai_RB;
	this->Delete_User = false;
}

template<typename Grid_type>
Distribution<Grid_type>::Distribution(const Grid_type* Grid, double mu, double fi, double VolFrac)
	: dist_sai_RB_Bool(1),
	Droplet_Grid_Radius_(nullptr),
	Droplet_Grid_Volume_(nullptr),
	Droplet_Grid_Volume_DL_(nullptr),
	Density_VBased_DL_(nullptr),
	Droplet_Grid_Diameter_(nullptr),
	Volume_Density_RBased_(nullptr),
	Volume_Density_VBased_(nullptr),
	Number_Density_RBased_(nullptr),
	Number_Density_VBased_(nullptr){
	this->Grid = Grid;

	double* dist_sai_RB_temp = new double[Grid->TotalPointsNumber];
	dist_sai_RB_temp[0] = 0;
	for (unsigned int i = 1; i < Grid->TotalPointsNumber; i++)
		dist_sai_RB_temp[i] = VolFrac * exp(-pow((log(Grid->Grid_Points_[i]) - mu), 2)/ (2 * pow(fi, 2))) / ( (Grid->Grid_Points_[i] * fi * pow(2 * pi_const, 0.5)));
	this->dist_sai_RB = dist_sai_RB_temp;
	this->Delete_User = false;
}

template<typename Grid_type>
Distribution<Grid_type>::Distribution(const Grid_type* Grid, const double* User_dist, const double* User_grid, unsigned int user_N)
	: dist_sai_RB_Bool(1),
	Droplet_Grid_Radius_(nullptr),
	Droplet_Grid_Volume_(nullptr),
	Droplet_Grid_Volume_DL_(nullptr),
	Density_VBased_DL_(nullptr),
	Droplet_Grid_Diameter_(nullptr),
	Volume_Density_RBased_(nullptr),
	Volume_Density_VBased_(nullptr),
	Number_Density_RBased_(nullptr),
	Number_Density_VBased_(nullptr) {
	this->Grid = Grid;

	Distribution::InputChecking(User_grid, User_dist, user_N);

	//this->dist_sai_RB = Distribution::Lagrange_Interpolation(User_grid, User_dist, Grid->Grid_Points_, user_N, Grid->TotalPointsNumber);
	this->dist_sai_RB = Distribution::PCHIP_Mono_Interpolation(User_grid, User_dist, Grid->Grid_Points_, user_N, Grid->TotalPointsNumber);
	this->Delete_User = false;
}

template<typename Grid_type>
Distribution<Grid_type>::Distribution(double mu, double fi, double VolFrac, unsigned int user_N)
	: dist_sai_RB_Bool(1),
	Droplet_Grid_Radius_(nullptr),
	Droplet_Grid_Volume_(nullptr),
	Droplet_Grid_Volume_DL_(nullptr),
	Density_VBased_DL_(nullptr),
	Droplet_Grid_Diameter_(nullptr),
	Volume_Density_RBased_(nullptr),
	Volume_Density_VBased_(nullptr),
	Number_Density_RBased_(nullptr),
	Number_Density_VBased_(nullptr) {

	double* User_grid_temp = new double[user_N];
	double* User_dist_temp = new double[user_N];
	double dx = 1.0 / (user_N - 1);
	User_grid_temp[0] = 0;
	User_dist_temp[0] = 0;
		for (unsigned int i = 1; i < user_N; i++){
		User_grid_temp[i] = i*dx;
		User_dist_temp[i] = VolFrac * exp(-pow((log(User_grid_temp[i]) - mu), 2) / (2 * pow(fi, 2))) / ((User_grid_temp[i] * fi * pow(2 * pi_const, 0.5)));
	}

	//this->User_grid = User_grid_temp;
	//this->User_dist = User_dist_temp;
	this->Int_99= Exp_Int_99(User_grid_temp, User_dist_temp, user_N);
	this->Int_01= Exp_Int_01(User_grid_temp, User_dist_temp, user_N);
	this->Mean= Exp_Ave(User_grid_temp, User_dist_temp, user_N);
	Distribution::clean_up(User_grid_temp);
	Distribution::clean_up(User_dist_temp);
	this->Delete_User = false;
	this->user_N = user_N;
}

template<typename Grid_type>
Distribution<Grid_type>::Distribution(double* User_dist, double* User_grid, unsigned int user_N)
	: dist_sai_RB_Bool(1),
	Droplet_Grid_Radius_(nullptr),
	Droplet_Grid_Volume_(nullptr),
	Droplet_Grid_Volume_DL_(nullptr),
	Density_VBased_DL_(nullptr),
	Droplet_Grid_Diameter_(nullptr),
	Volume_Density_RBased_(nullptr),
	Volume_Density_VBased_(nullptr),
	Number_Density_RBased_(nullptr),
	Number_Density_VBased_(nullptr) {

	Distribution::InputChecking(User_grid, User_dist, user_N);

	this->User_grid = User_grid;
	this->User_dist = User_dist;
	this->Int_99 = Exp_Int_99(User_grid, User_dist, user_N);
	this->Int_01 = Exp_Int_01(User_grid, User_dist, user_N);
	this->Mean = Exp_Ave(User_grid, User_dist, user_N);
	this->Delete_User = false;
	this->user_N = user_N;
}


template<typename Grid_type>
Distribution<Grid_type>::~Distribution() {
	if (dist_sai_RB_Bool)
		Distribution::clean_up(dist_sai_RB);

	if (Droplet_Grid_Radius_)
		Distribution::clean_up(Droplet_Grid_Radius_);

	if (Droplet_Grid_Volume_) 
		Distribution::clean_up(Droplet_Grid_Volume_);

	if (Droplet_Grid_Volume_DL_) 
		Distribution::clean_up(Droplet_Grid_Volume_DL_);

	if (Density_VBased_DL_) 
		Distribution::clean_up(Density_VBased_DL_);

	if (Droplet_Grid_Diameter_) 
		Distribution::clean_up(Droplet_Grid_Diameter_);

	if (Volume_Density_RBased_) 
		Distribution::clean_up(Volume_Density_RBased_);

	if (Volume_Density_VBased_) 
		Distribution::clean_up(Volume_Density_VBased_);

	if (Number_Density_RBased_)
		Distribution::clean_up(Number_Density_RBased_);

	if (Number_Density_VBased_)
		Distribution::clean_up(Number_Density_VBased_);

	//if (Delete_User) {
	//	Distribution::clean_up(User_grid);
	//	Distribution::clean_up(User_dist);
	//}

}

template<typename Grid_type>
double* Distribution<Grid_type>::Droplet_Grid_Radius() {
	double* _Droplet_Grid_Radius = new double[Grid->TotalPointsNumber];
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		_Droplet_Grid_Radius[i] = Grid->Grid_Points_[i] * Grid->R_max;
	}
	Droplet_Grid_Radius_ = _Droplet_Grid_Radius;
	return _Droplet_Grid_Radius;
}

template<typename Grid_type>
double* Distribution<Grid_type>::Droplet_Grid_Diameter() {
	double* _Droplet_Grid_Diameter = new double[Grid->TotalPointsNumber];
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		_Droplet_Grid_Diameter[i] = Grid->Grid_Points_[i] * Grid->R_max *2;
	}
	Droplet_Grid_Diameter_ = _Droplet_Grid_Diameter;
	return _Droplet_Grid_Diameter;
}

template<typename Grid_type>
double* Distribution<Grid_type>::Droplet_Grid_Volume() {
	double* _Droplet_Grid_Volume = new double[Grid->TotalPointsNumber];
	double coef = 4 * pi_const / 3;
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		double r = Grid->Grid_Points_[i] * Grid->R_max;
		_Droplet_Grid_Volume[i] = coef * pow(r, 3);
	}
	Droplet_Grid_Volume_ = _Droplet_Grid_Volume;
	return _Droplet_Grid_Volume;
}

template<typename Grid_type>
double* Distribution<Grid_type>::Droplet_Grid_Radius_DL() {     //do not need to fee the memory
	return Grid->Grid_Points_;
}

template<typename Grid_type>
double* Distribution<Grid_type>::Droplet_Grid_Volume_DL() {
	double* _Droplet_Grid_Volume_DL = new double[Grid->TotalPointsNumber];
	double vmax = 4 * pi_const * Grid->R_max * Grid->R_max * Grid->R_max / 3;
	double coef = 4 * pi_const / (3 * vmax);
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		double r = Grid->Grid_Points_[i] * Grid->R_max;
		_Droplet_Grid_Volume_DL[i] = coef * pow (r,3);
	}
	Droplet_Grid_Volume_DL_ = _Droplet_Grid_Volume_DL;
	return _Droplet_Grid_Volume_DL;
}


template<typename Grid_type>
double* Distribution<Grid_type>::Density_RBased_DL() {
	return dist_sai_RB;
}


template<typename Grid_type>
double* Distribution<Grid_type>::Volume_Density_RBased() {
	double* _Volume_Density_RBased = new double[Grid->TotalPointsNumber];
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		_Volume_Density_RBased[i] = dist_sai_RB[i] / Grid->R_max;
	}
	Volume_Density_RBased_ = _Volume_Density_RBased;
	return _Volume_Density_RBased;
}


template<typename Grid_type>
double* Distribution<Grid_type>::Number_Density_RBased() {
	double* _Number_Density_RBased = new double[Grid->TotalPointsNumber];
	_Number_Density_RBased[0] = 0;
	for (unsigned int i = 1; i < Grid->TotalPointsNumber; i++) {
		_Number_Density_RBased[i] = dist_sai_RB[i] / Grid->R_max / (4* pi_const *pow(Grid->R_max * Grid->Grid_Points_[i],3)/3);
	}
	Number_Density_RBased_ = _Number_Density_RBased;
	return _Number_Density_RBased;
}


template<typename Grid_type>
double* Distribution<Grid_type>::Density_VBased_DL() {
	double* _Density_RBased_DL = new double[Grid->TotalPointsNumber];
	_Density_RBased_DL[0] = 0;
	for (unsigned int i = 1; i < Grid->TotalPointsNumber; i++) {
		_Density_RBased_DL[i] = dist_sai_RB[i] / (3 * Grid->Grid_Points_[i] * Grid->Grid_Points_[i]);
	}
	Density_VBased_DL_ = _Density_RBased_DL;
	return _Density_RBased_DL;
}


template<typename Grid_type>
double* Distribution<Grid_type>::Volume_Density_VBased() {
	double* _Volume_Density_VBased = new double[Grid->TotalPointsNumber];
	_Volume_Density_VBased[0] = 0;
	for (unsigned int i = 1; i < Grid->TotalPointsNumber; i++) {
		_Volume_Density_VBased[i] = dist_sai_RB[i] / (3 * Grid->Grid_Points_[i] * Grid->Grid_Points_[i]) / Grid->R_max;
	}
	Volume_Density_VBased_ = _Volume_Density_VBased;
	return _Volume_Density_VBased;
}


template<typename Grid_type>
double* Distribution<Grid_type>::Number_Density_VBased() {
	double* _Number_Density_VBased = new double[Grid->TotalPointsNumber];
	_Number_Density_VBased[0] = 0;
	for (unsigned int i = 1; i < Grid->TotalPointsNumber; i++) {
		_Number_Density_VBased[i] = dist_sai_RB[i] / (3 * Grid->Grid_Points_[i] * Grid->Grid_Points_[i]) / Grid->R_max / (4 * pi_const * pow(Grid->R_max * Grid->Grid_Points_[i], 3) / 3);
	}
	Number_Density_VBased_ = _Number_Density_VBased;
	return _Number_Density_VBased;
}


template<typename Grid_type>
double Distribution<Grid_type>::moment__RBased(short int degree) {
	double _moment__RBased = 0;
	for (unsigned int i = 1; i < Grid->TotalPointsNumber; i++) {
		_moment__RBased += Grid->Quadrature_weight_[i]* pow(Grid->Grid_Points_[i],degree) *dist_sai_RB[i];
	}
	return _moment__RBased;
}


template<typename Grid_type>
double Distribution<Grid_type>::moment__VBased(short int degree) {
	return Distribution::moment__RBased(3*degree);
}


template<typename Grid_type>
double Distribution<Grid_type>::Sauter_Mean_Radius_DL() {
	return (Distribution::moment__RBased(0))/ (Distribution::moment__RBased(-1));
}


template<typename Grid_type>
double Distribution<Grid_type>::Sauter_Mean_Radius() {
	return Distribution::Sauter_Mean_Radius_DL() * Grid->R_max;
}


template<typename Grid_type>
double Distribution<Grid_type>::Total_Number_DL() {
	return (Distribution::moment__RBased(-3));
}


template<typename Grid_type>
double Distribution<Grid_type>::Total_Number() {
	return (Distribution::Total_Number_DL())/ (4 * pi_const * pow(Grid->R_max, 3) / 3);
}


template<typename Grid_type>
double Distribution<Grid_type>::Volume_Fraction() {
	return (Distribution::moment__RBased(0));
}


template<typename Grid_type>
double Distribution<Grid_type>::Average_Radius_RBased_DL() {
	return (Distribution::moment__RBased(1)) / (Distribution::moment__RBased(0));
}


template<typename Grid_type>
double Distribution<Grid_type>::Average_Radius_RBased() {
	return (Distribution::Average_Radius_RBased_DL()) * Grid->R_max;
}


template<typename Grid_type>
double Distribution<Grid_type>::Average_Volume_RBased_DL() {
	return pow(Distribution::Average_Radius_RBased_DL(),3);
}


template<typename Grid_type>
double Distribution<Grid_type>::Average_Volume_RBased() {
	return (4 * pi_const * pow(Distribution::Average_Radius_RBased_DL()* Grid->R_max, 3) / 3);
}


template<typename Grid_type>
double Distribution<Grid_type>::Average_Volume_VBased_DL() {
	return (Distribution::moment__RBased(3)) / (Distribution::moment__RBased(0));
}


template<typename Grid_type>
double Distribution<Grid_type>::Average_Volume_VBased() {
	return (Distribution::Average_Volume_VBased_DL()) * (4 * pi_const * pow(Grid->R_max, 3) / 3);
}


template<typename Grid_type>
double Distribution<Grid_type>::Average_Radius_VBased_DL() {
	return pow(Distribution::Average_Volume_VBased_DL(),1.0/3);
}


template<typename Grid_type>
double Distribution<Grid_type>::Average_Radius_VBased() {
	return Distribution::Average_Radius_VBased_DL() * Grid->R_max;
}


template<typename Grid_type>
double Distribution<Grid_type>::Standard_Deviation_RBased_DL() {
	double integal = 0;
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		integal += Grid->Quadrature_weight_[i] * pow(Grid->Grid_Points_[i] - Distribution::Average_Radius_RBased_DL(),2) * dist_sai_RB[i];
	}
	return pow(integal / Distribution::Volume_Fraction(), 0.5);
}


template<typename Grid_type>
double Distribution<Grid_type>::Standard_Deviation_RBased() {
	return Distribution::Standard_Deviation_RBased_DL()* Grid->R_max;
}


template<typename Grid_type>
double Distribution<Grid_type>::Standard_Deviation_VBased_DL() {
	double integal = 0;
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		integal += Grid->Quadrature_weight_[i] * pow(pow(Grid->Grid_Points_[i],3) - Distribution::Average_Volume_VBased_DL(), 2) * dist_sai_RB[i];
	}
	return pow(integal / Distribution::Volume_Fraction(), 0.5);
}


template<typename Grid_type>
double Distribution<Grid_type>::Standard_Deviation_VBased() {
	return Distribution::Standard_Deviation_VBased_DL()* (4 * pi_const * pow(Grid->R_max, 3) / 3);
}


template<typename Grid_type>
double Distribution<Grid_type>::peak_x() {
	double max_sai = 0;
	double max_x = 0;
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		if (dist_sai_RB[i] > max_sai){
			max_sai = dist_sai_RB[i];
			max_x = Grid->Grid_Points_[i];
		}
	}
	return max_x;
}


template<typename Grid_type>
double Distribution<Grid_type>::tail_x(double criteria) {
	double peak = Distribution::peak_x();
	for (unsigned int i=0; i < Grid->TotalPointsNumber; i++) {
		if (dist_sai_RB[i] > criteria / 100 * peak) {
			if (i != 0)
				return Grid->Grid_Points_[i-1];
			else
				return Grid->Grid_Points_[i];
		}
	}
}


template<typename Grid_type>
double Distribution<Grid_type>::head_x(double criteria) {
	double peak = Distribution::peak_x();
	for (int i= Grid->TotalPointsNumber -1; i >=0; i--) {
		if (dist_sai_RB[i] > criteria / 100 * peak) {
			if (i != Grid->TotalPointsNumber -1)
				return Grid->Grid_Points_[i+1];
			else
				return Grid->Grid_Points_[i];
		}
	}
}


template<typename Grid_type>
void Distribution<Grid_type>::clean_up(double arr[]) {
	delete[] arr;
}


template<typename Grid_type>
double Distribution<Grid_type>::zerotail(const double* x, const double* y, const unsigned int &n) {
	// the function either returns the largest x component with zero value in correspondent y component or -1 if there is no zero value on tail 
	double _zerotail = -1;
	for (int i = 0; i < n; i++) {
		if (y[i] == 0)
			_zerotail = x[i];
		else
			break;
	}
	return _zerotail;
}


template<typename Grid_type>
double Distribution<Grid_type>::zerohead(const double* x, const double* y, const unsigned int& n) {
	// the function either returns the smallest x component with zero value in correspondent y component or 2 if there is no zero value on head
	double _zerohead = 2;
	for (int i = n - 1; i >= 0; i--) {
		if (y[i] == 0)
			_zerohead = x[i];
		else
			break;
	}
	return _zerohead;
}


template<typename Grid_type>
double* Distribution<Grid_type>::Lagrange_Interpolation(const double* x, const double* y, const double* a, const unsigned int &n, const unsigned int &n1) {
	/*
	x: pointer to the independent variable array (as raw data)
	y: pointer to the dependent variable array (as raw data)
	a: pointer to interpolation points as independent variable
	n: size of 'x' and 'y' array
	n1: size of 'a' array
	the function returns the interpolated values for dependent variable correspondent to 'a' array
	*/
	double  s, t, k;
	int i, j;
	double* b = new double[n1];
	double _zerotail = zerotail(x, y, n);
	double _zerohead = zerohead(x, y, n);

	for (int w = 0; w < n1; w++) {
		if (a[w]< _zerotail || a[w] > _zerohead)  //it is added to avoid osilation
			k = 0;
		else {
			s = 1; t = 1; k = 0;
			for (i = 0; i < n; i++)
			{
				s = 1;
				t = 1;
				for (j = 0; j < n; j++)
				{
					if (j != i)
					{
						s = s * (a[w] - x[j]);
						t = t * (x[i] - x[j]);
					}
				}
				k = k + ((s / t) * y[i]);
			}
			if (k < 0)
				k = 0;     //added since probability density can not be negative.
		}
		b[w] = k;
	}
	return b;
}


template<typename Grid_type>
double* Distribution<Grid_type>::PCHIP_Mono_Interpolation(const double* x, const double* y, const double* a, const unsigned int& n, const unsigned int& n1, double Inter_Para) {
	
	unsigned int i, index;
	double alpha, beta, to;
	double* dy_dx = new double[n];
	double* m = new double[n];
	double* c_coef = new double[n - 1];
	double* d_coef = new double[n - 1];
	double* yhat = new double[n1];


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
	dy_dx[i] = (y[n-1] - y[n - 2]) / (x[n-1] - x[n - 2]);
	m[i] = dy_dx[i];
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
		index = Binary_Search(x, a[i], n);
		double xx = a[i] - x[index];
		yhat[i] = y[index] + xx * (m[index] + xx * (c_coef[index] + xx * d_coef[index]));
	}

	Distribution::clean_up(dy_dx);
	Distribution::clean_up(m);
	Distribution::clean_up(c_coef);
	Distribution::clean_up(d_coef);

	return yhat;
}

template<typename Grid_type>
unsigned int Distribution<Grid_type>::Binary_Search(const double* x, const double xhat, const unsigned int n) {
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
void Distribution<Grid_type>::InputChecking(const double* User_grid, const double* User_dist, const unsigned int &user_N) {
	/*
	if User_grid is outside 0-1:									returns 1   for   (error1)
	if User_dist is negative   :									returns 1   for   (error2)
	if User_grid is not sorted or duplicated points :				returns 1   for   (warning1)
	if User_dist has more than 1 peak:								returns 1   for   (warning2)
	*/

	short int sum = 0;
	for (int i = 0; i < user_N; i++) {
		if ((User_grid[i] >= User_grid[i + 1]) && i != user_N - 1)
			Logger::Instance().Add_Warning("Distribution.W1				User-specified grid is not sorted or there are duplicated points.\n");//log_W1 = 1;

		if (User_grid[i] < 0 || User_grid[i]>1){
			Logger::Instance().Add_Error("Distribution.E1				User-specified grid is outside 0 - 1.\n");//log_E1 = 1;
			//exit(-1);
		}

		if (User_dist[i] < 0){
			Logger::Instance().Add_Error("Distribution.E2				User-specified distribution is negative.\n");//log_E2 = 1;
			//exit(-2);
		}

		if ((User_dist[i] > User_dist[i + 1]) && (User_dist[i] > User_dist[i - 1]) && i != 0 && i != user_N - 1)
			sum++;
	}
	if (sum > 1)
		Logger::Instance().Add_Warning("Distribution.W2				User-specified distribution has more than 1 peak.\n");//log_W2 = 1;
}

template<typename Grid_type>
void Distribution<Grid_type>::Normalize(double VF0) {
	double coef = VF0/Volume_Fraction();
	for (unsigned int i = 0; i < Grid->TotalPointsNumber; i++) {
		dist_sai_RB[i] = dist_sai_RB[i] * coef;
	}
}

template<typename Grid_type>
double Distribution<Grid_type>::Exp_Trapz(const double* x, const double* y, unsigned int N) {
	double I = 0.0;
	for (unsigned int i = 1; i < N; i++) {
		I += (x[i]- x[i-1])*(y[i] + y[i - 1])/2;
	}
	return I;
}

template<typename Grid_type>
double Distribution<Grid_type>::Exp_Int_99(const double* x, const double* y, unsigned int N) {
	unsigned int i = 1;
	double INT = Exp_Trapz(x, y, N);
	double I99 = 0;
	while (I99 < INT*0.99){
		I99 += (x[i] - x[i - 1]) * (y[i] + y[i - 1]) / 2;
		i += 1;
	}
	if (i == N)
		Logger::Instance().Add_Warning("Distribution.W3				More points on right of the distribution head is recommended for the user specified data.\n");
	return x[i];
}

template<typename Grid_type>
double Distribution<Grid_type>::Exp_Int_01(const double* x, const double* y, unsigned int N) {
	unsigned int i;
	double INT = Exp_Trapz(x, y, N);
	double I01 = 0;
	for (i = 1; i < N; i++) {
		I01 += (x[i] - x[i - 1]) * (y[i] + y[i - 1]) / 2;
		if (I01 > INT * 0.01)
			break;
	}
	if ((i - 1) == N)
		Logger::Instance().Add_Warning("Distribution.W4				More points on left of the distribution head is recommended for the user specified data.\n");

	return x[i-1];
}

template<typename Grid_type>
double Distribution<Grid_type>::Exp_Ave(const double* x, const double* y, unsigned int N) {
	double I_0 = Exp_Trapz(x, y, N);
	double I_1 = 0.0;
	for (unsigned int i = 1; i < N; i++) {
		I_1 += (x[i] - x[i - 1]) * (y[i]* x[i] + y[i - 1]* x[i - 1]) / 2;
	}
	return I_1/ I_0;
}

template<typename Grid_type>
void Distribution<Grid_type>::Set_Grid(const Grid_type* Grid) {
	this->Grid = Grid;
	//this->dist_sai_RB = Distribution::Lagrange_Interpolation(User_grid, User_dist, Grid->Grid_Points_, user_N, Grid->TotalPointsNumber);
	this->dist_sai_RB = Distribution::PCHIP_Mono_Interpolation(User_grid, User_dist, Grid->Grid_Points_, user_N, Grid->TotalPointsNumber);
}