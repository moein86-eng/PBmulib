#include "Grid_OC_Declaration.h"

Grid_OC::Grid_OC(unsigned int colloc_points, unsigned int left, unsigned int right, double R_max, double t_c)
	: _1st_deriv_bool(false), _2nd_deriv_bool(false) {
	this->colloc_points= colloc_points;
	this->left= left;
	this->right= right;
	this->TotalPointsNumber = colloc_points + left + right;
	this->Tol = 1.0e-9;
	this->Grid_Points_ = Grid_OC::jcobi_roots();
	this->Quadrature_weight_ = Grid_OC::Quadrature_weight();
	this->R_max = R_max;
	this->t_c = t_c;
}

Grid_OC::~Grid_OC() {
	//std::cout << "Grid_OC object destroyed" << endl;
	Grid_OC::clean_up(Grid_Points_);

	Grid_OC::clean_up(Quadrature_weight_);

	if(_1st_deriv_bool)
		Grid_OC::clean_up(_1st_derivatic_weight_);

	if(_2nd_deriv_bool)
		Grid_OC::clean_up(_2nd_derivatic_weight_);
}

double* Grid_OC::dif() {
	double* dif1 = new double[TotalPointsNumber];
	double* dif2 = new double[TotalPointsNumber];
	double* dif3 = new double[TotalPointsNumber];
	double* dif_T = new double[TotalPointsNumber *3];
	for (unsigned int i = 0; i < TotalPointsNumber; i++) {
		dif1[i] = 1.0;
		dif2[i] = 0.0;
		dif3[i] = 0.0;
		for (unsigned int j = 0; j < TotalPointsNumber; j++) {
			if (j!=i) {
				double y = Grid_Points_[i] - Grid_Points_[j];
				dif3[i] = y * dif3[i] + 3.0 * dif2[i];
				dif2[i] = y * dif2[i] + 2.0 * dif1[i];
				dif1[i] = y * dif1[i];
			}
		}
		dif_T[i] = dif1[i];
		dif_T[i+ TotalPointsNumber] = dif2[i];
		dif_T[i+ TotalPointsNumber*2] = dif3[i];
	}

	delete[] dif3, dif2, dif1;

	return dif_T;
}


double* Grid_OC::jcobi_roots(double alpha, double beta) {
	double* dif1 = new double[TotalPointsNumber];
	double* dif2 = new double[TotalPointsNumber];
	double* Grid_Points_ = new double[TotalPointsNumber];
	double ab = alpha + beta;
	double ad = beta - alpha;
	double ap = beta * alpha;
	dif1[0] = (ad / (ab + 2.0) + 1.0) / 2.0;
	dif2[0] = 0;
	if (colloc_points >= 2) {
		for (unsigned int i = 1; i < colloc_points; i++) {
			unsigned int z1 = i ;
			double z = ab + 2.0 * z1;
			dif1[i]= (ab * ad / z / (z + 2.0) + 1.0) / 2.0;
			if (i == 1) {
				dif2[i] = (ab + ap + z1) / z / z / (z + 1.0);
			}
			else{
				z = z * z;
				double y = z1 * (ab + z1);
				y= y * (ap + y);
				dif2[i]=y/z/(z-1.0);
			}
		}
	}

	// ROOT DETERMINATION BY NEWTON METHOD WITH SUPPRESSION OF PREVIOUSLY DETERMINED ROOTS
	double x = 0;
	Grid_Points_[0] = 1e66;
	for (unsigned int i = left; i < colloc_points+left; i++) {
		double z = 1;
		while (abs(z) > this->Tol) {
			double xd = 0.0;
			double xn = 1.0;
			double xd1 = 0.0;
			double xn1 = 0.0;
			for (unsigned int j = 0; j < colloc_points; j++) {
				double xp = (dif1[j] - x) * xn - dif2[j] * xd;
				double xp1 = (dif1[j] - x) * xn1 - dif2[j] * xd1 - xn;
				xd = xn;
				xd1 = xn1;
				xn = xp;
				xn1 = xp1;
			}
			double zc = 1.0;
			z = xn / xn1;
			if (i != 0) {
				for (unsigned int j = 1; j < i+1;j++) {
					zc = zc-z/(x-Grid_Points_[j-1]);
				}
			}
			z = z / zc;
			x = x - z;
		}
		Grid_Points_[i] = x;
		x = x + 0.0001;
	}
	if (left == 1)
		Grid_Points_[0] = 0;
	if (right == 1)
		Grid_Points_[TotalPointsNumber - 1] = 1.0;
	delete[] dif2, dif1;
	return Grid_Points_;
}


double* Grid_OC::dfopr(const unsigned int &i, const unsigned int &id, double dif[]) {
	//id shall be 1, 2 or 3

	double* vect = new double[TotalPointsNumber];
	if (id != 3) {
		for (unsigned int j = 0; j < TotalPointsNumber; j++) {
			if (j == i) {
				if (id == 1)
					vect[i] = dif[i+ TotalPointsNumber] / dif[i] / 2.0;
				else 
					vect[i] = dif[i+2* TotalPointsNumber] / dif[i] / 3.0;
			}
			else {
				double y = Grid_Points_[i] - Grid_Points_[j];
				vect[j] = dif[i] / dif[j] / y;
				if(id==2)
					vect[j] = vect[j] * (dif[i + TotalPointsNumber] / dif[i] - 2.0 / y);
			}
		}
	}
	else {
		double y = 0;
		for (unsigned int j = 0; j < TotalPointsNumber; j++) {
			double x = Grid_Points_[j];
			double ax = x * (1.0 - x);
			if (left == 0)
				ax = ax / x / x;
			if (right == 0)
				ax = ax / (1.0 - x) / (1.0 - x);
			vect[j] = ax / dif[j] / dif[j];
			y = y + vect[j];
		}
		for (unsigned int j = 0; j < TotalPointsNumber; j++) {
			vect[j] = vect[j] / y;
		}
	}
	return vect;
}


double* Grid_OC::Quadrature_weight() {
	double* dif_ = Grid_OC::dif();
	double* q = Grid_OC::dfopr(-1, 3, dif_);
	Grid_OC::clean_up(dif_);
	return q;
}


double** Grid_OC::_1st_derivatic_weight() {
	double* dif_ = Grid_OC::dif();
	double** A = 0;
	A = new double* [TotalPointsNumber];

	for (unsigned int i = 0; i < TotalPointsNumber; i++) {
		A[i] = new double[TotalPointsNumber];
		A[i] = Grid_OC::dfopr(i, 1, dif_);
	}
	Grid_OC::clean_up(dif_);
	_1st_deriv_bool = true;
	_1st_derivatic_weight_ = A;
	return A;
}


double** Grid_OC::_2nd_derivatic_weight() {
	double* dif_ = Grid_OC::dif();
	double** B = 0;
	B = new double* [TotalPointsNumber];

	for (unsigned int i = 0; i < TotalPointsNumber; i++) {
		B[i] = new double[TotalPointsNumber];
		B[i] = Grid_OC::dfopr(i, 2, dif_);
	}
	Grid_OC::clean_up(dif_);
	_2nd_deriv_bool = true;
	_2nd_derivatic_weight_ = B;
	return B;
}


void Grid_OC::clean_up(double arr[]) {
	delete[] arr;
}


void Grid_OC::clean_up(double **arr) {
	for (unsigned int i = 0; i < TotalPointsNumber; i++) {
		delete[] arr[i];
	}
	delete[] arr;
}