#include "Grid_FE_Declaration.h"
#include "Grid_OC_Declaration.h"

Grid_FE::Grid_FE(const double* ElementBoundaryNodes, const unsigned int* ElementPointsNumber, unsigned int ElementNumber, bool _1st_derivatic_bool, bool _2nd_derivatic_bool, double R_max, double t_c) {
	this->ElementBoundaryNodes=ElementBoundaryNodes;
	this->ElementPointsNumber= ElementPointsNumber;
	this->ElementNumber = ElementNumber;	
	this->R_max = R_max;
	this->t_c = t_c;
	this->Boundary_Delete = false;
	this->_1st_derivatic_bool = _1st_derivatic_bool;
	this->_2nd_derivatic_bool = _2nd_derivatic_bool;

	weights();
}

Grid_FE::Grid_FE(double Int_01, double Mean, double Int_99, const unsigned int* ElementPointsNumber, unsigned int ElementNumber, bool _1st_derivatic_bool, bool _2nd_derivatic_bool, double R_max, double t_c) {

	double* ElementBoundaryNodes_temp = new double[ElementNumber+1];
	ElementBoundaryNodes_temp[0] = 0.0;
	switch (ElementNumber) {
	case 0:
		//warning
		exit;
		break;
	case 1:
		ElementBoundaryNodes_temp[1] = 1;
		break;
	case 2:
		ElementBoundaryNodes_temp[1] = Int_99;
		ElementBoundaryNodes_temp[2] = 1;
		break;
	case 3:
		ElementBoundaryNodes_temp[1] = Mean;
		ElementBoundaryNodes_temp[2] = Int_99;
		ElementBoundaryNodes_temp[3] = 1;
		break;
	case 4:
		ElementBoundaryNodes_temp[1] = Int_01;
		ElementBoundaryNodes_temp[2] = Mean;
		ElementBoundaryNodes_temp[3] = Int_99;
		ElementBoundaryNodes_temp[4] = 1;
		break;
	default:
		ElementBoundaryNodes_temp[1] = Int_01;
		ElementBoundaryNodes_temp[2] = Mean;
		double dx = (log10(1) - log10(Int_99)) / (ElementNumber - 3);
		double x= log10(Int_99);
		for (unsigned int i = 3; i <= ElementNumber; i++) {
			ElementBoundaryNodes_temp[i] = pow(10, x);
			x += dx ;
		}
	}
	this->ElementBoundaryNodes = ElementBoundaryNodes_temp;

	//std::cout << std::endl << std::endl;
	//for (int i = 0; i < ElementNumber+1; i++)
	//	std::cout << ElementBoundaryNodes[i] << "         ";
	//std::cout << std::endl << std::endl;
	this->ElementPointsNumber = ElementPointsNumber;
	this->ElementNumber = ElementNumber;
	this->R_max = R_max;
	this->t_c = t_c;
	this->Boundary_Delete = true;
	this->_1st_derivatic_bool = _1st_derivatic_bool;
	this->_2nd_derivatic_bool = _2nd_derivatic_bool;
	weights();
}


Grid_FE::~Grid_FE() {
	//std::cout << "Grid_FE object destroyed" << std::endl;
	Grid_FE::clean_up(Grid_Points_);

	Grid_FE::clean_up(Quadrature_weight_);

	if (_1st_derivatic_bool)
		Grid_FE::clean_up(_1st_derivatic_weight_);
	else
		delete[] _1st_derivatic_weight_;


	if (_2nd_derivatic_bool)
		Grid_FE::clean_up(_2nd_derivatic_weight_);
	else
		delete[] _2nd_derivatic_weight_;

	if (Boundary_Delete)
		delete[] ElementBoundaryNodes;
}

void Grid_FE::clean_up(double arr[]) {
	delete[] arr;
}

void Grid_FE::clean_up(double** arr) {
	for (unsigned int i = 0; i < TotalPointsNumber; i++) {
		delete[] arr[i];
	}
	delete[] arr;
}

void Grid_FE::weights() {
	int TotalPointsNumber = 1 - ElementNumber;
	for (unsigned int i = 0; i < ElementNumber; i++) {
		TotalPointsNumber += ElementPointsNumber[i];
	}
	this->TotalPointsNumber = TotalPointsNumber;

	double** _1st_derivatic_weight_temp;
	_1st_derivatic_weight_temp = new double* [TotalPointsNumber];
	if (_1st_derivatic_bool) {
		for (unsigned int i = 0; i < TotalPointsNumber; i++) {
			_1st_derivatic_weight_temp[i] = new double[TotalPointsNumber];
			std::fill(_1st_derivatic_weight_temp[i], _1st_derivatic_weight_temp[i] + TotalPointsNumber, 0);
		}
	}

	double** _2nd_derivatic_weight_temp;
	_2nd_derivatic_weight_temp = new double* [TotalPointsNumber];
	if (_2nd_derivatic_bool) {
		for (unsigned int i = 0; i < TotalPointsNumber; i++) {
			_2nd_derivatic_weight_temp[i] = new double[TotalPointsNumber];
			std::fill(_2nd_derivatic_weight_temp[i], _2nd_derivatic_weight_temp[i] + TotalPointsNumber, 0);
		}
	}

	double* Quadrature_weight_temp = new double[TotalPointsNumber];
	double* Grid_Points_temp = new double[TotalPointsNumber];

	double* xx0;
	double* qq0;
	double** AA0;
	double** BB0;

	bool derivatic_bool = _1st_derivatic_bool || _2nd_derivatic_bool;

	unsigned int in = 0;
	unsigned int out = 0;
	for (unsigned int i = 0; i < ElementNumber; i++) {
		if (i == 0) {
			Grid_OC mygrid(ElementPointsNumber[i] - 2, 1, 1);
			xx0 = mygrid.Grid_Points_;
			qq0 = mygrid.Quadrature_weight_;

			//xx0[0] = 1e-60;
			out = ElementPointsNumber[0];
			for (unsigned int j = in; j < out; j++) {
				Grid_Points_temp[j] = xx0[j - in] * (ElementBoundaryNodes[i + 1] - ElementBoundaryNodes[i]) + ElementBoundaryNodes[i];
				Quadrature_weight_temp[j] = qq0[j - in] * (ElementBoundaryNodes[i + 1] - ElementBoundaryNodes[i]);
				//std::cout << xx0[0] << std::endl;
				if (derivatic_bool) {
					for (unsigned int w = in; w < out; w++) {
						if (_1st_derivatic_bool) {
							AA0 = mygrid._1st_derivatic_weight();
							_1st_derivatic_weight_temp[j][w] = AA0[j - in][w - in] / (ElementBoundaryNodes[i + 1] - ElementBoundaryNodes[i]);
						}
						if (_2nd_derivatic_bool) {
							BB0 = mygrid._2nd_derivatic_weight();
							_2nd_derivatic_weight_temp[j][w] = BB0[j - in][w - in] / pow((ElementBoundaryNodes[i + 1] - ElementBoundaryNodes[i]), 2);
						}
					}
				}
			}
			this->Grid_Points_ = Grid_Points_temp;
			this->Quadrature_weight_ = Quadrature_weight_temp;
			this->_1st_derivatic_weight_ = _1st_derivatic_weight_temp;
			this->_2nd_derivatic_weight_ = _2nd_derivatic_weight_temp;
		}

		else {
			Grid_OC mygrid(ElementPointsNumber[i] - 2, 0, 1);
			xx0 = mygrid.Grid_Points_;
			qq0 = mygrid.Quadrature_weight_;
			in = out;
			out += ElementPointsNumber[i] - 1;
			for (unsigned int j = in; j < out; j++) {
				Grid_Points_temp[j] = xx0[j - in] * (ElementBoundaryNodes[i + 1] - ElementBoundaryNodes[i]) + ElementBoundaryNodes[i];
				Quadrature_weight_temp[j] = qq0[j - in] * (ElementBoundaryNodes[i + 1] - ElementBoundaryNodes[i]);
				//std::cout << xx0[0] << std::endl;
				if (derivatic_bool) {
					for (unsigned int w = in; w < out; w++) {
						if (_1st_derivatic_bool) {
							AA0 = mygrid._1st_derivatic_weight();
							_1st_derivatic_weight_temp[j][w] = AA0[j - in][w - in] / (ElementBoundaryNodes[i + 1] - ElementBoundaryNodes[i]);
						}
						if (_2nd_derivatic_bool) {
							BB0 = mygrid._2nd_derivatic_weight();
							_2nd_derivatic_weight_temp[j][w] = BB0[j - in][w - in] / pow((ElementBoundaryNodes[i + 1] - ElementBoundaryNodes[i]), 2);
						}
					}
				}
			}
			this->Grid_Points_ = Grid_Points_temp;
			this->Quadrature_weight_ = Quadrature_weight_temp;
			this->_1st_derivatic_weight_ = _1st_derivatic_weight_temp;
			this->_2nd_derivatic_weight_ = _2nd_derivatic_weight_temp;
		}
	}
}