#pragma warning(disable: 4996)
//#pragma warning(suppress : 4996)

#include "pch.h"


/*typedef std::vector< double > state_type;
const double gam = 0.15;
void harmonic_oscillator(const state_type& x, state_type& dxdt, const double t) {
	dxdt[0] = x[1];
	dxdt[1] = -x[0] - gam * x[1];
}

class harm_osc {

	double m_gam;

public:
	harm_osc(double gam) : m_gam(gam) { }

	void operator() (const state_type& x, state_type& dxdt, const double t)
	{
		dxdt[0] = x[1];
		dxdt[1] = -x[0] - m_gam * x[1];
	}
};*/


using namespace std;


int main()
{
	//cout << std::fixed;
	{
#if 0
    // Grid_OC///////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << "Grid_OC class---------------------------------------------------" << endl;
	Grid_OC mygrid(4,1,1);
	//Grid_OC mygrid(4, 1, 1, 200e-6);
	//Grid_OC mygrid(4);

	cout << "roots: "<<endl;
    double* ptr_root = mygrid.Grid_Points_;
    for (int i = 0; i < mygrid.TotalPointsNumber; i++) {
        cout << ptr_root[i] << std::setw(14);
    }
    cout << "\n\n";

	
	cout << "quadrature weight: " << endl;
    double* ptr_q = mygrid.Quadrature_weight_;
    for (int i = 0; i < 6; i++) {
        cout << ptr_q[i] << std::setw(14);
    }
    cout << "\n\n";

	
	cout << "1st derivative weight: " << endl;
    double** ptr_A = mygrid._1st_derivatic_weight();
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            cout << ptr_A[i][j] << std::setw(14);
        }
        cout << "\n";
    }
    cout << "\n\n";
	

	cout << "2nd derivative weight: " << endl;
    double** ptr_B = mygrid._2nd_derivatic_weight();
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            cout << ptr_B[i][j] << std::setw(14);
        }
        cout << "\n";
    }
    cout << "\n\n";

	//double** ptr_A2 = mygrid._1st_derivatic_weight_;
	//double** ptr_B2 = mygrid._2nd_derivatic_weight_;

#endif
	
	
#if 0
    // Distribution///////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout << "Distribution class---------------------------------------------------" << endl;


    Grid_OC mygrid1(18,1,1,100e-6);
    double* ptr_root1 = mygrid1.Grid_Points_;
    double* ptr_q1 = mygrid1.Quadrature_weight_;
    double ptr_sai[20] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, .4466, 1.9116, 0.0797, 0.0003, 0, 0, 0, 0, 0, 0 };
	
	/*
	ptr_root1[0] = -0.0001;
	ptr_root1[1] = 0.0001;
	ptr_root1[2] = 0.0001;
	ptr_sai[0] = -0.001;
	ptr_sai[3] = 0.1;
	*/

	Grid_OC mygrid2(38, 1, 1, 100e-6);
	double* ptr_root2 = mygrid2.Grid_Points_;
	double* ptr_q2 = mygrid2.Quadrature_weight_;

	double ElementBoundaryNodes0[] = { 0, 0.3, 0.6, 1 };
	unsigned int ElementPointsNumber0[] = { 4, 20, 15 };
	unsigned int ElementNumber0 = 3;
	Grid_FE myFEGrid0(ElementBoundaryNodes0, ElementPointsNumber0, ElementNumber0, true, true, 100e-6);

	//Distribution<Grid_OC> MyDistribution(&mygrid1, ptr_sai);
	Distribution<Grid_OC> MyDistribution(&mygrid1, -0.5, 0.06, 0.2);
	//Distribution<Grid_OC> MyDistribution(&mygrid2, ptr_sai, ptr_root1, 20);
	//Distribution<Grid_FE> MyDistribution(&myFEGrid0, ptr_sai, ptr_root1, 20);
	
	//MyDistribution.Normalize(0.2);
	//std::cout << std::endl << std::endl << MyDistribution.Exp_Trapz(ptr_root1, ptr_sai, 20) << std::endl << std::endl;
	//std::cout << std::endl << std::endl << MyDistribution.Exp_Int_99(ptr_root1, ptr_sai, 20) << std::endl << std::endl;
	//std::cout << std::endl << std::endl << MyDistribution.Exp_Int_01(ptr_root1, ptr_sai, 20) << std::endl << std::endl;
	//std::cout << std::endl << std::endl << MyDistribution.Exp_Ave(ptr_root1, ptr_sai, 20) << std::endl << std::endl;

    double* rr = MyDistribution.Droplet_Grid_Radius();
    double* dd = MyDistribution.Droplet_Grid_Diameter();
    double* vv = MyDistribution.Droplet_Grid_Volume();
	double* rr_bar = MyDistribution.Droplet_Grid_Radius_DL();
    double* vv_bar = MyDistribution.Droplet_Grid_Volume_DL();
	double* sai_R = MyDistribution.Density_RBased_DL();
    double* fv_R = MyDistribution.Volume_Density_RBased();
    double* fn_R = MyDistribution.Number_Density_RBased();
	double* sai_v = MyDistribution.Density_VBased_DL();
	double* fv_V = MyDistribution.Volume_Density_VBased();
	double* fn_V = MyDistribution.Number_Density_VBased();

	/*
	double* rr2 = MyDistribution.Droplet_Grid_Radius_;
	double* dd2 = MyDistribution.Droplet_Grid_Diameter_;
	double* vv2 = MyDistribution.Droplet_Grid_Volume_;
	double* rr_bar2 = MyDistribution.R_grid;
	double* vv_bar2 = MyDistribution.Droplet_Grid_Volume_DL_;
	double* sai_R2 = MyDistribution.dist_sai_RB;
	double* fv_R2 = MyDistribution.Volume_Density_RBased_;
	double* fn_R2 = MyDistribution.Number_Density_RBased_;
	double* sai_v2 = MyDistribution.Density_VBased_DL_;
	double* fv_V2 = MyDistribution.Volume_Density_VBased_;
	double* fn_V2 = MyDistribution.Number_Density_VBased_;
    */

	std::cout << "ARR:   " << endl;
    for(int i=0; i< 40;i++)
	    std::cout << fn_V[i]<< std::setw(14);
	std::cout << endl << endl;

	std::cout << "moment__RBased(3):              " << MyDistribution.moment__RBased(3)<<endl;
	std::cout << "moment__VBased(-1):             " << MyDistribution.moment__VBased(-1)<<endl;
	std::cout << "Sauter_Mean_Radius_DL:          " << MyDistribution.Sauter_Mean_Radius_DL() << endl;
	std::cout << "Sauter_Mean_Radius:             " << MyDistribution.Sauter_Mean_Radius() << endl;
	std::cout << "Total_Number_DL:                " << MyDistribution.Total_Number_DL() << endl;
	std::cout << "Total_Number:                   " << MyDistribution.Total_Number() << endl;
	std::cout << "Volume_Fraction:                " << MyDistribution.Volume_Fraction() << endl;
	std::cout << "Average_Radius_RBased:          " << MyDistribution.Average_Radius_RBased() << endl;
	std::cout << "Average_Radius_RBased_DL:       " << MyDistribution.Average_Radius_RBased_DL() << endl;
	std::cout << "Average_Radius_RBased:          " << MyDistribution.Average_Radius_RBased() << endl;
	std::cout << "Average_Volume_RBased_DL:       " << MyDistribution.Average_Volume_RBased_DL() << endl;
	std::cout << "Average_Volume_RBased:          " << MyDistribution.Average_Volume_RBased() << endl;
	std::cout << "Average_Volume_VBased_DL:       " << MyDistribution.Average_Volume_VBased_DL() << endl;
	std::cout << "Average_Volume_VBased:          " << MyDistribution.Average_Volume_VBased() << endl;
	std::cout << "Average_Radius_VBased_DL:       " << MyDistribution.Average_Radius_VBased_DL() << endl;
	std::cout << "Average_Radius_VBased:          " << MyDistribution.Average_Radius_VBased() << endl;
	std::cout << "Standard_Deviation_RBased_DL:   " << MyDistribution.Standard_Deviation_RBased_DL() << endl;
	std::cout << "Standard_Deviation_RBased:      " << MyDistribution.Standard_Deviation_RBased() << endl;
	std::cout << "Standard_Deviation_VBased_DL:   " << MyDistribution.Standard_Deviation_VBased_DL() << endl;
	std::cout << "Standard_Deviation_VBased:      " << MyDistribution.Standard_Deviation_VBased() << endl;
	std::cout << "peak_x:                         " << MyDistribution.peak_x() << endl;
	std::cout << "head_x:                         " << MyDistribution.head_x() << endl;
	std::cout << "tail_x:                         " << MyDistribution.tail_x() << endl;
	std::cout << endl << endl;


#endif


#if 0
	// Checking automatic grid///////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout << "Checking automatic grid---------------------------------------------------" << endl;
	
	Grid_OC mygrid1(18, 1, 1, 100e-6);
	double* ptr_root1 = mygrid1.Grid_Points_;
	double ptr_sai[20] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, .4466, 1.9116, 0.0797, 0.0003, 0, 0, 0, 0, 0, 0 };
	
	//Distribution<Grid_FE> MyDistribution(-0.5, 0.06, 0.2,100);
	Distribution<Grid_FE> MyDistribution(ptr_sai, ptr_root1, 20);
	
	std::cout << "I99:                   " << MyDistribution.Int_99 << endl;
	std::cout << "I01:                   " << MyDistribution.Int_01 << endl;
	std::cout << "mean:                  " << MyDistribution.Mean << endl;
	
	std::cout << "user_grid:   " << endl;
	for (int i = 0; i < 20; i++)
		std::cout << MyDistribution.User_grid[i] << std::setw(14);
	std::cout << endl << endl;
	
	for (int i = 0; i < 20; i++)
		std::cout << MyDistribution.User_dist[i] << std::setw(14);
	std::cout << endl << endl;
	
	double ElementBoundaryNodes0[] = { 0, 1 };
	unsigned int ElementPointsNumber0[] = { 15 };
	unsigned int ElementNumber0 = 1;
	Grid_FE myFEGrid0(ElementBoundaryNodes0, ElementPointsNumber0, ElementNumber0, true, true, 100e-6);

	//unsigned int ElementPointsNumber0[] = { 20,10,5,10};
	//unsigned int ElementNumber0 = 4;
	//Grid_FE myFEGrid0(MyDistribution.Int_01, MyDistribution.Mean, MyDistribution.Int_99, ElementPointsNumber0, ElementNumber0, true, true, 100e-6);
	
	std::cout << "Grid points:   " << endl;
	for (int i = 0; i < myFEGrid0.TotalPointsNumber; i++)
		std::cout << myFEGrid0.Grid_Points_[i] << std::setw(14);
	std::cout << endl;
	std::cout << endl << endl;

	std::cout << "Quadrature weight:   " << endl;
	for (int i = 0; i < myFEGrid0.TotalPointsNumber; i++)
		std::cout << myFEGrid0.Quadrature_weight_[i] << std::setw(14);
	std::cout << endl;
	std::cout << endl << endl;

	std::cout << "1st derivative weight: " << endl;
	double** ptr_A_FE = myFEGrid0._1st_derivatic_weight_;
	int nn = myFEGrid0.TotalPointsNumber;
	for (int i = 0; i < nn; i++) {
		for (int j = 0; j < nn; j++) {
			std::cout << ptr_A_FE[i][j] << std::setw(14);
		}
		std::cout << "\n";
	}
	std::cout << "\n\n";


	std::cout << "2nd derivative weight: " << endl;
	double** ptr_B_FE = myFEGrid0._2nd_derivatic_weight_;
	int mm = myFEGrid0.TotalPointsNumber;
	for (int i = 0; i < mm; i++) {
		for (int j = 0; j < nn; j++) {
			std::cout << ptr_B_FE[i][j] << std::setw(14);
		}
		std::cout << "\n";
	}
	std::cout << "\n\n";

	MyDistribution.Set_Grid(&myFEGrid0);
	
	double* rr = MyDistribution.Droplet_Grid_Radius();
	double* dd = MyDistribution.Droplet_Grid_Diameter();
	double* vv = MyDistribution.Droplet_Grid_Volume();
	double* rr_bar = MyDistribution.Droplet_Grid_Radius_DL();
	double* vv_bar = MyDistribution.Droplet_Grid_Volume_DL();
	double* sai_R = MyDistribution.Density_RBased_DL();
	double* fv_R = MyDistribution.Volume_Density_RBased();
	double* fn_R = MyDistribution.Number_Density_RBased();
	double* sai_v = MyDistribution.Density_VBased_DL();
	double* fv_V = MyDistribution.Volume_Density_VBased();
	double* fn_V = MyDistribution.Number_Density_VBased();
	/*
	double* rr2 = MyDistribution.Droplet_Grid_Radius_;
	double* dd2 = MyDistribution.Droplet_Grid_Diameter_;
	double* vv2 = MyDistribution.Droplet_Grid_Volume_;
	double* rr_bar2 = MyDistribution.R_grid;
	double* vv_bar2 = MyDistribution.Droplet_Grid_Volume_DL_;
	double* sai_R2 = MyDistribution.dist_sai_RB;
	double* fv_R2 = MyDistribution.Volume_Density_RBased_;
	double* fn_R2 = MyDistribution.Number_Density_RBased_;
	double* sai_v2 = MyDistribution.Density_VBased_DL_;
	double* fv_V2 = MyDistribution.Volume_Density_VBased_;
	double* fn_V2 = MyDistribution.Number_Density_VBased_;
	*/

	std::cout << "ARR:   " << endl;
	for (int i = 0; i < 20; i++)
		std::cout << sai_R[i] << std::setw(14);
	std::cout << endl << endl;

	std::cout << "moment__RBased(3):              " << MyDistribution.moment__RBased(3) << endl;
	std::cout << "moment__VBased(-1):             " << MyDistribution.moment__VBased(-1) << endl;
	std::cout << "Sauter_Mean_Radius_DL:          " << MyDistribution.Sauter_Mean_Radius_DL() << endl;
	std::cout << "Sauter_Mean_Radius:             " << MyDistribution.Sauter_Mean_Radius() << endl;
	std::cout << "Total_Number_DL:                " << MyDistribution.Total_Number_DL() << endl;
	std::cout << "Total_Number:                   " << MyDistribution.Total_Number() << endl;
	std::cout << "Volume_Fraction:                " << MyDistribution.Volume_Fraction() << endl;
	std::cout << "Average_Radius_RBased:          " << MyDistribution.Average_Radius_RBased() << endl;
	std::cout << "Average_Radius_RBased_DL:       " << MyDistribution.Average_Radius_RBased_DL() << endl;
	std::cout << "Average_Radius_RBased:          " << MyDistribution.Average_Radius_RBased() << endl;
	std::cout << "Average_Volume_RBased_DL:       " << MyDistribution.Average_Volume_RBased_DL() << endl;
	std::cout << "Average_Volume_RBased:          " << MyDistribution.Average_Volume_RBased() << endl;
	std::cout << "Average_Volume_VBased_DL:       " << MyDistribution.Average_Volume_VBased_DL() << endl;
	std::cout << "Average_Volume_VBased:          " << MyDistribution.Average_Volume_VBased() << endl;
	std::cout << "Average_Radius_VBased_DL:       " << MyDistribution.Average_Radius_VBased_DL() << endl;
	std::cout << "Average_Radius_VBased:          " << MyDistribution.Average_Radius_VBased() << endl;
	std::cout << "Standard_Deviation_RBased_DL:   " << MyDistribution.Standard_Deviation_RBased_DL() << endl;
	std::cout << "Standard_Deviation_RBased:      " << MyDistribution.Standard_Deviation_RBased() << endl;
	std::cout << "Standard_Deviation_VBased_DL:   " << MyDistribution.Standard_Deviation_VBased_DL() << endl;
	std::cout << "Standard_Deviation_VBased:      " << MyDistribution.Standard_Deviation_VBased() << endl;
	std::cout << "peak_x:                         " << MyDistribution.peak_x() << endl;
	std::cout << "head_x:                         " << MyDistribution.head_x() << endl;
	std::cout << "tail_x:                         " << MyDistribution.tail_x() << endl;
	std::cout << endl << endl;
	

#endif


#if 0
	// Grid_FE///////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout << "Grid_FE class---------------------------------------------------" << endl;

	double ElementBoundaryNodes[] = {0, 0.3, 0.6, 1};
	unsigned int ElementPointsNumber[] = {3, 4, 4};
	unsigned int ElementNumber = 3;
	//Grid_FE myFEGrid(ElementBoundaryNodes, ElementPointsNumber, ElementNumber, true, true);
	Grid_FE myFEGrid(ElementBoundaryNodes, ElementPointsNumber, ElementNumber, true, true, 200e-6);
	//Grid_FE myFEGrid(ElementBoundaryNodes, ElementPointsNumber, ElementNumber);

	std::cout << "Grid points:   " << endl;
	for (int i = 0; i < myFEGrid.TotalPointsNumber;i++)
		std::cout << myFEGrid.Grid_Points_[i] << std::setw(14);
	std::cout << endl;
	std::cout << endl << endl;

	std::cout << "Quadrature weight:   " << endl;
	for (int i = 0; i < myFEGrid.TotalPointsNumber;i++)
		std::cout << myFEGrid.Quadrature_weight_[i] << std::setw(14);
	std::cout << endl;
	std::cout << endl << endl;

	
	std::cout << "1st derivative weight: " << endl;
	double** ptr_A_FE = myFEGrid._1st_derivatic_weight_;
	int nn = myFEGrid.TotalPointsNumber;
	for (int i = 0; i < nn; i++) {
		for (int j = 0; j <nn; j++) {
			std::cout << ptr_A_FE[i][j] << std::setw(14);
		}
		std::cout << "\n";
	}
	std::cout << "\n\n";


	std::cout << "2nd derivative weight: " << endl;
	double** ptr_B_FE = myFEGrid._2nd_derivatic_weight_;
	int mm = myFEGrid.TotalPointsNumber;
	for (int i = 0; i < mm; i++) {
		for (int j = 0; j < nn; j++) {
			std::cout << ptr_B_FE[i][j] << std::setw(14);
		}
		std::cout << "\n";
	}
	std::cout << "\n\n";
	
#endif
	
#if 0
	// Daughter///////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout << "Daughter class---------------------------------------------------" << endl;
	int N_d = 8;
	Grid_OC mygrid1(N_d-2, 1, 1);

	Daughter<Grid_OC> MyDaughter(&mygrid1, "2.4/vj*exp(-4.5*(2*vi-vj)^2/vj^2)", "vj", "vi", "volume_DL");
	//Daughter<Grid_OC> MyDaughter(&mygrid1, "7.2/rj^3*exp(-4.5*(2*ri^3-rj^3)^2/rj^6)*ri^2", "rj", "ri", "radius_DL");
	//Daughter<Grid_OC> MyDaughter(&mygrid1, "3*ri^2/rj^3", "rj", "ri", "radius_DL");
	//Daughter<Grid_OC> MyDaughter(&mygrid1, "Model1");
	//Daughter<Grid_OC> MyDaughter(&mygrid1);

	cout << "Formula str       :" << MyDaughter.Get_Formula() << endl;
	cout << "Variable str      :" << MyDaughter.Get_Variable() << endl;
	cout << "Variable Type str :" << MyDaughter.Get_VarType() << endl;
	cout << "Integration       :" << MyDaughter.Get_Int() << endl;

	std::cout << "Daughter distribution matrix: " << endl;
	for (int j = 0; j < N_d; j++) {
		for (int i = 0; i < N_d; i++) {
			std::cout << MyDaughter.beta[j][i] << std::setw(14);
		}
		std::cout << "\n";
	}
	std::cout << "\n\n";

	std::cout << "shidfted daughter distribution matrix: " << endl;
	for (int j = 0; j < N_d; j++) {
		for (int i = 0; i < N_d; i++) {
			std::cout << MyDaughter.beta_Shifted[j][i] << std::setw(14);
		}
		std::cout << "\n";
	}
	std::cout << "\n\n";

	

#endif

#if 0
	// Fluid///////////////////////////////////////////////////////////////////////////////////////////////////////////

	Fluid Fluid1, Fluid2;

	Fluid1.Set_Density(995);
	Fluid1.Set_Viscosity(1e-3);

	Fluid2.Set_Density(927);
	Fluid2.Set_Viscosity(156e-3);


	std::cout << "Density_Fluid1:     " << Fluid1.Density << endl;
	std::cout << "Viscosity_Fluid1:     " << Fluid1.Viscosity << endl;

	std::cout << "Density_Fluid2:     " << Fluid2.Density << endl;
	std::cout << "Viscosity_Fluid2:     " << Fluid2.Viscosity << endl;

	std::cout << "\n\n";

#endif

#if 0
	// MP_Fluid///////////////////////////////////////////////////////////////////////////////////////////////////////////

	MP_Fluid MyFluid(&Fluid1, &Fluid2, 0.3);

	MyFluid.Set_Tur_Ener_Diss_Rate_Av(0.1815);
	//MyFluid.Set_Tur_Ener_Diss_Rate_Av(0.02, 1e-4);
	
	MyFluid.Set_Surface_Tension(0.0140);
	


	
	std::cout << "Density_Dispersed:        "		<< MyFluid.Fluid_D->Density << endl;
	std::cout << "Density_Continuous:       "		<< MyFluid.Fluid_C->Density << endl;
	std::cout << "Viscosity_Dispersed:      "		<< MyFluid.Fluid_D->Viscosity << endl;
	std::cout << "Viscosity_Continuous:     "		<< MyFluid.Fluid_C->Viscosity << endl;
	std::cout << "Surface_Tension:          "		<< MyFluid.Surface_Tension << endl;
	std::cout << "Tur_Ener_Diss_Rate_Av:    "		<< MyFluid.Tur_Ener_Diss_Rate_Av << endl;
	
	std::cout << "\n\n";
	
#endif
	
#if 0
	// Coalescence///////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout << "Coalescence class---------------------------------------------------" << endl;
	double Rm1 = 2e-4;
	int N_d2 = 8;
	Grid_OC mygrid1(N_d2 - 2, 1, 1, Rm1);
	double* ptr_root_d = mygrid1.Grid_Points_;

	std::string Formula_str = "4*2^(1/3)*kc1*e^(1/3)*(ri+rj)^2*(ri^(2/3)+rj^(2/3))^0.5*exp(-kc2*ro_c^0.5*e^(1/3)/2^(1/6)/to_wo^0.5*(1/2*(1/ri+1/rj)^-1)^(5/6))";
	std::string Variable_str_i = "ri";
	std::string Variable_str_j = "rj";
	std::string VarType_str = "radius";
	std::string Parameter_str[2] = { "kc1" , "kc2" };
	double Parameter_Val[2] = { 0.004 ,900 };
	unsigned int Parameter_number = 2;
	std::string Fluid_str[3] = { "ro_c", "to_wo", "e" };
	unsigned int Fluid_number = 3;
	std::string Fluid_ID[] = { "Density_Continuous", "Surface_Tension", "Tur_Ener_Diss_Rate_Av" };

	
	Coalescence<Grid_OC> MyCoalescence(&mygrid1,  &MyFluid, "Model1", Parameter_Val);
	//Coalescence<Grid_OC> MyCoalescence(&mygrid1, &MyFluid, Formula_str, Variable_str_i, Variable_str_j, VarType_str, Parameter_str, Parameter_Val, Parameter_number, Fluid_str, Fluid_ID, Fluid_number);

	for (int j = 0; j < N_d2; j++) {
		for (int i = 0; i < N_d2; i++)
			std::cout << MyCoalescence.C_r[j][i] << "     ";
		std::cout << endl;
	}
	std::cout << endl << endl << endl;
	
	for (int j = 0; j < N_d2; j++) {
		for (int i = 0; i < N_d2; i++)
			std::cout << MyCoalescence.C_r_Shifted[i][j] << "     ";
		std::cout << endl;
	}
	std::cout << endl;
	
	
	std::cout << "Formula str       :" << MyCoalescence.Get_Formula() << endl;
	std::cout << "Variable str      :" << MyCoalescence.Get_Variable() << endl;
	std::cout << "Variable Type str :" << MyCoalescence.Get_VarType() << endl;
	std::cout << "Parameter val :" << MyCoalescence.Get_Parameter_Val()[0] << endl;
	std::cout << "Parameter number :" << MyCoalescence.Get_Parameter_number() << endl;
	//std::cout << "Fluid val :" << MyCoalescence.Get_Fluid_Val()[0] << endl;
	std::cout << "Fluid number :" << MyCoalescence.Get_Fluid_number() << endl;
	std::cout << "Parameter str :" << MyCoalescence.Get_Parameter_str() << endl;
	std::cout << "Fluid str :" << MyCoalescence.Get_Fluid_str() << endl;

	std::cout << "\n\n";


#endif

#if 0
	// Breakage///////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout << "Breakage class---------------------------------------------------" << endl;
	double Rm1 = 2e-4;
	int N_d2 = 8;
	Grid_OC mygrid1(N_d2 - 2, 1, 1, Rm1);
	
	std::string Formula_str = "kb1*e^(1/3)/2^(2/3)/r^(2/3)*(ro_d/ro_c)^0.5*exp(-kb2*to_wo/ro_d/2^(5/3)/e^(2/3)/r^(5/3))";
	std::string Variable_str = "r";
	std::string VarType_str = "radius";
	std::string Parameter_str[2] = {"kb1" , "kb2"};
	double Parameter_Val[2] = { 3e-6 ,2e-4 };
	unsigned int Parameter_number=2;
	std::string Fluid_str[4] = {"ro_d" , "ro_c", "to_wo", "e"};
	unsigned int Fluid_number = 4;
	std::string Fluid_ID[] = { "Density_Dispersed", "Density_Continuous", "Surface_Tension", "Tur_Ener_Diss_Rate_Av"};
	
	
	//Breakage<Grid_OC> MyBreakage(&mygrid1, &MyFluid, "Model1", Parameter_Val );
	Breakage<Grid_OC> MyBreakage(&mygrid1, &MyFluid, Formula_str, Variable_str, VarType_str, Parameter_str, Parameter_Val, Parameter_number, Fluid_str, Fluid_ID, Fluid_number);   
	

	for (int i = 0; i < N_d2; i++) 
		std::cout << MyBreakage.B_f[i] << "     ";
	std::cout << endl << endl << endl;

	for (int j = 0; j < N_d2; j++) {
		for (int i = 0; i < N_d2; i++)
			std::cout << MyBreakage.B_f_Shifted[i][j] << "     ";
		std::cout << endl;
	}
	std::cout << endl;


	std::cout << "Formula str       :" << MyBreakage.Get_Formula() << endl;
	std::cout << "Variable str      :" << MyBreakage.Get_Variable() << endl;
	std::cout << "Variable Type str :" << MyBreakage.Get_VarType() << endl;
	std::cout << "Parameter val :" << MyBreakage.Get_Parameter_Val()[0] << endl;
	std::cout << "Parameter number :" << MyBreakage.Get_Parameter_number() << endl;
	std::cout << "Fluid number :" << MyBreakage.Get_Fluid_number() << endl;
	std::cout << "Parameter str :" << MyBreakage.Get_Parameter_str() << endl;
	std::cout << "Fluid str :" << MyBreakage.Get_Fluid_str() << endl;

	std::cout << "\n\n";


#endif


#if 0
	// Logger///////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout << "Logger class---------------------------------------------------" << endl;
	
	Logger::Instance().Print_Error();
	Logger::Instance().Print_Warning();
	Logger::Instance().Print_Info();
	/*
	std::cout << "Error_flag:      " << Logger::Instance().Get_Error_flag() << endl;
	std::cout << "Warning_flag:    " << Logger::Instance().Get_Warning_flag() << endl;
	std::cout << "Get_Info_flag:   " << Logger::Instance().Get_Info_flag() << endl;
	
	std::cout << "Error_str:       \n" << Logger::Instance().Get_Error_str();
	std::cout << "Warning_str:     \n" << Logger::Instance().Get_Warning_str();
	std::cout << "Get_Info_str:    \n" << Logger::Instance().Get_Info_str();
	
	string e1 = "Distribution.E1				User-specified grid is outside 0 - 1.\n";
	string e2 = "Distribution.E2				User-specified distribution is negative.\n";
	string w1 = "Distribution.W1				User-specified grid is not sorted or there are duplicated points.\n";
	string w2 = "Distribution.W2				User-specified distribution has more than 1 peak.\n";

	Logger::Instance().Add_Error(e1);
	Logger::Instance().Add_Error(e2);
	Logger::Instance().Add_Warning(w1);
	Logger::Instance().Add_Warning(w2);
	Logger::Instance().Add_Info("Welecome to To PBE Class Library.");

	Logger::Instance().Print_Error();
	Logger::Instance().Print_Warning();
	Logger::Instance().Print_Info();

	std::cout << "Error_flag:      " << Logger::Instance().Get_Error_flag() << endl;
	std::cout << "Warning_flag:    " << Logger::Instance().Get_Warning_flag() << endl;
	std::cout << "Get_Info_flag:   " << Logger::Instance().Get_Info_flag() << endl;

	std::cout << "Error_str:       \n" << Logger::Instance().Get_Error_str();
	std::cout << "Warning_str:     \n" << Logger::Instance().Get_Warning_str();
	std::cout << "Get_Info_str:    \n" << Logger::Instance().Get_Info_str();
	*/
#endif
	
	
#if 1
	// 0D_solver///////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout << "0D_solver (1st kernel)---------------------------------------------------" << endl;
	std::cout << "Grid generation..." << std::endl << std::endl;
	std::cout << "Inlet-initial distribution..." << std::endl << std::endl;

	typedef Grid_FE GRID_TYPE;
	
	double ptr_root1[20] = { 0.0, 0.0042, 0.0221, 0.0537, 0.0981, 0.1542, 0.2201, 0.2941, 0.3741, 0.4576, 0.5424, 0.6259, 0.7059, 0.7799, 0.8458, 0.9019, 0.9463, 0.9779, 0.9958, 1.0 };
	double ptr_sai[20] = { 0, 0, 0, 0, 0, 0, 0, 0, .4466, 1.9116, 0.0797, 0.0003, 0, 0, 0, 0, 0, 0, 0, 0 };

	Distribution<GRID_TYPE> MyInitialDistribution(ptr_sai, ptr_root1, 20);

	//unsigned int ElementPointsNumber[] = { 10,15,15,10,5,5,5 };
	//unsigned int ElementNumber = 7;
	unsigned int ElementPointsNumber[] = { 9,8,8,9 };
	unsigned int ElementNumber = 4;
	Grid_FE mygrid1(MyInitialDistribution.Int_01, MyInitialDistribution.Mean, MyInitialDistribution.Int_99, ElementPointsNumber, ElementNumber, false, false, 200e-6);
	
	MyInitialDistribution.Set_Grid(&mygrid1);
	MyInitialDistribution.Normalize(0.2);
	
	Distribution<GRID_TYPE> MyInletDistribution(&mygrid1, -1.5, 0.1, 0.3);
	MyInletDistribution.Normalize(0.3);
	
	

	cout << "Daughter distribution..." << std::endl << std::endl;
	Daughter<GRID_TYPE> MyDaughter(&mygrid1, "7.2/rj^3*exp(-4.5*(2*ri^3-rj^3)^2/rj^6)*ri^2", "rj", "ri", "radius_DL");

	cout << "Fluid-muliPhase fluid..." << std::endl << std::endl;
	Fluid Fluid1, Fluid2;
	
	Fluid1.Set_Density(995);
	Fluid1.Set_Viscosity(1e-3);
	
	Fluid2.Set_Density(927);
	Fluid2.Set_Viscosity(156e-3);
	
	MP_Fluid MyFluid(&Fluid1, &Fluid2, 0.2);
	
	MyFluid.Set_Tur_Ener_Diss_Rate_Av(0.1815);
	MyFluid.Set_Surface_Tension(0.0140);

	
	cout << "Coalescence..." << std::endl << std::endl;
	std::string Formula_str_Co = "4*2^(1/3)*kc1*e^(1/3)*(ri+rj)^2*(ri^(2/3)+rj^(2/3))^0.5*exp(-kc2*ro_c^0.5*e^(1/3)/2^(1/6)/to_wo^0.5*(1/2*(1/ri+1/rj)^-1)^(5/6))";
	std::string Variable_str_i_Co = "ri";
	std::string Variable_str_j_Co = "rj";
	std::string VarType_str_Co = "radius";
	std::string Parameter_str_Co[2] = { "kc1" , "kc2" };
	double Parameter_Val_Co[2] = { 8e-7 ,10 };
	unsigned int Parameter_number_Co = 2;
	std::string Fluid_str_Co[3] = { "ro_c", "to_wo", "e" };
	unsigned int Fluid_number_Co = 3;
	std::string Fluid_ID_Co[] = { "Density_Continuous", "Surface_Tension", "Tur_Ener_Diss_Rate_Av" };
	Coalescence<GRID_TYPE> MyCoalescence(&mygrid1, &MyFluid, Formula_str_Co, Variable_str_i_Co, Variable_str_j_Co, VarType_str_Co, Parameter_str_Co, Parameter_Val_Co, Parameter_number_Co, Fluid_str_Co, Fluid_ID_Co, Fluid_number_Co);
	
	
	cout << "Breakage..." << std::endl << std::endl;
	std::string Formula_str_Br = "kb1*e^(1/3)/2^(2/3)/r^(2/3)*(ro_d/ro_c)^0.5*exp(-kb2*to_wo/ro_d/2^(5/3)/e^(2/3)/r^(5/3))";
	std::string Variable_str_Br = "r";
	std::string VarType_str_Br = "radius";
	std::string Parameter_str_Br[2] = { "kb1" , "kb2" };
	double Parameter_Val_Br[2] = { 3e-6 ,2e-3 };
	unsigned int Parameter_number_Br = 2;
	std::string Fluid_str_Br[4] = { "ro_d" , "ro_c", "to_wo", "e" };
	unsigned int Fluid_number_Br = 4;
	std::string Fluid_ID_Br[] = { "Density_Dispersed", "Density_Continuous", "Surface_Tension", "Tur_Ener_Diss_Rate_Av" };
	Breakage<GRID_TYPE> MyBreakage(&mygrid1, &MyFluid, Formula_str_Br, Variable_str_Br, VarType_str_Br, Parameter_str_Br, Parameter_Val_Br, Parameter_number_Br, Fluid_str_Br, Fluid_ID_Br, Fluid_number_Br);
	
	PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyDaughter, &MyBreakage, &MyCoalescence);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyDaughter, &MyBreakage);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyCoalescence);
	
	double Residence_Time = 200;
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyInletDistribution, &MyDaughter, &MyBreakage, &MyCoalescence, Residence_Time);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyInletDistribution, &MyDaughter, &MyBreakage, Residence_Time);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyInletDistribution, &MyCoalescence, Residence_Time);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyInletDistribution, &mygrid1, Residence_Time);
	
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyDaughter, &MyBreakage, &MyCoalescence, 2.0);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyDaughter, &MyBreakage, 2.0);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyCoalescence, 0.5);
	
	//runge_kutta_cash_karp54 
	//runge_kutta_dopri5     
	//runge_kutta_fehlberg78 
	//STI

	
	MyPBE_0D_Solution.Set_Solver("ODEINT", "STI", 1e-4, 1e-4, 1.0, 1.0);    
	cout << "simulation running..."<< std::endl<<std::endl;
	MyPBE_0D_Solution.Simulate(2000);

	std::cout << "preprocessing time:    " << MyPBE_0D_Solution.PrePros_Time << endl;
	std::cout << "solver time:    " << MyPBE_0D_Solution.Execution_Time << endl;
	std::cout << "step number:    " << MyPBE_0D_Solution.Step_Number << endl;
	std::cout << "initial distribution:    " << std::endl;
	for (int i = 0; i < mygrid1.TotalPointsNumber; i++)
		std::cout << MyPBE_0D_Solution.Dist_Sol[0][i] << "    ";
	std::cout << std::endl << std::endl;

	std::cout << "final distribution:    " << std::endl;
	for (int i=0; i< mygrid1.TotalPointsNumber; i++)
		std::cout << MyPBE_0D_Solution.Dist_Sol[MyPBE_0D_Solution.Step_Number-1][i] << "    ";
	std::cout << std::endl << std::endl;

	std::cout << "Grid points:    " << std::endl;
	for (int i = 0; i < mygrid1.TotalPointsNumber; i++)
		std::cout << mygrid1.Grid_Points_[i] << "    ";
	std::cout << std::endl << std::endl;
	
	std::cout << "quadrature:    " << std::endl;
	for (int i = 0; i < mygrid1.TotalPointsNumber; i++)
		std::cout << mygrid1.Quadrature_weight_[i] << "    ";
	std::cout << std::endl << std::endl;
#endif

#if 0
	// 0D_solver///////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout << "0D_solver (2nd kernel)---------------------------------------------------" << endl;
	std::cout << "Grid generation..." << std::endl << std::endl;
	std::cout << "Inlet-initial distribution..." << std::endl << std::endl;

	typedef Grid_FE GRID_TYPE;

	double ptr_root1[20] = { 0.0, 0.0042, 0.0221, 0.0537, 0.0981, 0.1542, 0.2201, 0.2941, 0.3741, 0.4576, 0.5424, 0.6259, 0.7059, 0.7799, 0.8458, 0.9019, 0.9463, 0.9779, 0.9958, 1.0 };
	double ptr_sai[20] = { 0, 0, 0, 0, 0, 0, .4466, 1.9116, 0.0797, 0.0003, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	Distribution<GRID_TYPE> MyInitialDistribution(ptr_sai, ptr_root1, 20);

	unsigned int ElementPointsNumber[] = { 10,15,15,10,5,5,5 };
	unsigned int ElementNumber = 7;
	//unsigned int ElementPointsNumber[] = { 9,8,8,9 };
	//unsigned int ElementNumber = 4;
	Grid_FE mygrid1(MyInitialDistribution.Int_01, MyInitialDistribution.Mean, MyInitialDistribution.Int_99, ElementPointsNumber, ElementNumber, false, false, 400e-6);

	MyInitialDistribution.Set_Grid(&mygrid1);
	MyInitialDistribution.Normalize(0.2);

	Distribution<GRID_TYPE> MyInletDistribution(&mygrid1, -1.5, 0.1, 0.3);
	MyInletDistribution.Normalize(0.3);



	cout << "Daughter distribution..." << std::endl << std::endl;
	Daughter<GRID_TYPE> MyDaughter(&mygrid1, "3.0*ri^2/rj^3", "rj", "ri", "radius_DL");

	cout << "Fluid-muliPhase fluid..." << std::endl << std::endl;
	Fluid Fluid1, Fluid2;

	Fluid1.Set_Density(995);
	Fluid1.Set_Viscosity(1e-3);

	Fluid2.Set_Density(927);
	Fluid2.Set_Viscosity(156e-3);

	MP_Fluid MyFluid(&Fluid1, &Fluid2, 0.2);

	MyFluid.Set_Tur_Ener_Diss_Rate_Av(0.1815);
	MyFluid.Set_Surface_Tension(0.0140);

	cout << "Coalescence..." << std::endl << std::endl;
	std::string Formula_str_Co = "cons";
	std::string Variable_str_i_Co = "ri";
	std::string Variable_str_j_Co = "rj";
	std::string VarType_str_Co = "radius";
	std::string Parameter_str_Co[] = { "cons" };
	double Parameter_Val_Co[] = { 1e-14 };
	unsigned int Parameter_number_Co = 1;
	std::string Fluid_str_Co[] = { "ro_c" };
	unsigned int Fluid_number_Co = 1;
	std::string Fluid_ID_Co[] = { "Density_Continuous", "Surface_Tension", "Tur_Ener_Diss_Rate_Av" };
	Coalescence<GRID_TYPE> MyCoalescence(&mygrid1, &MyFluid, Formula_str_Co, Variable_str_i_Co, Variable_str_j_Co, VarType_str_Co, Parameter_str_Co, Parameter_Val_Co, Parameter_number_Co, Fluid_str_Co, Fluid_ID_Co, Fluid_number_Co);


	cout << "Breakage..." << std::endl << std::endl;
	std::string Formula_str_Br = "cons*r^3";
	std::string Variable_str_Br = "r";
	std::string VarType_str_Br = "radius";
	std::string Parameter_str_Br[] = { "cons" };
	double Parameter_Val_Br[] = { 9e8};
	unsigned int Parameter_number_Br = 1;
	std::string Fluid_str_Br[] = { "ro_d" };
	unsigned int Fluid_number_Br = 1;
	std::string Fluid_ID_Br[] = { "Density_Dispersed"};
	Breakage<GRID_TYPE> MyBreakage(&mygrid1, &MyFluid, Formula_str_Br, Variable_str_Br, VarType_str_Br, Parameter_str_Br, Parameter_Val_Br, Parameter_number_Br, Fluid_str_Br, Fluid_ID_Br, Fluid_number_Br);

	PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyDaughter, &MyBreakage, &MyCoalescence);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyDaughter, &MyBreakage);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyCoalescence);

	double Residence_Time = 1000;
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyInletDistribution, &MyDaughter, &MyBreakage, &MyCoalescence, Residence_Time);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyInletDistribution, &MyDaughter, &MyBreakage, Residence_Time);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyInletDistribution, &MyCoalescence, Residence_Time);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyInletDistribution, &mygrid1, Residence_Time);

	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyDaughter, &MyBreakage, &MyCoalescence, 2.0);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyDaughter, &MyBreakage, 2.0);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyCoalescence, 0.5);

	//runge_kutta_cash_karp54 
	//runge_kutta_dopri5   
	//runge_kutta_fehlberg78 
	//STI
	MyPBE_0D_Solution.Set_Solver("ODEINT", "STI", 1e-4, 1e-4, 1.0, 1.0);
	cout << "simulation running..." << std::endl << std::endl;
	MyPBE_0D_Solution.Simulate(2000);

	std::cout << "preprocessing time:    " << MyPBE_0D_Solution.PrePros_Time << endl;
	std::cout << "solver time:    " << MyPBE_0D_Solution.Execution_Time << endl;
	std::cout << "step number:    " << MyPBE_0D_Solution.Step_Number << endl;
	std::cout << "initial distribution:    " << std::endl;
	for (int i = 0; i < mygrid1.TotalPointsNumber; i++)
		std::cout << MyPBE_0D_Solution.Dist_Sol[0][i] << "    ";
	std::cout << std::endl << std::endl;

	std::cout << "final distribution:    " << std::endl;
	for (int i = 0; i < mygrid1.TotalPointsNumber; i++)
		std::cout << MyPBE_0D_Solution.Dist_Sol[MyPBE_0D_Solution.Step_Number - 1][i] << "    ";
	std::cout << std::endl << std::endl;

	std::cout << "Grid points:    " << std::endl;
	for (int i = 0; i < mygrid1.TotalPointsNumber; i++)
		std::cout << mygrid1.Grid_Points_[i] << "    ";
	std::cout << std::endl << std::endl;

	std::cout << "quadrature:    " << std::endl;
	for (int i = 0; i < mygrid1.TotalPointsNumber; i++)
		std::cout << mygrid1.Quadrature_weight_[i] << "    ";
	std::cout << std::endl << std::endl;

#endif

#if 0
	// 0D_solver///////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout << "0D_solver (3rd kernel)---------------------------------------------------" << endl;
	std::cout << "Grid generation..." << std::endl << std::endl;
	std::cout << "Inlet-initial distribution..." << std::endl << std::endl;

	typedef Grid_FE GRID_TYPE;

	double ptr_root1[20] = { 0.0, 0.0042, 0.0221, 0.0537, 0.0981, 0.1542, 0.2201, 0.2941, 0.3741, 0.4576, 0.5424, 0.6259, 0.7059, 0.7799, 0.8458, 0.9019, 0.9463, 0.9779, 0.9958, 1.0 };
	double ptr_sai[20] = { 0, 0, 0,.4466, 1.9116, 0.0797, 0.0003, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	Distribution<GRID_TYPE> MyInitialDistribution(ptr_sai, ptr_root1, 20);

	unsigned int ElementPointsNumber[] = { 10,15,15,10,5,5,5 };
	unsigned int ElementNumber = 7;
	//unsigned int ElementPointsNumber[] = { 9,8,8,9 };
	//unsigned int ElementNumber = 4;
	Grid_FE mygrid1(MyInitialDistribution.Int_01, MyInitialDistribution.Mean, MyInitialDistribution.Int_99, ElementPointsNumber, ElementNumber, false, false, 800e-6);

	MyInitialDistribution.Set_Grid(&mygrid1);
	MyInitialDistribution.Normalize(0.2);

	Distribution<GRID_TYPE> MyInletDistribution(&mygrid1, -1.5, 0.1, 0.3);
	MyInletDistribution.Normalize(0.3);


	cout << "Daughter distribution..." << std::endl << std::endl;
	Daughter<GRID_TYPE> MyDaughter(&mygrid1, "45*ri^2*(ri/rj)^3*(1-(ri/rj)^3)/rj^3/2^(4/3)", "rj", "ri", "radius_DL");
	
	cout << "Fluid-muliPhase fluid..." << std::endl << std::endl;
	Fluid Fluid1, Fluid2;

	Fluid1.Set_Density(995);
	Fluid1.Set_Viscosity(1e-3);

	Fluid2.Set_Density(927);
	Fluid2.Set_Viscosity(156e-3);

	MP_Fluid MyFluid(&Fluid1, &Fluid2, 0.2);

	MyFluid.Set_Tur_Ener_Diss_Rate_Av(0.1815);
	MyFluid.Set_Surface_Tension(0.0140);

	cout << "Coalescence..." << std::endl << std::endl;
								 
	std::string Formula_str_Co = "4*2^(1/3)*kc1*e^(1/3)*(ri+rj)^2*(ri^(2/3)+rj^(2/3))^0.5*exp(-kc2*ro_c^0.5*e^(1/3)/2^(1/6)/to_wo^0.5*(1/2*(1/ri+1/rj)^-1)^4)";
	std::string Variable_str_i_Co = "ri";
	std::string Variable_str_j_Co = "rj";
	std::string VarType_str_Co = "radius";
	std::string Parameter_str_Co[2] = { "kc1" , "kc2" };
	double Parameter_Val_Co[2] = { 1e-5 ,0.8 };
	unsigned int Parameter_number_Co = 2;
	std::string Fluid_str_Co[3] = { "ro_c", "to_wo", "e" };
	unsigned int Fluid_number_Co = 3;
	std::string Fluid_ID_Co[] = { "Density_Continuous", "Surface_Tension", "Tur_Ener_Diss_Rate_Av" };
	Coalescence<GRID_TYPE> MyCoalescence(&mygrid1, &MyFluid, Formula_str_Co, Variable_str_i_Co, Variable_str_j_Co, VarType_str_Co, Parameter_str_Co, Parameter_Val_Co, Parameter_number_Co, Fluid_str_Co, Fluid_ID_Co, Fluid_number_Co);

	cout << "Breakage..." << std::endl << std::endl;
	std::string Formula_str_Br = "kb1*e^(1/3)/2^(2/3)/r^(2/3)*(ro_d/ro_c)^0.5*exp(-kb2*to_wo/ro_d/2^(5/3)/e^(2/3)/r^(5/3))";
	std::string Variable_str_Br = "r";
	std::string VarType_str_Br = "radius";
	std::string Parameter_str_Br[2] = { "kb1" , "kb2" };
	double Parameter_Val_Br[2] = { 3e-6 ,1e-3 };
	unsigned int Parameter_number_Br = 2;
	std::string Fluid_str_Br[4] = { "ro_d" , "ro_c", "to_wo", "e" };
	unsigned int Fluid_number_Br = 4;
	std::string Fluid_ID_Br[] = { "Density_Dispersed", "Density_Continuous", "Surface_Tension", "Tur_Ener_Diss_Rate_Av" };
	Breakage<GRID_TYPE> MyBreakage(&mygrid1, &MyFluid, Formula_str_Br, Variable_str_Br, VarType_str_Br, Parameter_str_Br, Parameter_Val_Br, Parameter_number_Br, Fluid_str_Br, Fluid_ID_Br, Fluid_number_Br);

	PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyDaughter, &MyBreakage, &MyCoalescence);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyDaughter, &MyBreakage);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyCoalescence);

	double Residence_Time = 500;
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyInletDistribution, &MyDaughter, &MyBreakage, &MyCoalescence, Residence_Time);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyInletDistribution, &MyDaughter, &MyBreakage, Residence_Time);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyInletDistribution, &MyCoalescence, Residence_Time);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyInletDistribution, &mygrid1, Residence_Time);

	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyDaughter, &MyBreakage, &MyCoalescence, 2.0);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyDaughter, &MyBreakage, 2.0);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyCoalescence, 2.0);

	//runge_kutta_cash_karp54 
	//runge_kutta_dopri5     
	//runge_kutta_fehlberg78 
	//STI

	MyPBE_0D_Solution.Set_Solver("ODEINT", "STI", 1e-4, 1e-4, 1.0, 1.0);
	cout << "simulation running..." << std::endl << std::endl;
	MyPBE_0D_Solution.Simulate(2000);

	std::cout << "preprocessing time:    " << MyPBE_0D_Solution.PrePros_Time << endl;
	std::cout << "solver time:    " << MyPBE_0D_Solution.Execution_Time << endl;
	std::cout << "step number:    " << MyPBE_0D_Solution.Step_Number << endl;
	std::cout << "initial distribution:    " << std::endl;
	for (int i = 0; i < mygrid1.TotalPointsNumber; i++)
		std::cout << MyPBE_0D_Solution.Dist_Sol[0][i] << "    ";
	std::cout << std::endl << std::endl;

	std::cout << "final distribution:    " << std::endl;
	for (int i = 0; i < mygrid1.TotalPointsNumber; i++)
		std::cout << MyPBE_0D_Solution.Dist_Sol[MyPBE_0D_Solution.Step_Number - 1][i] << "    ";
	std::cout << std::endl << std::endl;

	std::cout << "Grid points:    " << std::endl;
	for (int i = 0; i < mygrid1.TotalPointsNumber; i++)
		std::cout << mygrid1.Grid_Points_[i] << "    ";
	std::cout << std::endl << std::endl;

	std::cout << "quadrature:    " << std::endl;
	for (int i = 0; i < mygrid1.TotalPointsNumber; i++)
		std::cout << mygrid1.Quadrature_weight_[i] << "    ";
	std::cout << std::endl << std::endl;

#endif

#if 0
	// 0D_solver///////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout << "0D_solver (4th kernel)---------------------------------------------------" << endl;
	std::cout << "Grid generation..." << std::endl << std::endl;
	std::cout << "Inlet-initial distribution..." << std::endl << std::endl;

	typedef Grid_FE GRID_TYPE;

	double ptr_root1[] = { 0.0, 0.0042, 0.0221, 0.0537, 0.0981, 0.1542, 0.2201, 0.2941, 0.3741, 0.4576, 0.5424, 0.6259, 0.7059, 0.7799, 0.8458, 0.9019, 0.9463, 0.9779, 0.9958, 1.0 };
	double ptr_sai[] = { 0, 0, 0, 0,.04, 0.5, 0.8, 1.6, .3, .01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	Distribution<GRID_TYPE> MyInitialDistribution(ptr_sai, ptr_root1, 20);

	unsigned int ElementPointsNumber[] = { 10,15,15,10,5,5,5 };
	unsigned int ElementNumber = 7;
	//unsigned int ElementPointsNumber[] = { 9,8,8,9 };
	//unsigned int ElementNumber = 4;
	Grid_FE mygrid1(MyInitialDistribution.Int_01, MyInitialDistribution.Mean, MyInitialDistribution.Int_99, ElementPointsNumber, ElementNumber, false, false, 800e-6);

	MyInitialDistribution.Set_Grid(&mygrid1);
	MyInitialDistribution.Normalize(0.2);

	Distribution<GRID_TYPE> MyInletDistribution(&mygrid1, -1.5, 0.1, 0.3);
	MyInletDistribution.Normalize(0.3);



	cout << "Daughter distribution..." << std::endl << std::endl;
	
	Daughter<GRID_TYPE> MyDaughter(&mygrid1, "2*7/rj^3/(2*pi)^.5*exp(-(2*7)^2*(ri^3-0.5*rj^3)^2/rj^6/2)*3*ri^2", "rj", "ri", "radius_DL");

	cout << "Fluid-muliPhase fluid..." << std::endl << std::endl;
	Fluid Fluid1, Fluid2;

	Fluid1.Set_Density(995);
	Fluid1.Set_Viscosity(1e-3);

	Fluid2.Set_Density(927);
	Fluid2.Set_Viscosity(156e-3);

	MP_Fluid MyFluid(&Fluid1, &Fluid2, 0.2);

	MyFluid.Set_Tur_Ener_Diss_Rate_Av(0.1815);
	MyFluid.Set_Surface_Tension(0.0140);

	cout << "Coalescence..." << std::endl << std::endl;

	std::string Formula_str_Co = "kc1*0.309*(4*3.14/3)^.5*(ro_c*e/mu_c)^.5*(ri+rj)^3";
	std::string Variable_str_i_Co = "ri";
	std::string Variable_str_j_Co = "rj";
	std::string VarType_str_Co = "radius";
	std::string Parameter_str_Co[] = { "kc1" };
	double Parameter_Val_Co[] = { 9e-5 };
	unsigned int Parameter_number_Co = 1;
	std::string Fluid_str_Co[] = { "ro_c", "to_wo", "e", "mu_c" };
	unsigned int Fluid_number_Co = 4;
	std::string Fluid_ID_Co[] = { "Density_Continuous", "Surface_Tension", "Tur_Ener_Diss_Rate_Av", "Viscosity_Continuous" };
	Coalescence<GRID_TYPE> MyCoalescence(&mygrid1, &MyFluid, Formula_str_Co, Variable_str_i_Co, Variable_str_j_Co, VarType_str_Co, Parameter_str_Co, Parameter_Val_Co, Parameter_number_Co, Fluid_str_Co, Fluid_ID_Co, Fluid_number_Co);

	cout << "Breakage..." << std::endl << std::endl;
	std::string Formula_str_Br = "kb1*e^(1/3)/r^(2/3)*(erfc((1.5*kb2/(2*ro_c*e^(2/3)*r^(5/3)/to_wo))^.5)+2/3.14^.5*(1.5*kb2/(2*ro_c*e^(2/3)*r^(5/3)/to_wo))^.5*exp(-(1.5*kb2/(2*ro_c*e^(2/3)*r^(5/3)/to_wo))))";
	std::string Variable_str_Br = "r";
	std::string VarType_str_Br = "radius";
	std::string Parameter_str_Br[] = { "kb1" , "kb2" };
	double Parameter_Val_Br[] = { 3e-6 ,1e-3 };
	unsigned int Parameter_number_Br = 2;
	std::string Fluid_str_Br[] = { "ro_c", "to_wo", "e" };
	unsigned int Fluid_number_Br = 3;
	std::string Fluid_ID_Br[] = { "Density_Continuous", "Surface_Tension", "Tur_Ener_Diss_Rate_Av" };
	Breakage<GRID_TYPE> MyBreakage(&mygrid1, &MyFluid, Formula_str_Br, Variable_str_Br, VarType_str_Br, Parameter_str_Br, Parameter_Val_Br, Parameter_number_Br, Fluid_str_Br, Fluid_ID_Br, Fluid_number_Br);

	PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyDaughter, &MyBreakage, &MyCoalescence);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyDaughter, &MyBreakage);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyCoalescence);

	double Residence_Time = 1000;
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyInletDistribution, &MyDaughter, &MyBreakage, &MyCoalescence, Residence_Time);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyInletDistribution, &MyDaughter, &MyBreakage, Residence_Time);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyInletDistribution, &MyCoalescence, Residence_Time);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyInletDistribution, &mygrid1, Residence_Time);

	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyDaughter, &MyBreakage, &MyCoalescence, 2.0);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyDaughter, &MyBreakage, 2.0);
	//PBE_0D_Solution<GRID_TYPE> MyPBE_0D_Solution(&MyInitialDistribution, &MyCoalescence, 2);

	//runge_kutta_cash_karp54 
	//runge_kutta_dopri5     
	//runge_kutta_fehlberg78 
	//STI

	MyPBE_0D_Solution.Set_Solver("ODEINT", "STI", 1e-4, 1e-4, 1.0, 1.0);
	cout << "simulation running..." << std::endl << std::endl;
	MyPBE_0D_Solution.Simulate(2000);

	std::cout << "preprocessing time:    " << MyPBE_0D_Solution.PrePros_Time << endl;
	std::cout << "solver time:    " << MyPBE_0D_Solution.Execution_Time << endl;
	std::cout << "step number:    " << MyPBE_0D_Solution.Step_Number << endl;
	std::cout << "initial distribution:    " << std::endl;
	for (int i = 0; i < mygrid1.TotalPointsNumber; i++)
		std::cout << MyPBE_0D_Solution.Dist_Sol[0][i] << "    ";
	std::cout << std::endl << std::endl;

	std::cout << "final distribution:    " << std::endl;
	for (int i = 0; i < mygrid1.TotalPointsNumber; i++)
		std::cout << MyPBE_0D_Solution.Dist_Sol[MyPBE_0D_Solution.Step_Number - 1][i] << "    ";
	std::cout << std::endl << std::endl;

	std::cout << "Grid points:    " << std::endl;
	for (int i = 0; i < mygrid1.TotalPointsNumber; i++)
		std::cout << mygrid1.Grid_Points_[i] << "    ";
	std::cout << std::endl << std::endl;

	std::cout << "quadrature:    " << std::endl;
	for (int i = 0; i < mygrid1.TotalPointsNumber; i++)
		std::cout << mygrid1.Quadrature_weight_[i] << "    ";
	std::cout << std::endl << std::endl;

	/*
	std::cout << "breakage frequency:    " << std::endl;
	for (int i = 0; i < mygrid1.TotalPointsNumber; i++)
		std::cout << MyBreakage.B_f[i] << "    ";
	std::cout << std::endl << std::endl;
	*/

#endif

#if 0
/*
	state_type x = {0, 1};
	runge_kutta4< state_type > stepper;
	integrate_const(stepper, harmonic_oscillator, x, 0.0, 10.0, 0.01);
	for (int i = 0; i < 2; i++)
		std::cout << x[i] << "           ";
	std::cout << std::endl;
	

	state_type x = { 0, 1 };
	const double dt = 0.01;
	runge_kutta4< state_type > stepper;
	for (double t = 0.0; t < 10.0; t += dt){
		stepper.do_step(harmonic_oscillator, x, t, dt);
		std::cout << x[0] << "     " << x[1] << endl;
	}
	
	
	state_type x = { 0, 1 };
	runge_kutta4< state_type > stepper;
	integrate_const(stepper, [](const state_type& x, state_type& dxdt, double t) {
		dxdt[0] = x[1]; dxdt[1] = -x[0] - gam * x[1]; }
	, x, 0.0, 10.0, 0.01);
	for (int i = 0; i < 2; i++)
		std::cout << x[i] << "           ";
	std::cout << std::endl;

	

	state_type x = { 0, 1 };
	typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
	typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
	controlled_stepper_type controlled_stepper;
	integrate_adaptive(controlled_stepper, harmonic_oscillator, x, 0.0, 10.0, 0.01);
	for (int i = 0; i < 2; i++)
		std::cout << x[i] << "           ";
	std::cout << std::endl;
	

	state_type x = { 0, 1 };
	typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
	integrate_adaptive(make_controlled< error_stepper_type >(1.0e-10, 1.0e-6),
		harmonic_oscillator, x, 0.0, 10.0, 0.01);
	for (int i = 0; i < 2; i++)
		std::cout << x[i] << "           ";
	std::cout << std::endl;
	

	state_type x = { 0, 1 };
	typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
	integrate_adaptive(make_controlled(1.0e-10, 1.0e-6, error_stepper_type()),
		harmonic_oscillator, x, 0.0, 10.0, 0.01);
	for (int i = 0; i < 2; i++)
		std::cout << x[i] << "           ";
	std::cout << std::endl;
	*/

#endif	


#if 1

	//saving binary file
	
	std::string address_str = "D:\\my realm\\PhD\\C++\\rev2 final\\ongoing\\output";
	//std::string address_str = "C:\\moein laptop\\C++ rev2\\on-going\\output";

	
	/*
	const int size_a = 8;
	double a[size_a] = {8.0/3,3.0/2,4,7,9,11,34,78.2};
	const int size_a = 40;
	*/

	std::string fin_str = address_str + std::string("\\") + std::string("out_res.dat");
	std::cout << fin_str << endl << endl << endl;


	ofstream f1(fin_str, ios::binary);
	if (!f1)
		exit(-1);

	cout << "saving results..." << endl << endl;

	double a = (mygrid1.TotalPointsNumber) * 1.0;
	double b = (MyPBE_0D_Solution.Step_Number) * 1.0;
	f1.write((char*)&(a), sizeof(double));
	f1.write((char*)&(b), sizeof(double));

	for (int i = 0; i < mygrid1.TotalPointsNumber; i++) {
		f1.write((char*)&(mygrid1.Grid_Points_[i]), sizeof(double));
	}

	for (int i = 0; i < mygrid1.TotalPointsNumber; i++) {
		f1.write((char*)&(mygrid1.Quadrature_weight_[i]), sizeof(double));
	}

	for (int i = 0; i < MyPBE_0D_Solution.Step_Number; i++) {
		f1.write((char*)&(MyPBE_0D_Solution.Time_Sol[i]), sizeof(double));
	}

	for (int i = 0; i < MyPBE_0D_Solution.Step_Number; i++) {
		for(int j=0; j< mygrid1.TotalPointsNumber;j++){
			f1.write((char*)&(MyPBE_0D_Solution.Dist_Sol[i][j]), sizeof(double));
		}
	}

	for (int i = 0; i < mygrid1.TotalPointsNumber; i++) {
		f1.write((char*)&(MyInletDistribution.dist_sai_RB[i]), sizeof(double));
	}
	
	f1.close();

#endif

	}     //end of scope for dynamic arrays for cheking the destructors


	Logger::Instance().Print_Error();
	Logger::Instance().Print_Warning();
	Logger::Instance().Print_Info();

	//std::cout << 1.0 / 1e-400 << '    ' << 1e-400 / 1e-400 << std::endl;
	
	cin.get();
    return 0;
}