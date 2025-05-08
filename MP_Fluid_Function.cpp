#include "MP_Fluid_Declaration.h"
#include "Logger_Declaration.h"

MP_Fluid::MP_Fluid(const Fluid* Fluid1, const Fluid* Fluid2, double Vol_Frac, std::string Phase_Inversion_Model_Str) {
	if ((Vol_Frac > 0) && (Vol_Frac < 1))
		this->Vol_Frac = Vol_Frac;
	else
		Logger::Instance().Add_Error("MP_Fluid.E1	     			Specified volume fraction is out of range.\n");

	if (Fluid1->Temperature != Fluid2->Temperature)
		Logger::Instance().Add_Warning("MP_Fluid.W1	     			Specified fluids do not have same temperatures.\n");

	if (Phase_Inversion_Model_Str == "Model1") {
		Phase_Inversion_Model1();
		if (Continuous_Flag == 1) {
			this->Fluid_C = Fluid1;
			this->Fluid_D = Fluid2;
		}
		else {
			this->Fluid_C = Fluid2;
			this->Fluid_D = Fluid1;
		}
	}

	else
		Logger::Instance().Add_Error("MP_Fluid.E2	     			Specified Phase_Inversion model is unknown.\n");

	//Set_Density_Mix();        
	//Set_Viscosity_Mix();
}

void MP_Fluid::Set_Tur_Ener_Diss_Rate_Av(double Tur_Ener_Diss_Rate_Av) { this->Tur_Ener_Diss_Rate_Av = Tur_Ener_Diss_Rate_Av; }
void MP_Fluid::Set_Tur_Ener_Diss_Rate_Av(double Power, double Vol) {
	Tur_Ener_Diss_Rate_Av = Power / (Fluid_C->Density * Vol);
}

void MP_Fluid::Set_Surface_Tension(double Surface_Tension) { this->Surface_Tension = Surface_Tension; }

void MP_Fluid::Phase_Inversion_Model1() {
	if (Vol_Frac > 0.5)
		Continuous_Flag = 1;
	else
		Continuous_Flag = 2;
}