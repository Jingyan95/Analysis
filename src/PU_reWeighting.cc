#include "PU_reWeighting.h"

PU::PU(){}

double PU::PU_2016(int NumTrueInteraction, TString type_str){
if (NumTrueInteraction < 0 || NumTrueInteraction > 75 ) {return 1;}
else {
	if (type_str == "nominal") {return data_nominal_2016[NumTrueInteraction];}
	else if (type_str == "up") {return data_up_2016[NumTrueInteraction];}
	else if (type_str == "down") {return data_down_2016[NumTrueInteraction];}
	else {std::cout<<"Error pu string!"<<std::endl; return 1.0;}
}
}
		
double PU::PU_2017(int NumTrueInteraction, TString type_str){
if (NumTrueInteraction < 0 || NumTrueInteraction > 100 ) {return 1;}
else {
	if (type_str == "nominal") {return data_nominal_2017[NumTrueInteraction];}
	else if (type_str == "up") {return data_up_2017[NumTrueInteraction];}
	else if (type_str == "down") {return data_down_2017[NumTrueInteraction];}
	else {std::cout<<"Error pu string!"<<std::endl; return 1.0;}
}
}
		
double PU::PU_2018(int NumTrueInteraction, TString type_str){
if (NumTrueInteraction < 0 || NumTrueInteraction > 100 ) {return 1;}
else {
	if (type_str == "nominal") {return data_nominal_2018[NumTrueInteraction];}
	else if (type_str == "up") {return data_up_2018[NumTrueInteraction];}
	else if (type_str == "down") {return data_down_2018[NumTrueInteraction];}
	else {std::cout<<"Error pu string!"<<std::endl; return 1.0;}
}
}
				
