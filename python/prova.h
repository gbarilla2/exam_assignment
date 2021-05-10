#include <iostream>
#include "Math/Vector4D.h"
#include "ROOT/RVec.hxx"
    
using namespace ROOT::VecOps; 
using rvec_f = RVec<float>;
const auto Z_mass = 91.2;


//This function return the invariant mass of 2e 2mu
float invariant_mass_2el2mu(rvec_f el_pt, rvec_f el_eta, rvec_f el_phi, rvec_f el_mass, rvec_f mu_pt, rvec_f mu_eta,
                                                rvec_f mu_phi, rvec_f mu_mass)
    {
    ROOT::Math::PtEtaPhiMVector p1(mu_pt[0], mu_eta[0], mu_phi[0], mu_mass[0]);
    ROOT::Math::PtEtaPhiMVector p2(mu_pt[1], mu_eta[1], mu_phi[1], mu_mass[1]);
    ROOT::Math::PtEtaPhiMVector p3(el_pt[0], el_eta[0], el_phi[0], el_mass[0]);
    ROOT::Math::PtEtaPhiMVector p4(el_pt[1], el_eta[1], el_phi[1], el_mass[1]);
    
	return (p1 + p2 + p3 + p4).M();
    }

//This function return the mass of Z bosons ordered through the lepton pair that closest to Z mass
RVec<float> calculation_Z_mass_2el2mu(rvec_f el_pt, rvec_f el_eta, rvec_f el_phi, rvec_f el_mass, rvec_f mu_pt, rvec_f mu_eta,
                                                rvec_f mu_phi, rvec_f mu_mass)
    {
	RVec<float> Z_mass_vec(2);
    ROOT::Math::PtEtaPhiMVector p1(mu_pt[0], mu_eta[0], mu_phi[0], mu_mass[0]);
    ROOT::Math::PtEtaPhiMVector p2(mu_pt[1], mu_eta[1], mu_phi[1], mu_mass[1]);
    ROOT::Math::PtEtaPhiMVector p3(el_pt[0], el_eta[0], el_phi[0], el_mass[0]);
    ROOT::Math::PtEtaPhiMVector p4(el_pt[1], el_eta[1], el_phi[1], el_mass[1]);
    auto Muons_m = (p1+p2).M();
	auto Electrons_m = (p3+p4).M();
	if (std::abs(Muons_m-Z_mass)>std::abs(Electrons_m-Z_mass)){
		Z_mass_vec[0] = Electrons_m;
		Z_mass_vec[1] = Muons_m; 
	}
	else{
		Z_mass_vec[1] = Electrons_m;
		Z_mass_vec[0] = Muons_m;
	}
		
	
	return Z_mass_vec;
    }


//This function return the invariant mass of 4l
float invariant_mass_4l(rvec_f pt, rvec_f eta, rvec_f phi, rvec_f mass)
	{
	ROOT::Math::PtEtaPhiMVector p1(pt[0], eta[0], phi[0], mass[0]);
	ROOT::Math::PtEtaPhiMVector p2(pt[1], eta[1], phi[1], mass[1]);
	ROOT::Math::PtEtaPhiMVector p3(pt[2], eta[2], phi[2], mass[2]);
	ROOT::Math::PtEtaPhiMVector p4(pt[3], eta[3], phi[3], mass[3]);
	
	return (p1+p2+p3+p4).M();
	}
	
	
//This function return the mass of Z bosons ordered through the lepton pair that closest to Z mass
RVec<float> calculation_Z_mass_4l(rvec_f pt, rvec_f eta, rvec_f phi, rvec_f mass, rvec_f charge)
    {
	RVec<float> Z_mass_vec(2);
	RVec<float> pos_charge(2);
	RVec<float> neg_charge(2);
	RVec<float> Z_masses_temp(4);
	RVec<float> temp_masses(4);
	int i=0, j=0, index_min=-1;
	ROOT::Math::PtEtaPhiMVector p1(pt[0], eta[0], phi[0], mass[0]);
	ROOT::Math::PtEtaPhiMVector p2(pt[1], eta[1], phi[1], mass[1]);
	ROOT::Math::PtEtaPhiMVector p3(pt[2], eta[2], phi[2], mass[2]);
	ROOT::Math::PtEtaPhiMVector p4(pt[3], eta[3], phi[3], mass[3]);
	RVec<ROOT::Math::PtEtaPhiMVector> p{p1,p2,p3,p4};
	//This save the positive and negative index
	for (int k=0; k<4 ;k++) {
		if( charge[k] ==1 ){
			pos_charge[i]=k;
			i++;
		}
		else{
			neg_charge[j]=k;
			j++;
		}
	}
	//This save the Z masses possible combination with two lepton of the opposite charge
	Z_masses_temp[0]=(p[pos_charge[0]]+p[neg_charge[0]]).M();
	Z_masses_temp[1]=(p[pos_charge[1]]+p[neg_charge[1]]).M();	
	Z_masses_temp[2]=(p[pos_charge[0]]+p[neg_charge[1]]).M();
	Z_masses_temp[3]=(p[pos_charge[1]]+p[neg_charge[0]]).M();
	
	temp_masses[0]=std::abs(Z_masses_temp[0]-Z_mass);
	temp_masses[1]=std::abs(Z_masses_temp[1]-Z_mass);
	temp_masses[2]=std::abs(Z_masses_temp[2]-Z_mass);
	temp_masses[3]=std::abs(Z_masses_temp[3]-Z_mass);
	//This give the index of the best mass for Z
	index_min = ArgMin(temp_masses);
	if (index_min == 0 || index_min == 1){
		Z_mass_vec[0] = Z_masses_temp[index_min];
		Z_mass_vec[1] = Z_masses_temp[(index_min+1)%2]; 
	}
	else{
		Z_mass_vec[0] = Z_masses_temp[index_min];
		Z_mass_vec[1] = Z_masses_temp[(index_min+1)%2+2];
	}
	return Z_mass_vec;
 }	
	
	
	