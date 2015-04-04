#pragma once
#define _USE_MATH_DEFINES
#include <math.h>


//FUNDAMENTAL CONSTANTS
const double e = 1.60218e-19;   //charge of an electron (Coulombs)
const double hbar = 1.05457e-34;	//h-bar (m^2 kg / s)
const double m_e = 9.10938e-31;	// electron mass in kg 


// GERMANIUM CONSTANTS
const double m_l = 1.57 * m_e;	 // in units of me, from Fischetti 1991. 
const double m_t = 0.081 * m_e;	 // from Fischetti 1991. 
const double m_G = 0.035 * m_e;	 // from Fischetti 1991. 
const double m_lX = 1.353 * m_e; // from Fischetti 1991. 
const double m_tX = 0.288 * m_e; // from Fischetti 1991. 
const double rho = 5320; //kg / m^3 (density of germanium)
const double v_sl = 5.4 * pow(10,3); // m/s longitudonal speed of sound
const double v_st = 3.2 * pow(10,3); // m/s average transverse speed of sound
const double Theta_d = -4.43 * e; //dilational deformation potential
const double Theta_u = 16.8 * e; // uniaxial shear deformation potential
const double Theta_o = 8.811e-9; // optical deformation potential

//SIMULATION SPECIFIC CONSTANTS
//const double Gamma_0 = 4 * M_PI * pow(10,10); //4 * pi * 10^10, pow(a,b) = a^b
const double N_q = 0; // average thermal occupancy of phonon mode
const double w_intervalley = 3.63865e13;
const double w_optical = 5.62737e13; //5.62067951867e13;
const double w_slow = 1.57093e13;
const double hw_slow = hbar * w_slow;
const double hw_intervalley = hbar * w_intervalley;
const double hw_optical = hbar * w_optical; // 5.9274e-21; // intravalley optical phonon energy

