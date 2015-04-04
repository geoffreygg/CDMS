#pragma once
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include "armadillo"
using namespace arma;

class Electron
{


public:
	//CONSTRUCTOR AND DESTRUCTOR
	Electron();
	~Electron(void);


	//PROPERTIES
	vec kVector; // k-vector (momentum vector)
	vec kVector_spherical_lab;
	vec kVector_cartesian_valley;
	vec kVector_spherical_valley;
	vec groupVelocityVector; // v_g = hbar  [m*^-1]  k  (bloch wave group velocity)
	vec groupVelocityVector_spherical_lab;
	vec groupVelocityVector_cartesian_valley;
	vec groupVelocityVector_spherical_valley;
	vec positionVector; // position vector
	vec positionVector_spherical_lab;
	vec positionVector_cartesian_valley;
	vec positionVector_spherical_valley;
	double kineticEnergy; // kinetic energy of the electron 
	vec AfterEkVector;
	vec AfterEkVector_spherical_lab;
	vec AfterEkVector_cartesian_valley;
	vec AfterEkVector_spherical_valley; 
	vec AfterEgroupVelocityVector; 
	vec AfterEgroupVelocityVector_spherical_lab;
	vec AfterEgroupVelocityVector_cartesian_valley;
	vec AfterEgroupVelocityVector_spherical_valley;
	vec AfterEpositionVector;
	vec AfterEpositionVector_spherical_lab;
	vec AfterEpositionVector_cartesian_valley;
	vec AfterEpositionVector_spherical_valley;
	double AfterEkineticEnergy; 
	vec AfterScatterkVector; 
	vec AfterScatterkVector_spherical_lab;
	vec AfterScatterkVector_cartesian_valley;
	vec AfterScatterkVector_spherical_valley;
	vec AfterScattergroupVelocityVector; 
	vec AfterScattergroupVelocityVector_spherical_lab;
	vec AfterScattergroupVelocityVector_cartesian_valley;
	vec AfterScattergroupVelocityVector_spherical_valley;
	double AfterScatterkineticEnergy; 
	vec E_velocity;
	vec Initial_velocity;
	vec Scatter_velocity;

	//BOOK-KEEPING FUNCTIONS
	void PrintToConsole();

};

