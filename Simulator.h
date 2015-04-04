#pragma once

#include <random> 
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/random/mersenne_twister.hpp>

#include "armadillo"

#include <iostream>
#include <ostream>
#include <time.h>

#include "UsefulMath.h"
#include "Constants.h"
#include "Electron.h"
#include "Crystal.h"

class Simulator
{
public:
	Simulator(void);
	~Simulator(void);

	void RunElectronSimulation();
	double GetRandomNumber(); // Returns a Random Number Between 0 and 1
	int GetRandomValley();
	arma::vec EFIeldVector; //Electric Field 


private:
	//PROGRAM VALUES
	int numberOfIterations;
	int numberOfScatters;
	int Random_Valley_Number;
	int Random_old_Valley_Number;
	int physical_scatter_check;
	int Phys_self_scatter_check;
	double globalTime;
	double Gamma_0;
	double Gamma_0_new;
	double Gamma_0_angular;
	double Gamma_Fraction;
	double sumation;
	double ratio;
	double time_between_physical_scatters;
	double LessThan;
	double GreaterThan;
	double sumLessThan;
	double sumGreaterThan;
	//PHYSICAL VALUES
	arma::vec Test_Config;
	
	//RANDOM NUMBER GENERATION OBJECTS
	std::uniform_real_distribution<double> u01; //Specify uniform distribution between 0 and 1
	boost::mt19937 gen; //Mersenne Twister Random Number Generator

	//TOOLS FOR RECORDING DATA IN A TEXT FILE
	std::ofstream dataFileStream;
	std::string dataFilePath;

	//PHYSICS FUNCTIONS
	void ApplyElectricField(Electron* electron, Crystal* crystal, double timeStep, int valley_number);
	void Scatter(Electron* electron, Crystal* crystal, int valley_number, int stepNumber, double timeStep);

	//BOOK KEEPING FUNCTIONS
	void RecordData(Electron* electron, int stepNumber, double timeStep, int valley_number);
	void Initializer(Electron* electron);
	void Converter(Electron* electron, Crystal* crystal, int valley_number);
	const std::string GetCurrentDateTime();

};

