#include "Simulator.h"

using namespace std;

Simulator::Simulator(void)
{
	//LOAD PROPERTIES FROM CONFIGURATION FILE
	try 
	{
		boost::property_tree::ptree pt;
		boost::property_tree::ini_parser::read_ini("config.ini", pt);

		Test_Config << pt.get<double>("Simulator.Test x")    //x
					<< pt.get<double>("Simulator.Test y")    //y
					<< pt.get<double>("Simulator.Test z");   //z

		dataFilePath = pt.get<string>("Program.Data File Path");

		numberOfIterations = pt.get<int>("Program.Number Of Iterations");
	}
	catch (exception e)
	{
		cout << "Problem reading config file:" << endl;
		cout << e.what() << endl;
	}
	
	dataFileStream.open(dataFilePath, std::ofstream::app);


	//Seed Mersenne Twister Random Number Generator With Current System Time
	int seed = time(0);
	gen.seed(seed);

	//Initialize Simulation Parameters
	numberOfScatters = 0;
	globalTime = 0;


	//Print Information to Data text file
	//for (int i = 0; i < 3; i++)
	//	dataFileStream << "////////////////////////////////////////////" << endl;
	//dataFileStream << endl;
	//dataFileStream << "Simulation Begun: " << GetCurrentDateTime() << endl << endl;
	//dataFileStream << "Mersenne Twister Seed: " << seed << endl;
	//dataFileStream << "Electric Field Vector:" << endl;
	//dataFileStream << EFIeldVector << endl;

	//cout << "Beginning Simulation" << endl;
}

Simulator::~Simulator(void)
{

	//dataFileStream << "END DATA" << endl;
	//dataFileStream << "Total Number Of Physical Scatters: " << numberOfScatters << endl;
	//dataFileStream << "Simulation Ended: " << GetCurrentDateTime() << endl << endl;
	//dataFileStream << "Simulation Ended Successfully" << endl;
	dataFileStream.close(); 
	//cout << "End Simulation" << endl;
}


//SIMULATIONS
void Simulator::RunElectronSimulation()
{
	//cout << "Running electron simulation" << endl;

	Crystal* crystal = new Crystal();
	Electron* electron = new Electron();

	Gamma_0 = pow(10,12);
	Gamma_0_angular = pow(10,12);
	Gamma_0_new = Gamma_0;
	time_between_physical_scatters = 0.;
	Phys_self_scatter_check = 0.;
	physical_scatter_check = 0.;
	Random_Valley_Number = GetRandomValley();
	//ITERATE OVER THE REQUESTED NUMBER OF CYCLES
	for (int i = 1; i <= numberOfIterations; i++)
	{
		if (i % 1000 == 0)
			cout << "Step: " << i << endl;

		double t_c = - log( GetRandomNumber() ) / Gamma_0; //in c++, log = ln (natural logarithm, base e).

		ApplyElectricField(electron, crystal, t_c, Random_Valley_Number); 

		globalTime += t_c;

		Scatter(electron, crystal, Random_Valley_Number, i, t_c);

	}
	cout << "\n ratio: " << Phys_self_scatter_check << " " << physical_scatter_check << endl;
	electron->~Electron();

	//cout << "Electron Simulation Complete" << endl;
}

