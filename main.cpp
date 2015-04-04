#include <iostream>

//My Classes
#include "Simulator.h"

using namespace std;

int main(int argc, char* argv[])
{

	for (int i = 1; i <= 10000; i++)
	{
		Simulator* simulator = new Simulator();
		simulator->EFIeldVector << atoi(argv[1]) << atoi(argv[2]) << atoi(argv[3]);
		simulator->RunElectronSimulation();
		simulator->~Simulator();
	}

	//WAIT FOR USER TO PRESS ENTER BEFORE ENDING PROGRAM
	cout << "PRESS ENTER TO QUIT\n";
	cin.get(); //wait for user to press enter before quitting
	return 0;
} 