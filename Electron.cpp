#include "Electron.h"
#include "Constants.h"

#include <iostream>
#include <fstream>

using namespace std;

Electron::Electron(void)
{	
	//LOAD INITIAL PROPERTIES FROM CONFIGURATION FILE
	try 
	{
		boost::property_tree::ptree pt;
		boost::property_tree::ini_parser::read_ini("config.ini", pt);

		kVector << pt.get<double>("Electron.k x cartesian lab")    //x
				<< pt.get<double>("Electron.k y cartesian lab")    //y
				<< pt.get<double>("Electron.k z cartesian lab");   //z

		kVector_spherical_lab << pt.get<double>("Electron.k x spherical lab")
						  	  << pt.get<double>("Electron.k y spherical lab")
						  	  << pt.get<double>("Electron.k z spherical lab");

		kVector_cartesian_valley << pt.get<double>("Electron.k x cartesian valley")
							 	 << pt.get<double>("Electron.k y cartesian valley")
							 	 << pt.get<double>("Electron.k z cartesian valley");

		kVector_spherical_valley << pt.get<double>("Electron.k x spherical valley")
							 	 << pt.get<double>("Electron.k y spherical valley")
							 	 << pt.get<double>("Electron.k z spherical valley");

		positionVector	<< pt.get<double>("Electron.Position x cartesian lab")    //x
						<< pt.get<double>("Electron.Position y cartesian lab")    //y
						<< pt.get<double>("Electron.Position z cartesian lab");   //z 

		positionVector_spherical_lab << pt.get<double>("Electron.Position x spherical lab")
						  			 << pt.get<double>("Electron.Position y spherical lab")
						  			 << pt.get<double>("Electron.Position z spherical lab");

		positionVector_cartesian_valley << pt.get<double>("Electron.Position x cartesian valley")
							 			<< pt.get<double>("Electron.Position y cartesian valley")
							 			<< pt.get<double>("Electron.Position z cartesian valley");

		positionVector_spherical_valley << pt.get<double>("Electron.Position x spherical valley")
							 			<< pt.get<double>("Electron.Position y spherical valley")
							 			<< pt.get<double>("Electron.Position z spherical valley");

		groupVelocityVector << pt.get<double>("Electron.Group velocity x cartesian lab")
							<< pt.get<double>("Electron.Group velocity y cartesian lab")
							<< pt.get<double>("Electron.Group velocity z cartesian lab");

		groupVelocityVector_spherical_lab << pt.get<double>("Electron.Group velocity x spherical lab")
						  				  << pt.get<double>("Electron.Group velocity y spherical lab")
						  				  << pt.get<double>("Electron.Group velocity z spherical lab");

		groupVelocityVector_cartesian_valley << pt.get<double>("Electron.Group velocity x cartesian valley")
							 				 << pt.get<double>("Electron.Group velocity y cartesian valley")
							 				 << pt.get<double>("Electron.Group velocity z cartesian valley");

		groupVelocityVector_spherical_valley << pt.get<double>("Electron.Group velocity x spherical valley")
							 				 << pt.get<double>("Electron.Group velocity y spherical valley")
							 				 << pt.get<double>("Electron.Group velocity z spherical valley");

		kineticEnergy =  pt.get<double>("Electron.Kinetic Energy");

		AfterEkVector << pt.get<double>("Electron.k x cartesian lab after E")    //x
					  << pt.get<double>("Electron.k y cartesian lab after E")    //y
					  << pt.get<double>("Electron.k z cartesian lab after E");   //z

		AfterEkVector_spherical_lab << pt.get<double>("Electron.k x spherical lab after E")
									<< pt.get<double>("Electron.k y spherical lab after E")
									<< pt.get<double>("Electron.k z spherical lab after E");

		AfterEkVector_cartesian_valley << pt.get<double>("Electron.k x cartesian valley after E")
						   			   << pt.get<double>("Electron.k y cartesian valley after E")
						   			   << pt.get<double>("Electron.k z cartesian valley after E");

		AfterEkVector_spherical_valley << pt.get<double>("Electron.k x spherical valley after E")
						   			   << pt.get<double>("Electron.k y spherical valley after E")
						   			   << pt.get<double>("Electron.k z spherical valley after E");				

		AfterEpositionVector	<< pt.get<double>("Electron.Position x cartesian lab after E")    
								<< pt.get<double>("Electron.Position y cartesian lab after E")    
								<< pt.get<double>("Electron.Position z cartesian lab after E"); 

		AfterEpositionVector_spherical_lab << pt.get<double>("Electron.Position x spherical lab after E")
										   << pt.get<double>("Electron.Position y spherical lab after E")
										   << pt.get<double>("Electron.Position z spherical lab after E");

		AfterEpositionVector_cartesian_valley << pt.get<double>("Electron.Position x cartesian valley after E")
						      				  << pt.get<double>("Electron.Position y cartesian valley after E")
						      				  << pt.get<double>("Electron.Position z cartesian valley after E");

		AfterEpositionVector_spherical_valley << pt.get<double>("Electron.Position x spherical valley after E")
						      				  << pt.get<double>("Electron.Position y spherical valley after E")
						      				  << pt.get<double>("Electron.Position z spherical valley after E");

		AfterEgroupVelocityVector << pt.get<double>("Electron.Group velocity x cartesian lab after E")
								  << pt.get<double>("Electron.Group velocity y cartesian lab after E")
								  << pt.get<double>("Electron.Group velocity z cartesian lab after E");

		AfterEgroupVelocityVector_spherical_lab << pt.get<double>("Electron.Group velocity x spherical lab after E")
											    << pt.get<double>("Electron.Group velocity y spherical lab after E")
											    << pt.get<double>("Electron.Group velocity z spherical lab after E");

		AfterEgroupVelocityVector_cartesian_valley << pt.get<double>("Electron.Group velocity x cartesian valley after E")
						     					   << pt.get<double>("Electron.Group velocity y cartesian valley after E")
						     					   << pt.get<double>("Electron.Group velocity z cartesian valley after E");

		AfterEgroupVelocityVector_spherical_valley << pt.get<double>("Electron.Group velocity x spherical valley after E")
						     					   << pt.get<double>("Electron.Group velocity y spherical valley after E")
						     					   << pt.get<double>("Electron.Group velocity z spherical valley after E"); 

		AfterEkineticEnergy =  pt.get<double>("Electron.Kinetic Energy After E");

		AfterScatterkVector << pt.get<double>("Electron.k x cartesian lab after Scatter")
							<< pt.get<double>("Electron.k y cartesian lab after Scatter")
							<< pt.get<double>("Electron.k z cartesian lab after Scatter"); 

		AfterScatterkVector_spherical_lab << pt.get<double>("Electron.k x spherical lab after Scatter")
										  << pt.get<double>("Electron.k y spherical lab after Scatter")
										  << pt.get<double>("Electron.k z spherical lab after Scatter");

		AfterScatterkVector_cartesian_valley << pt.get<double>("Electron.k x cartesian valley after Scatter")
						     				 << pt.get<double>("Electron.k y cartesian valley after Scatter")
						     				 << pt.get<double>("Electron.k z cartesian valley after Scatter");

		AfterScatterkVector_spherical_valley << pt.get<double>("Electron.k x spherical valley after Scatter")
						     				 << pt.get<double>("Electron.k y spherical valley after Scatter")
						     				 << pt.get<double>("Electron.k z spherical valley after Scatter");  

		AfterScattergroupVelocityVector << pt.get<double>("Electron.Group velocity x cartesian lab after Scatter")
										<< pt.get<double>("Electron.Group velocity y cartesian lab after Scatter")
										<< pt.get<double>("Electron.Group velocity z cartesian lab after Scatter");

		AfterScattergroupVelocityVector_spherical_lab << pt.get<double>("Electron.Group velocity x spherical lab after Scatter")
													  << pt.get<double>("Electron.Group velocity y spherical lab after Scatter")
													  << pt.get<double>("Electron.Group velocity z spherical lab after Scatter");

		AfterScattergroupVelocityVector_cartesian_valley << pt.get<double>("Electron.Group velocity x cartesian valley after Scatter")
						     							 << pt.get<double>("Electron.Group velocity y cartesian valley after Scatter")
						     							 << pt.get<double>("Electron.Group velocity z cartesian valley after Scatter");

		AfterScattergroupVelocityVector_spherical_valley << pt.get<double>("Electron.Group velocity x spherical valley after Scatter")
						     							 << pt.get<double>("Electron.Group velocity y spherical valley after Scatter")
						     							 << pt.get<double>("Electron.Group velocity z spherical valley after Scatter");

		AfterScatterkineticEnergy =  pt.get<double>("Electron.Kinetic Energy After Scatter");

		E_velocity << 0.
			  << 0.
			  << 0.;
			  
		Scatter_velocity << 0.
			  << 0.
			  << 0.;

		Initial_velocity << 0.
			    << 0.
				<< 0.;

	}
	catch (exception e)
	{
		cout << "Problem reading config file:" << endl;
		cout << e.what() << endl;
	}

	cout << "Electron Created" << endl;
}


Electron::~Electron(void)
{
	cout << "Electron Destroyed" << endl;
}

//BOOK-KEEPING FUNCTIONS
void Electron::PrintToConsole()
{
	positionVector.print("Electron Position Vector: ");
	cout << endl;
	kVector.print("Electron k Vector: ");
	cout << endl;
}

