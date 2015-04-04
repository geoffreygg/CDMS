#include "Simulator.h"

using namespace std;

//
//PHYSICS FUNCTIONS
//
void Simulator::ApplyElectricField(Electron* electron, Crystal* crystal, double timeStep, int valley_number)
{
	mat inverseMassTensor = crystal->getInverseMassTensor(valley_number);

	electron->AfterEkVector = electron->kVector + e * EFIeldVector * timeStep / hbar;
	electron->AfterEgroupVelocityVector = hbar * inverseMassTensor * electron->AfterEkVector;
	electron->AfterEpositionVector = electron->positionVector + hbar * inverseMassTensor * electron->kVector * timeStep - (1./2.) * e * inverseMassTensor * EFIeldVector * timeStep * timeStep;
	electron->AfterEkineticEnergy = (double)(pow(hbar,2.) / 2.) * as_scalar(electron->AfterEkVector.t() * inverseMassTensor * electron->AfterEkVector);	
}

void Simulator::Scatter(Electron* electron, Crystal* crystal, int valley_number, int stepNumber, double timeStep)
{
	//VECTOR BASED MONTE CARLO SCATTERING PROCESS


	//ROTATE ELECTRON K VECTOR INTO THE VALLEY FRAME, AND CALL THIS k INITIAL. 
	arma::vec k_initial_vector_cartesian = crystal->getRotationMatrix(valley_number) * electron->AfterEkVector;

	//CONVERT k_i TO SPHERICAL COORDINATES
	arma::vec k_initial_vector_spherical = k_initial_vector_cartesian; 
	ConvertCartesianToSpherical(&k_initial_vector_spherical);

	//INITIAL COMPONENTS OF K VECTOR
	double k_i = k_initial_vector_spherical.at(0);
	double theta_i = k_initial_vector_spherical.at(1);
	double phi_i = k_initial_vector_spherical.at(2);
	

	//SELECT PHONON ANGLES
	//double theta_q = acos (1. - 2. * GetRandomNumber());   // isotropic angles for predicting rates
	double rand1 = GetRandomNumber();
	double theta_q = 0.29556502287075714*rand1 + 62.65847903378649*pow(rand1,2) - 1859.3439722987223*pow(rand1,3) + 
				32468.17885762252*pow(rand1,4) - 373482.7175293871*pow(rand1,5) + 3.039500315457711e6*pow(rand1,6) - 
				1.823753321756405e7*pow(rand1,7) + 8.285765654732832e7*pow(rand1,8) - 2.8997527271567386e8*pow(rand1,9) + 
				7.898852260950013e8*pow(rand1,10) - 1.683039544892932e9*pow(rand1,11) + 2.805250499399089e9*pow(rand1,12) - 
				3.639328852307488e9*pow(rand1,13) + 3.6347929486639495e9*pow(rand1,14) - 2.742224265546424e9*pow(rand1,15) + 
				1.5148328428419697e9*pow(rand1,16) - 5.816096377869084e8*pow(rand1,17) + 1.411916370376998e8*pow(rand1,18) - 
				1.763614580543882e7*pow(rand1,19) + 543755.4421268073*pow(rand1,20);
	if(theta_q > M_PI)
	{
		theta_q = M_PI - (0.0001)*GetRandomNumber();
	}
		
	if(theta_q < 0)
	{
		theta_q = 0.0001*GetRandomNumber();
	}
	
	double phi_q = 2. * M_PI * GetRandomNumber();

	//CREATE MATRIX STORING ALL POSSIBLE SCATTER TYPES AND THEIR PROPERTIES. 
	//COMPONENTS: a, b, c, q, dGamma
	//NUMBER OF ROWS = NUMBER OF SCATTERING PROCESSES TO CONSIDER
	
	arma::mat scatteringProcesses;	
	arma::mat scatteringProcess;
	double a, b, c, q, dGamma, z, x, y;

	double dGammaLongitudinal, dGammaTransverse, dGammaTranInterSlow1, dGammaTranInterSlow2, dGammaIntervalley1, dGammaIntervalley2, dGammaOptical1, dGammaOptical2;

///
///ACOUSTIC (LONGITUDINAL/TRANSVERSE) INTRAVALLEY EMISSIONS
///

	//
	// INTRAVALLEY ACOUSTIC LONGITUDINAL PHONON EMISSION	
	//

	a = pow(hbar, 2) * (m_l + m_t + (-m_l + m_t) * cos(2. * theta_q)) / (4. * m_l * m_t); //pow(x, y) = x^y
	b = hbar * v_sl - pow(hbar, 2) * k_i * cos(theta_i) * cos(theta_q) / m_l - pow(hbar, 2) * k_i * cos(phi_i - phi_q) * sin(theta_i) * sin(theta_q) / m_t;
	c = 0.;
	q = k_i*m_t*cos(theta_i)*cos(theta_q)*pow(m_t*pow(cos(theta_q),2) + m_l*pow(sin(theta_q),2),-1) + 
							pow(hbar,-1)*(-(m_l*m_t*v_sl*pow(m_t*pow(cos(theta_q),2) + m_l*pow(sin(theta_q),2),-1)) + 
							(pow(m_t*pow(cos(theta_q),2) + m_l*pow(sin(theta_q),2),-1)*
							pow(-4*(hbar*m_l*pow(k_i,2)*pow(sin(theta_i),2) - hbar*m_l*pow(k_i,2)*pow(cos(phi_i),2)*pow(sin(theta_i),2) - 
							hbar*m_l*pow(k_i,2)*pow(sin(phi_i),2)*pow(sin(theta_i),2))*
							(hbar*m_t*pow(cos(theta_q),2) + hbar*m_l*pow(sin(theta_q),2)) + 
							pow(2*m_l*m_t*v_sl - 2*hbar*k_i*m_t*cos(theta_i)*cos(theta_q) - 2*hbar*k_i*m_l*cos(phi_i - phi_q)*sin(theta_i)*sin(theta_q),2),0.5))/2.) + 
							k_i*m_l*cos(phi_i - phi_q)*pow(m_t*pow(cos(theta_q),2) + m_l*pow(sin(theta_q),2),-1)*sin(theta_i)*sin(theta_q);
	
	if (!(q > 0.))
	{
		q = 0.; //SET q TO ZERO IF NEGATIVE OR IMAGINARY
		a = 0.;
		b = 0.;
		c = 0.;
	}

	dGamma = 1/.4 * sin(theta_q) * pow( Theta_d + Theta_u * pow(cos(theta_q) , 2), 2) * (N_q + 1.) * pow(q, 3) / (8. * pow(M_PI,2) * rho * v_sl * abs(2. * a * q + b));

	if (!(dGamma > 0.))
	{
		dGamma = -12345.; // nonsense value for non-physical scatter
	}
	
	dGammaLongitudinal = dGamma;

	scatteringProcess 
		<< a << b << c << q << dGamma << 2. << endr;
	scatteringProcesses = join_cols(scatteringProcesses, scatteringProcess);

	//
	// INTRAVALLEY ACOUSTIC TRANSVERSE PHONON EMISSION	
	//

	a = pow(hbar, 2) * (m_l + m_t + (-m_l + m_t) * cos(2. * theta_q)) / (4. * m_l * m_t); 
	b = hbar * v_st - pow(hbar, 2) * k_i * cos(theta_i) * cos(theta_q) / m_l - pow(hbar, 2) * k_i * cos(phi_i - phi_q) * sin(theta_i) * sin(theta_q) / m_t;
	c = 0;
	q = k_i*m_t*cos(theta_i)*cos(theta_q)*pow(m_t*pow(cos(theta_q),2) + m_l*pow(sin(theta_q),2),-1) + 
							pow(hbar,-1)*(-(m_l*m_t*v_st*pow(m_t*pow(cos(theta_q),2) + m_l*pow(sin(theta_q),2),-1)) + 
							(pow(m_t*pow(cos(theta_q),2) + m_l*pow(sin(theta_q),2),-1)*
							pow(-4*(hbar*m_l*pow(k_i,2)*pow(sin(theta_i),2) - hbar*m_l*pow(k_i,2)*pow(cos(phi_i),2)*pow(sin(theta_i),2) - 
							hbar*m_l*pow(k_i,2)*pow(sin(phi_i),2)*pow(sin(theta_i),2))*
							(hbar*m_t*pow(cos(theta_q),2) + hbar*m_l*pow(sin(theta_q),2)) + 
							pow(2*m_l*m_t*v_st - 2*hbar*k_i*m_t*cos(theta_i)*cos(theta_q) - 2*hbar*k_i*m_l*cos(phi_i - phi_q)*sin(theta_i)*sin(theta_q),2),0.5))/2.) + 
							k_i*m_l*cos(phi_i - phi_q)*pow(m_t*pow(cos(theta_q),2) + m_l*pow(sin(theta_q),2),-1)*sin(theta_i)*sin(theta_q);
	
	if (!(q > 0))
	{
		q = 0; //SET q TO ZERO IF NEGATIVE OR IMAGINARY
		a = 0;
		b = 0;
		c = 0;
	}
		

	dGamma = 1/.4 * sin(theta_q) * pow( Theta_u * sin(theta_q) * cos(theta_q) , 2) * (N_q + 1.) * pow(q, 3) / (8. * pow(M_PI,2) * rho * v_st * abs(2. * a * q + b));

	if (!(dGamma > 0))
	{
		dGamma = -12345; // nonsense value for non-physical scatter
	}
	
	dGammaTransverse = dGamma;

	scatteringProcess 
		<< a << b << c << q << dGamma << 3 << endr;
	scatteringProcesses = join_cols(scatteringProcesses, scatteringProcess);

///
/// ACOUSTIC TRANSVERSE INTERVALLEY EMISSIONS
///

	//
	// ACOUSTIC TRANSVERSE (SLOW) INTERVALLEY PHONON EMISSION ROOT 1
	//

	a = 0;
	b = 0;
	c = 0;
	q =	pow(hbar,-2)*pow(m_t*pow(cos(theta_q),2) + m_l*pow(sin(theta_q),2),-1)*
			(pow(pow(hbar,2)*(2.*m_l*m_t*(-hbar*w_slow)*(m_t*pow(cos(theta_q),2) + m_l*pow(sin(theta_q),2)) + 
			pow(hbar,2)*pow(k_i,2)*pow(m_t*cos(theta_i)*cos(theta_q) + m_l*cos(phi_i - phi_q)*sin(theta_i)*sin(theta_q),2)),0.5) + 
			k_i*pow(hbar,2)*(m_t*cos(theta_i)*cos(theta_q) + m_l*cos(phi_i - phi_q)*sin(theta_i)*sin(theta_q)));

	if (q > 0)
	{
		dGamma = 1/.4 * (7./2.) * (pow(3.20436e-10,2)*pow(hbar,-2)*pow(M_PI,-2)*pow(q,2)*pow(rho,-1)*pow(w_slow,-1)*
					pow(pow(cos(theta_q)*(-(k_i*cos(theta_i)) + q*cos(theta_q))*pow(m_l,-1) + 
					pow(m_t,-1)*sin(theta_q)*(-(k_i*cos(phi_i - phi_q)*sin(theta_i)) + q*sin(theta_q)),2),-0.5)*sin(theta_q))/8.;
							
		if (!(dGamma > 0))
		{
			dGamma = -12345; // nonsense value for non-physical scatter
		}

		dGammaTranInterSlow1 = dGamma;

		scatteringProcess 
			<< a << b << c << q << dGamma << 4 << endr;
		scatteringProcesses = join_cols(scatteringProcesses, scatteringProcess);
	}

	//
	// ACOUSTIC TRANSVERSE (SLOW) INTERVALLEY PHONON EMISSION ROOT 2
	//

	a = 0;
	b = 0;
	c = 0;
	q =	pow(hbar,-2)*pow(m_t*pow(cos(theta_q),2) + m_l*pow(sin(theta_q),2),-1)*
			(-pow(pow(hbar,2)*(2*m_l*m_t*(-hbar*w_slow)*(m_t*pow(cos(theta_q),2) + m_l*pow(sin(theta_q),2)) + 
			pow(hbar,2)*pow(k_i,2)*pow(m_t*cos(theta_i)*cos(theta_q) + m_l*cos(phi_i - phi_q)*sin(theta_i)*sin(theta_q),2)),0.5) + 
			k_i*pow(hbar,2)*(m_t*cos(theta_i)*cos(theta_q) + m_l*cos(phi_i - phi_q)*sin(theta_i)*sin(theta_q)));						

	if (q > 0)
	{
		dGamma = 1/.4 * (7./2.) * (pow(3.20436e-10,2)*pow(hbar,-2)*pow(M_PI,-2)*pow(q,2)*pow(rho,-1)*pow(w_slow,-1)*
				pow(pow(cos(theta_q)*(-(k_i*cos(theta_i)) + q*cos(theta_q))*pow(m_l,-1) + 
				pow(m_t,-1)*sin(theta_q)*(-(k_i*cos(phi_i - phi_q)*sin(theta_i)) + q*sin(theta_q)),2),-0.5)*sin(theta_q))/8.;
		
		if (!(dGamma > 0))
		{
			dGamma = -12345; // nonsense value for non-physical scatter
		}

		dGammaTranInterSlow2 = dGamma;

		scatteringProcess 
			<< a << b << c << q << dGamma << 5 << endr;
		scatteringProcesses = join_cols(scatteringProcesses, scatteringProcess);
	}

///
/// OPTICAL INTERVALLEY EMISSIONS
///

	//
	//  INTERVALLEY OPTICAL PHONON EMISSION ROOT 1
	//

	a = 0;
	b = 0;
	c = 0;
	q = pow(hbar,-2)*pow(m_t*pow(cos(theta_q),2) + m_l*pow(sin(theta_q),2),-1)*
			(pow(pow(hbar,2)*(2.*m_l*m_t*(-hw_intervalley)*(m_t*pow(cos(theta_q),2) + m_l*pow(sin(theta_q),2)) + 
			pow(hbar,2)*pow(k_i,2)*pow(m_t*cos(theta_i)*cos(theta_q) + m_l*cos(phi_i - phi_q)*sin(theta_i)*sin(theta_q),2)),0.5) + 
			k_i*pow(hbar,2)*(m_t*cos(theta_i)*cos(theta_q) + m_l*cos(phi_i - phi_q)*sin(theta_i)*sin(theta_q)));

	if (q > 0)
	{
		dGamma = 1/.4 * (7./2.) * (pow(4.80653199e-9,2)*pow(hbar,-2)*pow(M_PI,-2)*pow(q,2)*pow(rho,-1)*pow(w_intervalley,-1)*
					pow(pow(cos(theta_q)*(-(k_i*cos(theta_i)) + q*cos(theta_q))*pow(m_l,-1) + 
					pow(m_t,-1)*sin(theta_q)*(-(k_i*cos(phi_i - phi_q)*sin(theta_i)) + q*sin(theta_q)),2),-0.5)*sin(theta_q))/8.;
								
		if (!(dGamma > 0))
		{
			dGamma = -12345; // nonsense value for non-physical scatter
		}
	
		dGammaIntervalley1 = dGamma;

		scatteringProcess 
			<< a << b << c << q << dGamma << 6 << endr;
		scatteringProcesses = join_cols(scatteringProcesses, scatteringProcess);
	}

	//
	// INTERVALLEY OPTICAL PHONON EMISSION ROOT 2
	//

	a = 0;
	b = 0;
	c = 0;
	q = pow(hbar,-2)*pow(m_t*pow(cos(theta_q),2) + m_l*pow(sin(theta_q),2),-1)*
			(-pow(pow(hbar,2)*(2*m_l*m_t*(-hw_intervalley)*(m_t*pow(cos(theta_q),2) + m_l*pow(sin(theta_q),2)) + 
			pow(hbar,2)*pow(k_i,2)*pow(m_t*cos(theta_i)*cos(theta_q) + m_l*cos(phi_i - phi_q)*sin(theta_i)*sin(theta_q),2)),0.5) + 
			k_i*pow(hbar,2)*(m_t*cos(theta_i)*cos(theta_q) + m_l*cos(phi_i - phi_q)*sin(theta_i)*sin(theta_q)));
	
	if (q > 0)
	{

		dGamma = 1/.4 * (7./2.) * (pow(4.80653199e-9,2)*pow(hbar,-2)*pow(M_PI,-2)*pow(q,2)*pow(rho,-1)*pow(w_intervalley,-1)*
						pow(pow(cos(theta_q)*(-(k_i*cos(theta_i)) + q*cos(theta_q))*pow(m_l,-1) + 
						pow(m_t,-1)*sin(theta_q)*(-(k_i*cos(phi_i - phi_q)*sin(theta_i)) + q*sin(theta_q)),2),-0.5)*sin(theta_q))/8.;
							
		if (!(dGamma > 0))
		{
			dGamma = -12345; // nonsense value for non-physical scatter
		}
	
		dGammaIntervalley2 = dGamma;

		scatteringProcess 
			<< a << b << c << q << dGamma << 7 << endr;
		scatteringProcesses = join_cols(scatteringProcesses, scatteringProcess);
	
	}

///
/// INTRAVALLEY OPTICAL EMISSIONS
///

	//
	// INTRAVALLEY OPTICAL PHONON EMISSION ROOT 1
	//

	a = pow(hbar, 2) * (m_l + m_t + (-m_l + m_t) * cos(2 * theta_q)) / (4. * m_l * m_t); 
	b = -pow(hbar, 2) * k_i * ((m_t * cos(theta_i) * cos(theta_q) + m_l * cos(phi_i - phi_q) * sin(theta_i) * sin(theta_q))) / (m_l * m_t);
	c = hw_optical;
	q = pow(hbar,-2)*pow(m_t*pow(cos(theta_q),2) + m_l*pow(sin(theta_q),2),-1)*
			(pow(pow(hbar,2)*(2*m_l*m_t*(-1 * hw_optical)*(m_t*pow(cos(theta_q),2) + m_l*pow(sin(theta_q),2)) + 
			pow(hbar,2)*pow(k_i,2)*pow(m_t*cos(theta_i)*cos(theta_q) + m_l*cos(phi_i - phi_q)*sin(theta_i)*sin(theta_q),2)),0.5) + 
			k_i*pow(hbar,2)*(m_t*cos(theta_i)*cos(theta_q) + m_l*cos(phi_i - phi_q)*sin(theta_i)*sin(theta_q)));

	if (q > 0)
	{

		dGamma = 1/.4 * sin(theta_q) * hbar * pow(Theta_o, 2) * (N_q + 1.) * pow(q, 2) / (8. * pow(M_PI,2) * rho * hw_optical * abs(2. * a * q + b));

		if (!(dGamma > 0))
		{
			dGamma = -12345; // nonsense value for non-physical scatter
		}

		dGammaOptical1 = dGamma;

		scatteringProcess 
			<< a << b << c << q << dGamma << 8 << endr;
		scatteringProcesses = join_cols(scatteringProcesses, scatteringProcess);
	}

	//
	// INTRAVALLEY OPTICAL PHONON EMISSION ROOT 2
	//

	a = pow(hbar, 2) * (m_l + m_t + (-m_l + m_t) * cos(2 * theta_q)) / (4. * m_l * m_t); 
	b = -pow(hbar, 2) * k_i * ((m_t * cos(theta_i) * cos(theta_q) + m_l * cos(phi_i - phi_q) * sin(theta_i) * sin(theta_q))) / (m_l * m_t);
	c = hw_optical;
	q = pow(hbar,-2)*pow(m_t*pow(cos(theta_q),2) + m_l*pow(sin(theta_q),2),-1)*
			(-pow(pow(hbar,2)*(2*m_l*m_t*(-1 * hw_optical)*(m_t*pow(cos(theta_q),2) + m_l*pow(sin(theta_q),2)) + 
			pow(hbar,2)*pow(k_i,2)*pow(m_t*cos(theta_i)*cos(theta_q) + m_l*cos(phi_i - phi_q)*sin(theta_i)*sin(theta_q),2)),0.5) + 
			k_i*pow(hbar,2)*(m_t*cos(theta_i)*cos(theta_q) + m_l*cos(phi_i - phi_q)*sin(theta_i)*sin(theta_q)));

	if (q > 0)
	{
		dGamma = 1/.4 * sin(theta_q) * hbar * pow(Theta_o, 2) * (N_q + 1.) * pow(q, 2) / (8. * pow(M_PI,2) * rho * hw_optical * abs(2. * a * q + b));
						
		if (!(dGamma > 0))
		{
			dGamma = -12345; // nonsense value for non-physical scatter
		}

		dGammaOptical2 = dGamma;

		scatteringProcess 
			<< a << b << c << q << dGamma << 9 << endr;
		scatteringProcesses = join_cols(scatteringProcesses, scatteringProcess);

	}


	//SELF SCATTER

	double selfScatterRate = Gamma_0 / (4 * M_PI); 
	for (int i = 0; i < scatteringProcesses.n_rows; i++)
		if (!(scatteringProcesses.at(i, 4) == -12345))
			selfScatterRate -= scatteringProcesses.at(i, 4);


	scatteringProcess
		<< 0 << 0 << 0 << 0 << selfScatterRate << 1 << endr;
	scatteringProcesses = join_cols(scatteringProcess, scatteringProcesses);
		

	//Create reduced matrix that only includes physical scatters
	arma::mat ReducedGammaMatrix;
	int n = scatteringProcesses.n_rows;
	for (int i = 0; i < n; i++)
	{
		if (!(scatteringProcesses.at(i, 4) == -12345))
		{
			ReducedGammaMatrix = join_cols(ReducedGammaMatrix, scatteringProcesses.row(i));
		}
	}

	if (selfScatterRate > 0)
	{

		//PICK RANDOM SCATTERING PROCESS
		double random = GetRandomNumber();
		int n_reduced = ReducedGammaMatrix.n_rows;
		int width = scatteringProcesses.n_cols;
		int selectedScatteringProcess = 0;


		//ITERATE OVER ALL SCATERRING PROCESSES TO PICK A PROCESS

		LessThan = 0;
		GreaterThan = 0;

		for (int i = 0; i < n_reduced; i++)
		{
			sumLessThan = 0;
			sumGreaterThan = 0;

			for (int j = 0; j < i;     j++)
				sumLessThan += ReducedGammaMatrix.at(j+1, 4);
			for (int j = 0; j < i + 1; j++)
				if (j == n_reduced - 1)
				{
					sumGreaterThan += ReducedGammaMatrix.at(0,4);
				}
				else
					sumGreaterThan += ReducedGammaMatrix.at(j+1, 4);

			LessThan = sumLessThan / Gamma_0;
			GreaterThan = sumGreaterThan / Gamma_0;

			if (LessThan <= random && random < GreaterThan)
			{
				selectedScatteringProcess = i;
				break;
			}
		}
	
		arma::vec q_vector_spherical;
		q_vector_spherical
			<< ReducedGammaMatrix.at(selectedScatteringProcess, 3) << theta_q << phi_q;


		double q_magnitude = q_vector_spherical.at(0);
		if (q_magnitude != 0) //If the scatter process is a physical scatter
			numberOfScatters = numberOfScatters + 1;

		arma::vec q_vector_cartesian = q_vector_spherical;
		ConvertSphericalToCartesian(&q_vector_cartesian);

		if (ReducedGammaMatrix.at(selectedScatteringProcess, 5) == 4 or ReducedGammaMatrix.at(selectedScatteringProcess, 5) == 5 or ReducedGammaMatrix.at(selectedScatteringProcess, 5) == 6 or ReducedGammaMatrix.at(selectedScatteringProcess, 5) == 7)
		{
			Random_old_Valley_Number = valley_number;
			valley_number = GetRandomValley();
			Random_Valley_Number = valley_number;
		} 

		electron->AfterScatterkVector = crystal->getRotationMatrix(valley_number + 4) * (k_initial_vector_cartesian - q_vector_cartesian);

		mat inverseMassTensor = crystal->getInverseMassTensor(valley_number);
		electron->AfterScattergroupVelocityVector = hbar * inverseMassTensor * electron->AfterScatterkVector;	
		electron->AfterScatterkineticEnergy = (double)(pow(hbar,2) / 2) * as_scalar(electron->AfterScatterkVector.t() * inverseMassTensor * electron->AfterScatterkVector);
		
		arma::vec q_vector_cartesian_labframe = crystal->getRotationMatrix(valley_number + 4) * q_vector_cartesian;
		arma::vec q_vector_spherical_lab = q_vector_cartesian_labframe;
		ConvertCartesianToSpherical(&q_vector_spherical_lab);
		arma::vec q_vector_cartesian_valley = crystal->getRotationMatrix(valley_number) * q_vector_cartesian_labframe;
		arma::vec q_vector_spherical_valley = q_vector_cartesian_valley;
		ConvertCartesianToSpherical(&q_vector_spherical_valley);


		//if (stepNumber == 10000)
		//{
				//TAKE DATA
				//dataFileStream << stepNumber << "\n" << endl;
				//dataFileStream << "\n" << stepNumber << "\n" << endl;
				//dataFileStream << "\n" << ReducedGammaMatrix.at(selectedScatteringProcess, 0) << "\n" << endl;
				//dataFileStream << "\n" << ReducedGammaMatrix.at(selectedScatteringProcess, 1) << "\n" << endl;
				//dataFileStream << "\n" << ReducedGammaMatrix.at(selectedScatteringProcess, 2) << "\n" << endl;
				//dataFileStream << "\n" << ReducedGammaMatrix.at(selectedScatteringProcess, 3) << "\n" << endl;
				//dataFileStream << ReducedGammaMatrix.at(selectedScatteringProcess, 5) << "\n" << endl;
				//dataFileStream << "\n" << q_vector_cartesian_labframe << endl;
				//dataFileStream << "\n" << q_vector_spherical_lab << endl;
				//dataFileStream << "\n" << q_vector_cartesian_valley << endl;	
				//dataFileStream << "\n" << q_vector_spherical_valley << endl;		
				//dataFileStream << "\n" << theta_q << "\n" << endl;
				//dataFileStream << "\n" << phi_q << "\n" << endl;
				//dataFileStream << "\n" << random << "\n" << endl;
				//dataFileStream << "\n" << Gamma_0 << "\n" << endl;
				//dataFileStream <<  selfScatterRate << "\n" << endl;
				//dataFileStream << "\n" << dGammaLongitudinal << "\n" << endl;
				//dataFileStream << "\n" << dGammaTransverse << "\n" << endl;
				//dataFileStream << "\n" << endl;
				//dataFileStream << "\n" << selectedScatteringProcess << endl;
				//dataFileStream << "\n" << globalTime << "\n" << endl;
				//dataFileStream << "\n" << electron->AfterEpositionVector / globalTime << endl;
		//}
				//QUESTION: WHAT IS THE TOTAL PHYSICAL PARTIAL RATE? - EVERYTHING BUT SELF SCATTER

		Converter(electron, crystal, valley_number);

		time_between_physical_scatters += timeStep;

			if (!(selectedScatteringProcess == 0)) 
			{
				physical_scatter_check += 1;
			}

			Phys_self_scatter_check += 1;
		
		//RecordData(electron, stepNumber, timeStep, valley_number);			
		Initializer(electron);			

		sumation = 0;
		ratio = 0;
		for (int i = 0; i < (n_reduced); i++)
		{
			sumation += ReducedGammaMatrix.at(i+1, 4);
		}
		ratio = 4 * M_PI * sumation / Gamma_0;
		if (ratio < 0.001)
			Gamma_0 = Gamma_0 * 0.9;
		if (ratio > 0.01)
			Gamma_0 =  Gamma_0 * 1.1;

	}
	else
	{
		globalTime -= timeStep;
		Gamma_0 = 40 * M_PI * Gamma_0 * 1.3;
	}
	if (stepNumber == 10000)
	{
		dataFileStream << "\n" << electron->AfterEpositionVector / globalTime << endl;
	}
}

//
//BOOK KEEPING FUNCTIONS
//

void Simulator::RecordData(Electron* electron, int stepNumber, double timeStep, int valley_number)
{
	dataFileStream << "\n" << timeStep << "\n" << endl;
	dataFileStream << "\n" << globalTime << "\n" << endl;
	dataFileStream << "\n" << valley_number << "\n" << endl;
	dataFileStream << "\n" << valley_number << "\n" << endl;

	dataFileStream << "\n" << electron->kVector << endl;
	dataFileStream << "\n" << electron->kVector_spherical_lab << endl;
	dataFileStream << "\n" << electron->kVector_cartesian_valley << endl;	
	dataFileStream << "\n" << electron->kVector_spherical_valley << endl;		

	dataFileStream << "\n" << electron->groupVelocityVector << endl;
	dataFileStream << "\n" << electron->groupVelocityVector_spherical_lab << endl;
	dataFileStream << "\n" << electron->groupVelocityVector_cartesian_valley << endl;	
	dataFileStream << "\n" << electron->groupVelocityVector_spherical_valley << endl;

	dataFileStream << "\n" << electron->positionVector << endl;
	dataFileStream << "\n" << electron->positionVector_spherical_lab << endl;
	dataFileStream << "\n" << electron->positionVector_cartesian_valley << endl;	
	dataFileStream << "\n" << electron->positionVector_spherical_valley << endl;

	dataFileStream << electron->kineticEnergy << "\n" << endl;

	dataFileStream << "\n" << electron->AfterEkVector << endl;
	dataFileStream << "\n" << electron->AfterEkVector_spherical_lab << endl;
	dataFileStream << "\n" << electron->AfterEkVector_cartesian_valley << endl;	
	dataFileStream << "\n" << electron->AfterEkVector_spherical_valley << endl;

	dataFileStream << "\n" << electron->AfterEgroupVelocityVector << endl;
	dataFileStream << "\n" << electron->AfterEgroupVelocityVector_spherical_lab << endl;
	dataFileStream << "\n" << electron->AfterEgroupVelocityVector_cartesian_valley << endl;	
	dataFileStream << "\n" << electron->AfterEgroupVelocityVector_spherical_valley << endl;	

	dataFileStream << "\n" << electron->AfterEpositionVector << endl;
	dataFileStream << "\n" << electron->AfterEpositionVector_spherical_lab << endl;
	dataFileStream << "\n" << electron->AfterEpositionVector_cartesian_valley << endl;	
	dataFileStream << "\n" << electron->AfterEpositionVector_spherical_valley << endl;

	dataFileStream << electron->AfterEkineticEnergy << "\n" << endl;

	dataFileStream << "\n" << electron->AfterScatterkVector << endl;
	dataFileStream << "\n" << electron->AfterScatterkVector_spherical_lab << endl;
	dataFileStream << "\n" << electron->AfterScatterkVector_cartesian_valley << endl;	
	dataFileStream << "\n" << electron->AfterScatterkVector_spherical_valley << endl;

	dataFileStream << "\n" << electron->AfterScattergroupVelocityVector << endl;
	dataFileStream << "\n" << electron->AfterScattergroupVelocityVector_spherical_lab << endl;
	dataFileStream << "\n" << electron->AfterScattergroupVelocityVector_cartesian_valley << endl;	
	dataFileStream << "\n" << electron->AfterScattergroupVelocityVector_spherical_valley << endl;

	dataFileStream << "\n" << electron->AfterEpositionVector << endl;
	dataFileStream << "\n" << electron->AfterEpositionVector_spherical_lab << endl;
	dataFileStream << "\n" << electron->AfterEpositionVector_cartesian_valley << endl;	
	dataFileStream << "\n" << electron->AfterEpositionVector_spherical_valley << endl;

	dataFileStream << electron->AfterScatterkineticEnergy << "\n" << endl;
}

void Simulator::Initializer(Electron* electron)
{
	electron->kVector = electron->AfterScatterkVector;
	electron->groupVelocityVector = electron->AfterScattergroupVelocityVector;
	electron->positionVector = electron->AfterEpositionVector;
	electron->kineticEnergy = electron->AfterScatterkineticEnergy;
}

void Simulator::Converter(Electron* electron, Crystal* crystal, int valley_number)
{
	electron->kVector_spherical_lab = electron->kVector;
	ConvertCartesianToSpherical(&electron->kVector_spherical_lab);
	electron->kVector_cartesian_valley = crystal->getRotationMatrix(valley_number) * electron->kVector; 
	//Note The BeforEfield vectors in the valley frame may be different from the previous after scatter vectors due to the potentiality that the random valley number is different between the two iterations. 
	electron->kVector_spherical_valley = electron->kVector_cartesian_valley;
	ConvertCartesianToSpherical(&electron->kVector_spherical_valley);

	electron->groupVelocityVector_spherical_lab = electron->groupVelocityVector;
	ConvertCartesianToSpherical(&electron->groupVelocityVector_spherical_lab);
	electron->groupVelocityVector_cartesian_valley = crystal->getRotationMatrix(valley_number) * electron->groupVelocityVector; 
	electron->groupVelocityVector_spherical_valley = electron->groupVelocityVector_cartesian_valley;
	ConvertCartesianToSpherical(&electron->groupVelocityVector_spherical_valley);

	electron->positionVector_spherical_lab = electron->positionVector;
	ConvertCartesianToSpherical(&electron->positionVector_spherical_lab);
	electron->positionVector_cartesian_valley = crystal->getRotationMatrix(valley_number) * electron->positionVector; 
	electron->positionVector_spherical_valley = electron->positionVector_cartesian_valley;
	ConvertCartesianToSpherical(&electron->positionVector_spherical_valley);

	electron->AfterEkVector_spherical_lab = electron->AfterEkVector;
	ConvertCartesianToSpherical(&electron->AfterEkVector_spherical_lab);
	electron->AfterEkVector_cartesian_valley = crystal->getRotationMatrix(valley_number) * electron->AfterEkVector; 
	electron->AfterEkVector_spherical_valley = electron->AfterEkVector_cartesian_valley;
	ConvertCartesianToSpherical(&electron->AfterEkVector_spherical_valley);

	electron->AfterEgroupVelocityVector_spherical_lab = electron->AfterEgroupVelocityVector;
	ConvertCartesianToSpherical(&electron->AfterEgroupVelocityVector_spherical_lab);
	electron->AfterEgroupVelocityVector_cartesian_valley = crystal->getRotationMatrix(valley_number) * electron->AfterEgroupVelocityVector; 
	electron->AfterEgroupVelocityVector_spherical_valley = electron->AfterEgroupVelocityVector_cartesian_valley;
	ConvertCartesianToSpherical(&electron->AfterEgroupVelocityVector_spherical_valley);

	electron->AfterEpositionVector_spherical_lab = electron->AfterEpositionVector;
	ConvertCartesianToSpherical(&electron->AfterEpositionVector_spherical_lab);
	electron->AfterEpositionVector_cartesian_valley = crystal->getRotationMatrix(valley_number) * electron->AfterEpositionVector; 
	electron->AfterEpositionVector_spherical_valley = electron->AfterEpositionVector_cartesian_valley;
	ConvertCartesianToSpherical(&electron->AfterEpositionVector_spherical_valley);

	electron->AfterScatterkVector_spherical_lab = electron->AfterScatterkVector;
	ConvertCartesianToSpherical(&electron->AfterScatterkVector_spherical_lab);
	electron->AfterScatterkVector_cartesian_valley = crystal->getRotationMatrix(valley_number) * electron->AfterScatterkVector; 
	electron->AfterScatterkVector_spherical_valley = electron->AfterScatterkVector_cartesian_valley;
	ConvertCartesianToSpherical(&electron->AfterScatterkVector_spherical_valley);

	electron->AfterScattergroupVelocityVector_spherical_lab = electron->AfterScattergroupVelocityVector;
	ConvertCartesianToSpherical(&electron->AfterScattergroupVelocityVector_spherical_lab);
	electron->AfterScattergroupVelocityVector_cartesian_valley = crystal->getRotationMatrix(valley_number) * electron->AfterScattergroupVelocityVector; 
	electron->AfterScattergroupVelocityVector_spherical_valley = electron->AfterScattergroupVelocityVector_cartesian_valley;
	ConvertCartesianToSpherical(&electron->AfterScattergroupVelocityVector_spherical_valley);
}

const std::string Simulator::GetCurrentDateTime()
{
	time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://www.cplusplus.com/reference/clibrary/ctime/strftime/
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);

    return buf;
}


//
//MERSENNE TWISTER RANDOM NUMBER GENERATOR
//RETURNS RANDOM NUMBER BETWEEN 0 and 1
//
double Simulator::GetRandomNumber()
{
	return u01(gen);
}

int Simulator::GetRandomValley()
{
	double test = GetRandomNumber() * 4;
	int random_valley_integer = 0;
	for (int i = 0; i < 4; i++)
	{
		if (i <= test && test <= (i + 1))
		{
			random_valley_integer = i;
			break;
		}
	}
	return(random_valley_integer);
}