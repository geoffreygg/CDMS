#include "Crystal.h"

Crystal::Crystal(void)
{
}

Crystal::~Crystal(void)
{
}

mat Crystal::getInverseMassTensor(int valleyNumber)
{
	mat inverseMassTensor(3,3);
	inverseMassTensor.fill(0);
		
	switch(valleyNumber)
	{
		case VALLEYTYPE_LVALLEY_1:
			inverseMassTensor	<< (1./m_l + 2./m_t)/3.	<< (1./m_l - 1./m_t)/3.	<< (1./m_l - 1./m_t)/3.		<< endr
								<< (1./m_l - 1./m_t)/3.		<< (1./m_l + 2./m_t)/3.		<< (1./m_l - 1./m_t)/3.		<< endr
								<< (1./m_l - 1./m_t)/3.		<< (1./m_l - 1./m_t)/3.		<< (1./m_l + 2./m_t)/3.;
			return inverseMassTensor;
			break;

		case VALLEYTYPE_LVALLEY_2: 
			inverseMassTensor	<< (1./m_l + 2./m_t)/3.		<< (-(1./m_l) + 1./m_t)/3.	<< (-(1./m_l) + 1./m_t)/3.	<< endr
								<< (-(1./m_l) + 1./m_t)/3.   << (1./m_l + 2./m_t)/3.		<< (1./m_l - 1./m_t)/3.		<< endr
								<< (-(1./m_l) + 1./m_t)/3.	<< (1./m_l - 1./m_t)/3.		<< (1./m_l + 2./m_t)/3.;
			return inverseMassTensor;
			break;

		case VALLEYTYPE_LVALLEY_3: 
			inverseMassTensor	<< (1./m_l + 2./m_t)/3.		<< (1./m_l - 1./m_t)/3.		<< (-(1./m_l) + 1./m_t)/3.	<< endr
								<< (1./m_l - 1./m_t)/3.		<< (1./m_l + 2./m_t)/3.		<< (-(1./m_l) + 1./m_t)/3.	<< endr
								<< (-(1./m_l) + 1./m_t)/3.	<< (-(1./m_l) + 1./m_t)/3.	<< (1./m_l + 2./m_t)/3.;
			return inverseMassTensor;
			break;

		case VALLEYTYPE_LVALLEY_4: 
			inverseMassTensor	<< (1./m_l + 2./m_t)/3.		<< (-(1./m_l) + 1./m_t)/3.	<< (1./m_l - 1./m_t)/3.		<< endr
								<< (-(1./m_l) + 1./m_t)/3.	<< (1./m_l + 2./m_t)/3.		<< (-(1./m_l) + 1./m_t)/3.	<< endr
								<< (1./m_l - 1./m_t)/3.		<< (-(1./m_l) + 1./m_t)/3.	<< (1./m_l + 2./m_t)/3.;
			return inverseMassTensor;
			break;

		case VALLEYTYPE_GAMMAVALLEY:  // Gamma valley 
			inverseMassTensor	<< 1./m_G		<< 0		<< 0		<< endr
								<< 0		<< 1./m_G		<< 0		<< endr
								<< 0		<< 0		<< 1./m_G;
			return inverseMassTensor;
			break;

		case VALLEYTYPE_XVALLEY_PLUS_X:  // X-valley: +x 
			inverseMassTensor	<< 1./m_lX	<< 0		<< 0		<< endr
								<< 0		<< 1./m_tX	<< 0		<< endr
								<< 0		<< 0		<< 1./m_tX;
			return inverseMassTensor;
			break;				

		case VALLEYTYPE_XVALEY_MINUS_X: // X-valley: -x 
			inverseMassTensor	<< 1./m_lX	<< 0		<< 0		<< endr
								<< 0		<< 1./m_tX	<< 0		<< endr
								<< 0		<< 0		<< 1./m_tX;
			return inverseMassTensor;
			break;			

		case VALLEYTYPE_XVALEY_PLUS_Y: // X-valley: +y  
			inverseMassTensor	<< 1./m_tX	<< 0		<< 0		<< endr
								<< 0		<< 1./m_lX	<< 0		<< endr
								<< 0		<< 0		<< 1./m_tX;	
			return inverseMassTensor;
			break;

		case VALLEYTYPE_XVALEY_MINUS_Y: // X-valley: -y  
			inverseMassTensor	<< 1./m_tX	<< 0		<< 0		<< endr
								<< 0		<< 1./m_lX	<< 0		<< endr
								<< 0		<< 0		<< 1./m_tX;	
			return inverseMassTensor;
			break;			

		case VALLEYTYPE_XVALEY_PLUS_Z: // X-valley: +z  
			inverseMassTensor	<< 1./m_tX	<< 0		<< 0		<< endr
								<< 0		<< 1./m_tX	<< 0		<< endr
								<< 0		<< 0		<< 1./m_lX;	
			return inverseMassTensor;
			break;

		case VALLEYTYPE_XVALEY_MINUS_Z: //X-valley: -z
			inverseMassTensor	<< 1./m_tX	<< 0		<< 0		<< endr
								<< 0		<< 1./m_tX	<< 0		<< endr
								<< 0		<< 0		<< 1./m_lX;
			return inverseMassTensor;
			break;

		
			

		//NO APPROPRIATE VALLEY WAS IDENTIFIED
		default:
			mat matZero(3,3); 
			matZero.fill(0);
			return matZero;
			break;

	}
}

mat Crystal::getRotationMatrix(int valleyNumber)
{
	mat rotationMatrix(3,3);
	rotationMatrix.fill(0);
		
	switch(valleyNumber)
	{
		case VALLEYTYPE_LVALLEY_1: 
			rotationMatrix		<< (3. + sqrt(3.)) / 6.	<< (-3. + sqrt(3.)) / 6.	<< 2. * sqrt(3.) / 6.		<< endr
								<< (-3. + sqrt(3.)) / 6.	<< (3. + sqrt(3.)) / 6.	<< 2. * sqrt(3.) / 6.		<< endr
								<< -1. / sqrt(3.)			<< -1. / sqrt(3.)			<< 1. / sqrt(3.);
			return rotationMatrix;
			break;

		case VALLEYTYPE_LVALLEY_2: 
			rotationMatrix		<< (3. + sqrt(3.)) / 6.	<< -(-3. + sqrt(3.)) / 6.	<< -2. * sqrt(3.) / 6.		<< endr
								<< -(-3. + sqrt(3.)) / 6.	<< (3. + sqrt(3.)) / 6.	<< 2. * sqrt(3.) / 6.		<< endr
								<< 1. / sqrt(3.)			<< -1. / sqrt(3.)			<< 1. / sqrt(3.);
			return rotationMatrix;
			break;

		case VALLEYTYPE_LVALLEY_3:
			rotationMatrix		<< (3. + sqrt(3.)) / 6.	<< (-3. + sqrt(3.)) / 6.	<< -2. * sqrt(3.) / 6.		<< endr
								<< (-3. + sqrt(3.)) / 6.	<< (3. + sqrt(3.)) / 6.	<< -2. * sqrt(3.) / 6.		<< endr
								<< 1. / sqrt(3.)			<< 1. / sqrt(3.)			<< 1. / sqrt(3.);
			return rotationMatrix;
			break;

		case VALLEYTYPE_LVALLEY_4: 
			rotationMatrix		<< (3. + sqrt(3.)) / 6.	<< -(-3. + sqrt(3.)) / 6.  << 2. * sqrt(3.) / 6.		<< endr
								<< -(-3. + sqrt(3.)) / 6.	<< (3. + sqrt(3.)) / 6.	<< -2. * sqrt(3.) / 6.		<< endr
								<< -1. / sqrt(3.)			<< 1. / sqrt(3.)			<< 1. / sqrt(3.);
			return rotationMatrix;
			break;

		case VALLEYTYPE_LVALLEY_INVERSE_1:
			rotationMatrix		<< (3. + sqrt(3.)) / 6.	<< (-3. + sqrt(3.)) / 6.	<< -2. * sqrt(3.) / 6.		<< endr
								<< (-3. + sqrt(3.)) / 6.	<< (3. + sqrt(3.)) / 6.	<< -2. * sqrt(3.) / 6.		<< endr
								<< sqrt(3.)			<< sqrt(3.)			<< sqrt(3.);
			return rotationMatrix;
			break;

		case VALLEYTYPE_LVALLEY_INVERSE_2:
			rotationMatrix		<< (3. + sqrt(3.)) / 6.	<< -(-3. + sqrt(3.)) / 6.	<< 2. * sqrt(3.) / 6.		<< endr
								<< -(-3. + sqrt(3.)) / 6.	<< (3. + sqrt(3.)) / 6.	<< -2. * sqrt(3.) / 6.		<< endr
								<< -sqrt(3.)			<< sqrt(3.)			<< sqrt(3.);
			return rotationMatrix;
			break;

		case VALLEYTYPE_LVALLEY_INVERSE_3:
			rotationMatrix		<< (3. + sqrt(3.)) / 6.	<< (-3. + sqrt(3.)) / 6.	<< 2. * sqrt(3.) / 6.		<< endr
								<< (-3. + sqrt(3.)) / 6.	<< (3. + sqrt(3.)) / 6.	<< 2. * sqrt(3.) / 6.		<< endr
								<< -sqrt(3.)			<< -sqrt(3.)			<< sqrt(3.);
			return rotationMatrix;
			break;

		case VALLEYTYPE_LVALLEY_INVERSE_4:
			rotationMatrix		<< (3. + sqrt(3.)) / 6.	<< -(-3. + sqrt(3.)) / 6.	<< -2. * sqrt(3.) / 6.		<< endr
								<< -(-3. + sqrt(3.)) / 6.	<< (3. + sqrt(3.)) / 6.	<< 2. * sqrt(3.) / 6.		<< endr
								<< sqrt(3.)			<< -sqrt(3.)			<< sqrt(3.);
			return rotationMatrix;
			break;

		
			

		//NO APPROPRIATE VALLEY WAS IDENTIFIED
		default:
			mat matZero(3,3); 
			matZero.fill(0);
			return matZero;
			break;

	}
}