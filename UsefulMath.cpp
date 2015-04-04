#include "UsefulMath.h"

//(x, y, z) -> (r, theta, phi), theta is polar, phi is azimuthal
void ConvertCartesianToSpherical(arma::vec* cartesianVectorPointer)
{
	double x = cartesianVectorPointer->at(0);
	double y = cartesianVectorPointer->at(1);
	double z = cartesianVectorPointer->at(2);

	double r = sqrt(x*x + y*y + z*z);

	if (r > 0)
	{
		cartesianVectorPointer->at(0) = r;
		cartesianVectorPointer->at(1) = acos(z/r);
		cartesianVectorPointer->at(2) = atan2(y,x); //atan2(y,x) gets arc tangent taking quadrant into account	
	}
	else
	{
		cartesianVectorPointer->at(0) = 0;
		cartesianVectorPointer->at(1) = 0;
		cartesianVectorPointer->at(2) = 0;	
	}
}

//(r, theta, phi) -> (x, y, z), theta is polar, phi is azimuthal
void ConvertSphericalToCartesian(arma::vec* sphericalVectorPointer)
{
	double r = sphericalVectorPointer->at(0);
	double theta = sphericalVectorPointer->at(1);
	double phi = sphericalVectorPointer->at(2);

	sphericalVectorPointer->at(0) = r * sin(theta) * cos(phi);
	sphericalVectorPointer->at(1) = r * sin(theta) * sin(phi);
	sphericalVectorPointer->at(2) = r * cos(theta);
}