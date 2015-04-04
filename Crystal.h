#pragma once

#include "armadillo"

#include "Constants.h"

using namespace arma;
class Crystal
{
public:
	enum ValleyType {
		VALLEYTYPE_LVALLEY_1,
        VALLEYTYPE_LVALLEY_2,
	    VALLEYTYPE_LVALLEY_3,
	    VALLEYTYPE_LVALLEY_4,
	    VALLEYTYPE_GAMMAVALLEY,
	    VALLEYTYPE_XVALLEY_PLUS_X,
	    VALLEYTYPE_XVALEY_MINUS_X,
    	VALLEYTYPE_XVALEY_PLUS_Y,
		VALLEYTYPE_XVALEY_MINUS_Y,
		VALLEYTYPE_XVALEY_PLUS_Z,
		VALLEYTYPE_XVALEY_MINUS_Z,
		VALLEYTYPE_LVALLEY_INVERSE_1,
		VALLEYTYPE_LVALLEY_INVERSE_2,
		VALLEYTYPE_LVALLEY_INVERSE_3,
		VALLEYTYPE_LVALLEY_INVERSE_4,
   };

	Crystal(void);
	~Crystal(void);

	mat getInverseMassTensor(int valleyNumber);
	mat getRotationMatrix(int valleyNumber);
	
};

