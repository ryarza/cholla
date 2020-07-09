#ifdef STARS

#ifndef STARS_H
#define STARS_H

#include<stdio.h>
#include <cmath>
#include"../global.h"

class Stellar
{
public:
	Real Mstar;
	Real Rstar;
	Real perturberMass;
	Real polyN;
	Real periBeta;

	Real posBhx;
	Real posBhy;
	Real posBhz;

	Real posStx;
	Real posSty;
	Real posStz;

	Real velBhx;
	Real velBhy;
	Real velBhz;

	Real velStx;
	Real velSty;
	Real velStz;

	Real accBhx;
	Real accBhy;
	Real accBhz;

	Real accStx;
	Real accSty;
	Real accStz;
};

#endif
#endif
