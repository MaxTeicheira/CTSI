#ifndef MYRANDOM_H
#define MYRANDOM_H

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

class MyRandom
{
	public:
		MyRandom();
		~MyRandom();
		double		uniform();
		double		normal();
		unsigned	randSeed();
		void		setSeed(int seed);

	protected:
		double	myConst;
};

#endif // MYRANDOM_H
