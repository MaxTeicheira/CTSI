//===================================================================
// Class authored by Yi Gu on Jan 27, 2011 based on implementations
// presented on http://eternallyconfuzzled.com/arts/jsw_art_rand.aspx
// and http://www.taygeta.com/random/gaussian.html.
//===================================================================
#include "MyRandom.h"

MyRandom::MyRandom()
{
	myConst = 1/(((double) RAND_MAX) + 1);
}

MyRandom::~MyRandom()
{
}

//===================================================================
// This function returns a uniformly distributed pseudo-random
// variable in the range [0, 1).
// Implementation is drawn from 
// http://eternallyconfuzzled.com/arts/jsw_art_rand.aspx
double MyRandom::uniform()
{
	return rand()*myConst;
}

//===================================================================
// This function returns a normally distributed pseudo-random
// variable with mean 0 and variance 1. Implementation is drawn from
// http://www.taygeta.com/random/gaussian.html.
double MyRandom::normal()
{
	double x1 = 0, x2 = 0, w = 0;

	do
	{
		x1 = 2*uniform() - 1;
		x2 = 2*uniform() - 1;
		w = x1*x1 + x2*x2;
	}
	while (w >= 1);

	w = sqrt((-2*log(w))/w);
	return x1*w;
}

//===================================================================
// This function returns a random seed value.
// Implementation is drawn from 
// http://eternallyconfuzzled.com/arts/jsw_art_rand.aspx
unsigned MyRandom::randSeed()
{
	time_t now = time(0);
	unsigned char *p = (unsigned char *)&now;
	unsigned seed = 0;
	size_t i;
	
	for (i = 0; i < sizeof now; i++)
		seed = seed*(UCHAR_MAX + 2U) + p[i];
    
	return seed;
}

//===================================================================
// This function sets the seed for rand()
void MyRandom::setSeed(int seed)
{
	srand(seed);
}
