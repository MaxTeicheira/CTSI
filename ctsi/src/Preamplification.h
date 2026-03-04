#ifndef PREAMP_H
#define PREAMP_H

#include <iostream>
#include <fstream>
#include <vector>

using std::vector;

class Preamplification{
	public:
		Preamplification();
		~Preamplification();
		double			getSamplingInterval();
		int				importPreamp(void* ctsiObjectPtr, int indent, void (*dispFuncPtr)(void*, int), char *preampFile);
		vector<double>	preamplify(vector<double> &inSig);

	protected:
		double			samplingInterval;
		vector<double>	ipCoefficients;
		vector<double>	opCoefficients;
};

#endif // PREAMP_H
