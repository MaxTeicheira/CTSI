#include "Preamplification.h"

using namespace std;

Preamplification::Preamplification()
{
	samplingInterval = 0;
}

Preamplification::~Preamplification()
{
}

//===================================================================
// This function returns the sampling interval.
double Preamplification::getSamplingInterval()
{
	return samplingInterval;
}

//===================================================================
// This function reads in the coefficients of the discrete filter
// that models the preamplifier. The filter is implemented as a
// difference equation.
int Preamplification::importPreamp(void* ctsiObjectPtr, int indent, void (*dispFuncPtr)(void*, int), char *preampFile)
{
	char caption[256];
	double numNumCoeffs = 0, numDenCoeffs = 0, temp = 0;
	ifstream inFile(preampFile);

	ipCoefficients.clear();
	opCoefficients.clear();

	(*dispFuncPtr)(ctsiObjectPtr, indent);
	cout << "Importing preamplifier specifications:" << endl;
	(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
	cout << "Opening " << preampFile << "..." << endl;

	// Check file status
	if (!inFile.good())
	{
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cerr << "Unable to open preamplifier coefficient file " << preampFile << endl;
		return 1;
	}
	else
	{
		inFile >> caption >> samplingInterval;
		inFile >> caption >> numNumCoeffs;
		inFile >> caption;
		for (int i = 0; i < numNumCoeffs; i++)
		{
			inFile >> temp;
			ipCoefficients.push_back(temp);
		}

		inFile >> caption >> numDenCoeffs;
		inFile >> caption;
		for (int i = 0; i < numDenCoeffs; i++)
		{
			inFile >> temp;
			opCoefficients.push_back(temp);
		}

		inFile.close();

		if ((numNumCoeffs == 0) || (numDenCoeffs == 0))
		{
			(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
			cerr << "Error: discrete filter model for preamplifier cannot have zero coefficients for input or output samples." << endl;
			return 1;
		}
		else
		{
			// Normalize all coefficients by the coefficient for y[n]
			for (int i = 0; i < (int) ipCoefficients.size(); i++)
				ipCoefficients[i] = ipCoefficients[i]/opCoefficients[0];

			for (int i = 1; i < (int) opCoefficients.size(); i++)
				opCoefficients[i] = opCoefficients[i]/opCoefficients[0];
			opCoefficients[0] = 1;
		}

		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "Done." << endl;
	}
	return 0;
}

//===================================================================
// This function preamplifies the input signal by applying digital
// filtering with the preamplifier difference equation coefficients.
vector<double> Preamplification::preamplify(vector<double> &inSig)
{
	vector<double> outSig;
	double currSample = 0;

	// Loop over all input samples
	for (int i = 0; i < (int) inSig.size(); i++)
	{
		currSample = 0;

		// Loop over x[n - j] coefficients
		for (int j = 0; j < (int) ipCoefficients.size(); j++)
		{
			if (j <= i)
                currSample += ipCoefficients[j]*inSig[i - j];
			else
				break;
		}

		// Loop over y[n - j] coefficients
		for (int j = 1; j < (int) opCoefficients.size(); j++)
		{
			if (j <= i)
                currSample -= opCoefficients[j]*outSig[i - j];
			else
				break;
		}

		outSig.push_back(currSample);
	}

	return outSig;
}
