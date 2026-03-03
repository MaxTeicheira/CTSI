#ifndef DETECTORSTRUCTS_H
#define DETECTORSTRUCTS_H

#include <vector>

using std::vector;

//===================================================================
// Structure storing detector specifications.
struct detectorSpecStruct{
	double L;								// Thickness of detector in cm
	double W;								// Width of detector in cm
	int NUM_ANODES;							// Number of anodes on detector
	int NUM_CATHODES;						// Number of cathodes on detector
	double ANODE_PITCH;						// Anode pitch in mm
	double ANODE_WIDTH;						// Anode width in um
	double CATHODE_PITCH;					// Cathode pitch in mm
	double CATHODE_WIDTH;					// Cathode width in um
	double BIAS;							// Anode potential relative to the cathode in V/cm
	double CZT_W_FACTOR;					// W factor (ionization energy) in eV
	double DIFFUSION_E;						// Electron diffusion in cm^2/s
	double DIFFUSION_H;						// Hole diffusion in cm^2/s
	double MU_E;							// Electron mobility in cm^2/Vs
	double MU_H;							// Hole mobility in cm^2/Vs
	double TAU_E;							// Electron lifetime in s
	double TAU_H;							// Hole lifetime in s
	double FANO;							// Fano factor used for blurring number of electron-hole pairs generated
	double EPSILON_R;						// Relative permitivity
	double TEMP;							// Absolute temperature in K
	double LAMBDA_E;						// Electron mean drift length in cm
	double LAMBDA_H;						// Hole mean drift length in cm
};

// Structure storing the 2D electric field within the detector.
struct vecField2D{
	vector<double> xPos;					// X grid points
	vector<double> zPos;					// Z grid points
	vector<vector<double> > xComp;			// X component value
	vector<vector<double> > zComp;			// Z component value
};

// Structure storing #D scalar weighting potential field within the detector.
struct scalarField3D{
	vector<double> xPos;					// X grid points
	vector<double> yPos;					// Y grid points
	vector<double> zPos;					// Z grid points
	vector<vector<vector<double> > > scalar;// Value at grid point
};

# endif //DETECTORSTRUCTS_H
