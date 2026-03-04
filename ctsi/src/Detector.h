#ifndef DETECTOR_H
#define DETECTOR_H

#include "configParamsStruct.h"
#include "detectorStructs.h"
#include "phiLogStruct.h"
#include "trailLogStruct.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using std::vector;

#define PI 3.14159265358979323          // Pi

class Detector
{
	public:
		Detector();
		~Detector();
		void			anodeCollection(struct trailLog &chargeTrailsIn, int &currSimStep, vector<int> &anCollectorOneHot, bool caution);
		void			cathodeCollection(struct trailLog &chargeTrailsIn, int &currSimStep, vector<int> &caCollectorOneHot, bool caution);
		struct			detectorSpecStruct&	getDetectorSpecs();
		void			getEField(struct configParamStruct &configParams, double const &xVal, const double &zVal, vector<double> &eFieldIn);
		int				getEventTrailsSize(vector<vector<struct trailLog> > &eventTrailsIn, int &numChargeElems, vector<vector<int> > &trailSizeE, vector<vector<int> > &trailSizeH, bool caution);
		double			getWeightingFn2D(double W, double xValue, double zValue);
		double			getWeightingFn3D(struct configParamStruct &configParams, double xValue, double yValue, double zValue, int electrode);
		vector<double>	getXOffset();
		vector<double>	getYOffset();
		int				importDetectorSpecs(void* ctsiObjectPtr, int indent, void (*dispFuncPtr)(void* , int), char *detectorSpecFile);
		int				importEField(void* ctsiObjectPtr, int indent, void (*dispFuncPtr)(void*, int), char *eFieldFile, struct configParamStruct &configParams);
		int				importPhiAnode(void* ctsiObjectPtr, int indent, void (*dispFuncPtr)(void*, int), char *weightPotentialFile);
		int				importPhiCathode(void* ctsiObjectPtr, int indent, void (*dispFuncPtr)(void*, int), char *weightPotentialFile);
		void			printEFieldToConsole();
		void			querryPhiMat();

	protected:
		// Methods
		void	calcOffset();
		double	coth(double x);
		int		importPhi(void* ctsiObjectPtr, int indent, void (*dispFuncPtr)(void*, int), char *weightPotentialFile, struct scalarField3D &phiW3,
						  double &phiMinX, double &phiMaxX, double &phiMinY, double &phiMaxY, double &phiMinZ, double &phiMaxZ);

		// Attributes
		struct	detectorSpecStruct	detectorSpecs;	// Structure storing detector and material specifications
		struct	vecField2D			eField;			// Structure storing 2D electric field
		struct	scalarField3D		phiW3Anode;		// Structure storing the 3D anode weighting potential
		struct	scalarField3D		phiW3Cathode;	// Structure storing the 3D cathode weighting potential
		vector<double>				xOffset;		// Vector storing the x offset of each anode strip's center line from the detector center
		vector<double>				yOffset;		// Vector storing the y offset of each cathode strip's center line from the detector center
		double						anPhiMinX;		// Minimum x position in phiW3Anode
		double						anPhiMaxX;		// Maximum x position in phiW3Anode
		double						anPhiMinY;		// Minimum y position in phiW3Anode
		double						anPhiMaxY;		// Maximum y position in phiW3Anode
		double						anPhiMinZ;		// Minimum z position in phiW3Anode
		double						anPhiMaxZ;		// Maximum z position in phiW3Anode
		double						caPhiMinX;		// Minimum x position in phiW3Cathode
		double						caPhiMaxX;		// Maximum x position in phiW3Cathode
		double						caPhiMinY;		// Minimum y position in phiW3Cathode
		double						caPhiMaxY;		// Maximum y position in phiW3Cathode
		double						caPhiMinZ;		// Minimum z position in phiW3Cathode
		double						caPhiMaxZ;		// Maximum z position in phiW3Cathode
};

#endif // DETECTOR_H
