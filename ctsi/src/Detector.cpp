#include "Detector.h"

#include <limits>

using namespace std;

#define NUM_DET_SPEC_PARAMS_TO_IMPORT 19

Detector::Detector()
{
	phiW3Anode.scalar.clear();
	phiW3Anode.xPos.clear();
	phiW3Anode.yPos.clear();
	phiW3Anode.zPos.clear();

	phiW3Cathode.scalar.clear();
	phiW3Cathode.xPos.clear();
	phiW3Cathode.yPos.clear();
	phiW3Cathode.zPos.clear();

	eField.xPos.clear();
	eField.zPos.clear();
	eField.xComp.clear();
	eField.zComp.clear();

	xOffset.clear();
	yOffset.clear();

	detectorSpecs.L = 0;
	detectorSpecs.W = 0;
	detectorSpecs.NUM_ANODES = 0;
	detectorSpecs.NUM_CATHODES = 0;
	detectorSpecs.ANODE_PITCH = 0;
	detectorSpecs.ANODE_WIDTH = 0;
	detectorSpecs.CATHODE_PITCH = 0;
	detectorSpecs.CATHODE_WIDTH = 0;
	detectorSpecs.BIAS = 0;
	detectorSpecs.CZT_W_FACTOR = 0;
	detectorSpecs.DIFFUSION_E = 0;
	detectorSpecs.DIFFUSION_H = 0;
	detectorSpecs.MU_E = 0;
	detectorSpecs.MU_H = 0;
	detectorSpecs.TAU_E = 0;
	detectorSpecs.TAU_H = 0;
	detectorSpecs.FANO = 0;
	detectorSpecs.EPSILON_R = 0;
	detectorSpecs.TEMP = 0;
	detectorSpecs.LAMBDA_E = 0;
	detectorSpecs.LAMBDA_H = 0;

	anPhiMinX = 0;
	anPhiMaxX = 0;
	anPhiMinY = 0;
	anPhiMaxY = 0;
	anPhiMinZ = 0;
	anPhiMaxZ = 0;
	caPhiMinX = 0;
	caPhiMaxX = 0;
	caPhiMinY = 0;
	caPhiMaxY = 0;
	caPhiMinZ = 0;
	caPhiMaxZ = 0;
}

Detector::~Detector()
{
}

//===================================================================
// Given the final location of a charge element, this function sets
// the anCollectorOneHot flag vector to indicate which anode actually
// collected the charge. It is possible for none of the anodes to
// collect the charge.
void Detector::anodeCollection(struct trailLog &chargeTrailsIn, int &currSimStep, vector<int> &anCollectorOneHot, bool caution)
{
	double xValue = chargeTrailsIn.xPosE[currSimStep];

	// Logic:
	// The xValue*10 is to convert from cm to mm
	//
	// xValue*10 + detectorSpecs.ANODE_PITCH*0.5 shifts all interactions between n-pitch/2 and n+pitch/2
	// to between n + pitch. n is the center of the n th anode, and pitch/2 = 0.5 in this case. 
	// 
	// floor() sweeps everything in the range n + pitch to n. This essentially identifies which anode an
	// interaction would trigger.
	//
	// abs() takes the absolute value of the difference between an interaction and the position of the
	// center of the anode whose pitch the interaction fell in.
	//
	// detectorSpecs.ANODE_WIDTH*0.5 makes sure the charge landed on the electrode.
	if (abs(xValue*10 - floor(xValue*10 + detectorSpecs.ANODE_PITCH*0.5)) <= detectorSpecs.ANODE_WIDTH*0.5)
	{
		// Offsets the actual position to produce a 0-based vector index.
		xValue = xValue*10 + detectorSpecs.ANODE_PITCH*(((double) detectorSpecs.NUM_ANODES)/2);
		if (xValue < 0)
			xValue = 0;
		int anodeInd = (int) floor(xValue/detectorSpecs.ANODE_PITCH);

		// Account for discretization inaccuracies
		if ((anodeInd > detectorSpecs.NUM_ANODES - 1) && (anodeInd < detectorSpecs.NUM_ANODES))
			anodeInd = detectorSpecs.NUM_ANODES - 1;

		// Check index range
		if (caution)
		{
			if ((anodeInd < 0) || (anodeInd >= (int) anCollectorOneHot.size()))
				throw 1;

			if (abs(chargeTrailsIn.xPosE[currSimStep]-xOffset[anodeInd]) > detectorSpecs.ANODE_PITCH/20)
				throw 2;
		}

		anCollectorOneHot[anodeInd] = 1;
	}
}

//===================================================================
// This function determines the x offset of each anode strip and y
// offset of each cathode strip from the center of the detector, and
// stores the pffset in xOffset and yOffset (in cm).
void Detector::calcOffset()
{
	// The anchor point is essentially the x coordinate of the anode
	// strip with the most negative x coordinate, which, if there are
	// n anodes, is n/2-0.5 pitch distances negative of 0, where 0 is
	// the center.

	// Anodes
	double anchor = -(((double) detectorSpecs.NUM_ANODES)/2 - 0.5)*detectorSpecs.ANODE_PITCH;
	for (int i = 0; i < detectorSpecs.NUM_ANODES; i++)
		xOffset.push_back((anchor + i*detectorSpecs.ANODE_PITCH)*0.1);

	// Cathodes
	anchor = -(((double) detectorSpecs.NUM_CATHODES)/2 - 0.5)*detectorSpecs.CATHODE_PITCH;
	for (int i = 0; i < detectorSpecs.NUM_CATHODES; i++)
		yOffset.push_back((anchor + i*detectorSpecs.CATHODE_PITCH)*0.1);
}

//===================================================================
// Given the final location of a charge element, this function sets
// the caCollectorOneHot flag vector to indicate which cathode actually
// collected the charge.
void Detector::cathodeCollection(struct trailLog &chargeTrailsIn, int &currSimStep, vector<int> &caCollectorOneHot, bool caution)
{
	double yValue = chargeTrailsIn.yPosH[currSimStep];

	// Logic:
	// The yValue*10 is to convert from cm to mm
	//
	// yValue*10 + detectorSpecs.CATHODE_PITCH*0.5 shifts all interactions between n-pitch/2 and n+pitch/2
	// to between n + pitch. n is the center of the n th cathode, and pitch/2 = 0.5 in this case. 
	// 
	// floor() sweeps everything in the range n + pitch to n. This essentially identifies which cathode an
	// interaction would trigger.
	//
	// abs() takes the absolute value of the difference between an interaction and the position of the
	// center of the cathode whose pitch the interaction fell in.
	//
	// detectorSpecs.CATHODE_WIDTH*0.5 makes sure the charge landed on the electrode.
	if (abs(yValue*10 - floor(yValue*10 + detectorSpecs.CATHODE_PITCH*0.5)) <= detectorSpecs.CATHODE_WIDTH*0.5)
	{
		// Offsets the actual position to produce a 0-based vector index.
		yValue = yValue*10 + detectorSpecs.CATHODE_PITCH*(((double) detectorSpecs.NUM_CATHODES)/2);
		if (yValue < 0)
			yValue = 0;
		int cathodeInd = (int) floor(yValue/detectorSpecs.CATHODE_PITCH);

		// Account for discretization inaccuracies
		if ((cathodeInd > detectorSpecs.NUM_CATHODES - 1) && (cathodeInd < detectorSpecs.NUM_CATHODES))
			cathodeInd = detectorSpecs.NUM_CATHODES - 1;

		// Check index range
		if (caution)
		{
			if ((cathodeInd < 0) || (cathodeInd >= (int) caCollectorOneHot.size()))
				throw 3;

			if (abs(chargeTrailsIn.yPosH[currSimStep]-yOffset[cathodeInd]) > detectorSpecs.CATHODE_PITCH/20)
				throw 4;
		}
		caCollectorOneHot[cathodeInd] = 1;
	}
}

//===================================================================
// This function is the hyperbolic cotangent function
double Detector::coth(double x)
{
	return (exp(2*x)+1)/(exp(2*x)-1);
}

//===================================================================
// This function returns the detector specification structure
struct detectorSpecStruct& Detector::getDetectorSpecs()
{
	return detectorSpecs;
}

//===================================================================
// This function returns the interpolated electric field components
// in the x and z direction.
void Detector::getEField(struct configParamStruct &configParams, double const &xValIn, double const &zValIn, vector<double> &eFieldIn)
{
	// Convert from cm into mm
	double xVal = xValIn*10;
	double zVal = zValIn*10;

	// Correct for boundary cases
	if(xVal < -detectorSpecs.W*10*0.5)
		xVal = -detectorSpecs.W*10*0.5;
	else if(xVal > (detectorSpecs.W*10*0.5 - configParams.E_FIELD_GRID_SPACE_X))
		xVal = detectorSpecs.W*10*0.5 - configParams.E_FIELD_GRID_SPACE_X;
	if(zVal < 0)
		zVal = 0;
	else if(zVal > detectorSpecs.L*10)
		zVal = detectorSpecs.L*10;

	// Work out the indices of the E matrix that correspond to
	// the x and z locations
	int x_lower_ind = (int) floor((xVal + detectorSpecs.W*10*0.5)/configParams.E_FIELD_GRID_SPACE_X);
	int x_upper_ind = (int) ceil((xVal + detectorSpecs.W*10*0.5)/configParams.E_FIELD_GRID_SPACE_X);
	int z_lower_ind = (int) floor(zVal/configParams.E_FIELD_GRID_SPACE_Z);
	int z_upper_ind = (int) ceil(zVal/configParams.E_FIELD_GRID_SPACE_Z);

	double x_lower = eField.xPos[x_lower_ind];
	double x_upper = eField.xPos[x_upper_ind];
	double z_lower = eField.zPos[z_lower_ind];
	double z_upper = eField.zPos[z_upper_ind];

	// Find E field values at four corners
	// The electric field values are converted from V/m to V/cm
	double E11x = eField.xComp[x_lower_ind][z_lower_ind]*0.01;
	double E21x = eField.xComp[x_upper_ind][z_lower_ind]*0.01;
	double E22x = eField.xComp[x_upper_ind][z_upper_ind]*0.01;
	double E12x = eField.xComp[x_lower_ind][z_upper_ind]*0.01;
	double E11z = eField.zComp[x_lower_ind][z_lower_ind]*0.01;
	double E21z = eField.zComp[x_upper_ind][z_lower_ind]*0.01;
	double E22z = eField.zComp[x_upper_ind][z_upper_ind]*0.01;
	double E12z = eField.zComp[x_lower_ind][z_upper_ind]*0.01;

	// Bilinear interpolation
	double E1x, E2x, E1z, E2z, Ex_interp, Ez_interp;
	if(x_lower_ind == x_upper_ind) {
		E1x = E11x;
		E2x = E12x;
		E1z = E11z;
		E2z = E12z;
	}
	else {
		E1x = E11x + (E21x - E11x)*(xVal - x_lower)/(x_upper - x_lower);
		E2x = E12x + (E22x - E12x)*(xVal - x_lower)/(x_upper - x_lower);
		E1z = E11z + (E21z - E11z)*(xVal - x_lower)/(x_upper - x_lower);
		E2z = E12z + (E22z - E12z)*(xVal - x_lower)/(x_upper - x_lower);
	}

	if(z_lower_ind == z_upper_ind) {
		Ex_interp = E1x;
		Ez_interp = E1z;
	}
	else {
		Ex_interp = E1x + (E2x - E1x)*(zVal - z_lower)/(z_upper - z_lower);
		Ez_interp = E1z + (E2z - E1z)*(zVal - z_lower)/(z_upper - z_lower);
	}

	eFieldIn.push_back(Ex_interp);
	eFieldIn.push_back(Ez_interp);
}

//===================================================================
// This function finds the size of all trail logs in an event, stores
// them in trailSizeE and trailSizeH and returns the length of the
// longest trail (in terms of number of simulation steps).
int Detector::getEventTrailsSize(vector<vector<struct trailLog> > &eventTrailsIn, int &numChargeElems, vector<vector<int> > &trailSizeE, vector<vector<int> > &trailSizeH, bool caution)
{
	// Declare and initialize maximum trail length
	int numIntrxn = (int) eventTrailsIn.size();
	size_t maxSize = eventTrailsIn[0][0].xPosE.size();
	vector<int> tempE, tempH;
	size_t sizeE = 0, sizeH = 0;

	// Iterate over all interactions
	for(int intrxnIndex = 0; intrxnIndex < numIntrxn; intrxnIndex++)
	{
		tempE.clear();
		tempH.clear();

		// Iterate over all charge elements
		for(int chargeElem = 0; chargeElem < numChargeElems; chargeElem++)
		{
			sizeE = eventTrailsIn[intrxnIndex][chargeElem].xPosE.size();
			sizeH = eventTrailsIn[intrxnIndex][chargeElem].xPosH.size();

			tempE.push_back((int) sizeE);
			tempH.push_back((int) sizeH);

			// Run in cautious mode
			if (caution)
			{
				// While we're at it, check that all vectors belonging
				// to the same charge element has the same length.
				if ((sizeE != eventTrailsIn[intrxnIndex][chargeElem].yPosE.size()) ||
					(sizeE != eventTrailsIn[intrxnIndex][chargeElem].zPosE.size()) ||
					(sizeE != eventTrailsIn[intrxnIndex][chargeElem].qETrapped.size()) ||
					(sizeE != eventTrailsIn[intrxnIndex][chargeElem].qEMobile.size()))
					throw 0;

				if ((sizeH != eventTrailsIn[intrxnIndex][chargeElem].yPosH.size()) ||
					(sizeH != eventTrailsIn[intrxnIndex][chargeElem].zPosH.size()) ||
					(sizeH != eventTrailsIn[intrxnIndex][chargeElem].qHTrapped.size()) ||
					(sizeH != eventTrailsIn[intrxnIndex][chargeElem].qHMobile.size()))
					throw 0;
			}

			// Check for max
			if(sizeE > maxSize)
				maxSize = sizeE;
			if(sizeH > maxSize)
				maxSize = sizeH;
		}

		// Log
		trailSizeE.push_back(tempE);
		trailSizeH.push_back(tempH);
	}

	return (int) maxSize;
}

//===================================================================
// This function calculates the weighting function value
// for a specific x- and z-value (W is in cm). x is transverse to the
// electrode pitch, z is measured from the electrode plane.
double Detector::getWeightingFn2D(double W, double xValue, double zValue)
{
	double temp = 0;

	// Correcting for boundary case in weighting function
	// Boundary check on xValue is tricky because xValue is shifted,
	// so the bounds depend on which anode we are on.
	if(zValue == 0)
		zValue = 1e-6;
	else if (zValue > detectorSpecs.L)
		zValue = detectorSpecs.L;
    
	// Matlab equivalent code
	// myArg1Y = coth((pi/2/L)*(x - W/2))*tan(pi*zArray/2/L);
	double myArg1Y = coth((PI/2/detectorSpecs.L)*(xValue - W/2))*tan(PI*zValue/2/detectorSpecs.L);

	// myArg1X = sign(myArg1Y);
	double myArg1X = myArg1Y > 0 ? 1 : (myArg1Y < 0 ? -1 : 0);
	
	// myArg1Y = abs(myArg1Y);
	myArg1Y = abs(myArg1Y);
	
	// myArg2Y = coth((pi/2/L)*(x + W/2))*tan(pi*zArray/2/L);
	double myArg2Y = coth((PI/2/detectorSpecs.L)*(xValue + W/2))*tan(PI*zValue/2/detectorSpecs.L);
	
	// myArg2X = sign(myArg2Y);
	double myArg2X = myArg2Y > 0 ? 1 : (myArg2Y < 0 ? -1 : 0);
	
	// myArg2Y = abs(myArg2Y);
	myArg2Y = abs(myArg2Y);
	
	// phiW(i, :) = max(0, (atan2(myArg1Y, myArg1X) - atan2(myArg2Y, myArg2X))/pi);
	temp = (atan2(myArg1Y, myArg1X) - atan2(myArg2Y, myArg2X))/PI;
	double retVal = temp > 0 ? temp : 0;

	return retVal;
}

//===================================================================
// This function retrieves the weighting potential value for a given
// (x, y, z) position as stored in phiW3Anode and phiW3Cathode. It
// performs trilinear interpolation. The phi data structures exploit
// the symmetry of the electrodes, so they store just one quadrant of
// the full scalar function. electrode = 0 is anode, electrode = 1 is
// cathode.
double Detector::getWeightingFn3D(struct configParamStruct &configParams, double xValue, double yValue, double zValue, int electrode)
{
	double absXValue = 0, absYValue = 0, w1 = 0, w2 = 0;
	double xLo = 0, xHi = 0, yLo = 0, yHi = 0, zLo = 0, zHi = 0;
	double xLyLzL = 0, xLyLzH = 0, xLyHzL = 0, xLyHzH = 0, xHyLzL = 0, xHyLzH = 0, xHyHzL = 0, xHyHzH = 0;
	double loFaceupSide = 0, loFaceloSide = 0, hiFaceupSide = 0, hiFaceloSide = 0, loFaceInterp = 0, hiFaceInterp = 0;
	int xLoInd = 0, xHiInd = 0, yLoInd = 0, yHiInd = 0, zLoInd = 0, zHiInd = 0;

	// Note that the first points to the left and right of 0 in the
	// x and y directions are not separated by the normal step size.
	
	// zValue between 0 and anPhiMinZ is treated as being at anPhiMinZ

	// Make sure the input is within bound
	// 0 = anode, 1 = cathode
	if (electrode == 0)
	{
		xValue = (xValue < -anPhiMaxX) ? -anPhiMaxX : ((xValue > anPhiMaxX) ? anPhiMaxX : xValue);
		yValue = (yValue < -anPhiMaxY) ? -anPhiMaxY : ((yValue > anPhiMaxY) ? anPhiMaxY : yValue);
		zValue = (zValue < 0) ? 0 : ((zValue > anPhiMaxZ) ? anPhiMaxZ : zValue);

		absXValue = abs(xValue);
		absYValue = abs(yValue);

		//===================================================================
		// Determine x interval
		// Check if the x position is in [-anPhiMinX anPhiMinX]
		if (absXValue <= anPhiMinX)
		{
			xLoInd = 0;
			xHiInd = 0;
		}
		else
		{
			xLoInd = (int) floor((absXValue - anPhiMinX)/configParams.AN_PHI_GRID_SPACE_X);
			xHiInd = (int) ceil((absXValue - anPhiMinX)/configParams.AN_PHI_GRID_SPACE_X);
		}

		//===================================================================
		// Determine y interval
		// Check if the y position is in [-anPhiMinY anPhiMinY]
		if (absYValue <= anPhiMinY)
		{
			yLoInd = 0;
			yHiInd = 0;
		}
		else
		{
			yLoInd = (int) floor((absYValue - anPhiMinY)/configParams.AN_PHI_GRID_SPACE_Y);
			yHiInd = (int) ceil((absYValue - anPhiMinY)/configParams.AN_PHI_GRID_SPACE_Y);
		}

		//===================================================================
		// Determine z interval
		// Check if the z position is in [0 anPhiMinZ]
		if (abs(zValue) <= anPhiMinZ)
		{
			zLoInd = 0;
			zHiInd = 0;
		}
		else
		{
			zLoInd = (int) floor((abs(zValue) - anPhiMinZ)/configParams.AN_PHI_GRID_SPACE_Z);
			zHiInd = (int) ceil((abs(zValue) - anPhiMinZ)/configParams.AN_PHI_GRID_SPACE_Z);
		}

		// Cube corners
		xLo = phiW3Anode.xPos[xLoInd];
		xHi = phiW3Anode.xPos[xHiInd];
		yLo = phiW3Anode.yPos[yLoInd];
		yHi = phiW3Anode.yPos[yHiInd];
		zLo = phiW3Anode.zPos[zLoInd];
		zHi = phiW3Anode.zPos[zHiInd];

		// Weighting potential at cube corners
		xLyLzL = phiW3Anode.scalar[xLoInd][yLoInd][zLoInd];
		xLyLzH = phiW3Anode.scalar[xLoInd][yLoInd][zHiInd];
		xLyHzL = phiW3Anode.scalar[xLoInd][yHiInd][zLoInd];
		xLyHzH = phiW3Anode.scalar[xLoInd][yHiInd][zHiInd];
		xHyLzL = phiW3Anode.scalar[xHiInd][yLoInd][zLoInd];
		xHyLzH = phiW3Anode.scalar[xHiInd][yLoInd][zHiInd];
		xHyHzL = phiW3Anode.scalar[xHiInd][yHiInd][zLoInd];
		xHyHzH = phiW3Anode.scalar[xHiInd][yHiInd][zHiInd];

		//===================================================================
		// Interpolate in the x direction
		if ((xLoInd == 0) && (xHiInd == 0))
		{
			// In this case, as long as w1 + w2 = 1, it doesn't matter what
			// their values are otherwise because
			// xLyHzL = xHyHzL
			// xLyLzL = xHyLzL
			// xLyHzH = xHyHzH
			// xLyLzH = xHyLzH
			w1 = 0.5;
			w2 = 0.5;
		}
		else
		{
			w1 = (absXValue - xLo)/configParams.AN_PHI_GRID_SPACE_X;
			w2 = 1 - w1;
		}
		// Lower face upper side
		loFaceupSide = w2*xLyHzL + w1*xHyHzL;
		// Lower face lower side
		loFaceloSide = w2*xLyLzL + w1*xHyLzL;
		// Upper face upper side
		hiFaceupSide = w2*xLyHzH + w1*xHyHzH;
		// Upper face lower side
		hiFaceloSide = w2*xLyLzH + w1*xHyLzH;

		//===================================================================
		// Interpolate in the y direction
		if ((yLoInd == 0) && (yHiInd == 0))
		{
			// Again, it doesn't matter as long as w1 + w2 = 1.
			w1 = 0.5;
			w2 = 0.5;
		}
		else
		{
			w1 = (absYValue - yLo)/configParams.AN_PHI_GRID_SPACE_Y;
			w2 = 1 - w1;
		}
		loFaceInterp = w2*loFaceloSide + w1*loFaceupSide;
		hiFaceInterp = w2*hiFaceloSide + w1*hiFaceupSide;

		//===================================================================
		// Interpolate in the z direction
		if ((zLoInd == 0) && (zHiInd == 0))
		{
			// Again, it doesn't matter as long as w1 + w2 = 1.
			w1 = 0.5;
			w2 = 0.5;
		}
		else
		{
			w1 = (zValue - zLo)/configParams.AN_PHI_GRID_SPACE_Z;
			w2 = 1 - w1;
		}

		return w2*loFaceInterp + w1*hiFaceInterp;

	}
	// Cathode
	else
	{
		xValue = (xValue < -caPhiMaxX) ? -caPhiMaxX : ((xValue > caPhiMaxX) ? caPhiMaxX : xValue);
		yValue = (yValue < -caPhiMaxY) ? -caPhiMaxY : ((yValue > caPhiMaxY) ? caPhiMaxY : yValue);
		zValue = (zValue < -caPhiMaxZ) ? -caPhiMaxZ : ((zValue > caPhiMaxZ) ? caPhiMaxZ : zValue);

		absXValue = abs(xValue);
		absYValue = abs(yValue);

		//===================================================================
		// Determine x interval
		// Check if the x position is in [-caPhiMinX caPhiMinX]
		if (absXValue <= caPhiMinX)
		{
			xLoInd = 0;
			xHiInd = 0;
		}
		else
		{
			xLoInd = (int) floor((absXValue - caPhiMinX)/configParams.CA_PHI_GRID_SPACE_X);
			xHiInd = (int) ceil((absXValue - caPhiMinX)/configParams.CA_PHI_GRID_SPACE_X);
		}

		//===================================================================
		// Determine y interval
		// Check if the y position is in [-caPhiMinY caPhiMinY]
		if (absYValue <= caPhiMinY)
		{
			yLoInd = 0;
			yHiInd = 0;
		}
		else
		{
			yLoInd = (int) floor((absYValue - caPhiMinY)/configParams.CA_PHI_GRID_SPACE_Y);
			yHiInd = (int) ceil((absYValue - caPhiMinY)/configParams.CA_PHI_GRID_SPACE_Y);
		}

		//===================================================================
		// Determine z interval
		// Check if the z position is in [0 caPhiMinZ]
		if (abs(zValue) <= caPhiMinZ)
		{
			zLoInd = 0;
			zHiInd = 0;
		}
		else
		{
			zLoInd = (int) floor((abs(zValue) - caPhiMinZ)/configParams.CA_PHI_GRID_SPACE_Z);
			zHiInd = (int) ceil((abs(zValue) - caPhiMinZ)/configParams.CA_PHI_GRID_SPACE_Z);
		}

		// Cube corners
		xLo = phiW3Cathode.xPos[xLoInd];
		xHi = phiW3Cathode.xPos[xHiInd];
		yLo = phiW3Cathode.yPos[yLoInd];
		yHi = phiW3Cathode.yPos[yHiInd];
		zLo = phiW3Cathode.zPos[zLoInd];
		zHi = phiW3Cathode.zPos[zHiInd];

		// Weighting potential at cube corners
		xLyLzL = phiW3Cathode.scalar[xLoInd][yLoInd][zLoInd];
		xLyLzH = phiW3Cathode.scalar[xLoInd][yLoInd][zHiInd];
		xLyHzL = phiW3Cathode.scalar[xLoInd][yHiInd][zLoInd];
		xLyHzH = phiW3Cathode.scalar[xLoInd][yHiInd][zHiInd];
		xHyLzL = phiW3Cathode.scalar[xHiInd][yLoInd][zLoInd];
		xHyLzH = phiW3Cathode.scalar[xHiInd][yLoInd][zHiInd];
		xHyHzL = phiW3Cathode.scalar[xHiInd][yHiInd][zLoInd];
		xHyHzH = phiW3Cathode.scalar[xHiInd][yHiInd][zHiInd];

		//===================================================================
		// Interpolate in the x direction
		if ((xLoInd == 0) && (xHiInd == 0))
		{
			// In this case, as long as w1 + w2 = 1, it doesn't matter what
			// their values are otherwise because
			// xLyHzL = xHyHzL
			// xLyLzL = xHyLzL
			// xLyHzH = xHyHzH
			// xLyLzH = xHyLzH
			w1 = 0.5;
			w2 = 0.5;
		}
		else
		{
			w1 = (absXValue - xLo)/configParams.CA_PHI_GRID_SPACE_X;
			w2 = 1 - w1;
		}
		// Lower face upper side
		loFaceupSide = w2*xLyHzL + w1*xHyHzL;
		// Lower face lower side
		loFaceloSide = w2*xLyLzL + w1*xHyLzL;
		// Upper face upper side
		hiFaceupSide = w2*xLyHzH + w1*xHyHzH;
		// Upper face lower side
		hiFaceloSide = w2*xLyLzH + w1*xHyLzH;

		//===================================================================
		// Interpolate in the y direction
		if ((yLoInd == 0) && (yHiInd == 0))
		{
			// Again, it doesn't matter as long as w1 + w2 = 1.
			w1 = 0.5;
			w2 = 0.5;
		}
		else
		{
			w1 = (absYValue - yLo)/configParams.CA_PHI_GRID_SPACE_Y;
			w2 = 1 - w1;
		}
		loFaceInterp = w2*loFaceloSide + w1*loFaceupSide;
		hiFaceInterp = w2*hiFaceloSide + w1*hiFaceupSide;

		//===================================================================
		// Interpolate in the z direction
		if ((zLoInd == 0) && (zHiInd == 0))
		{
			// Again, it doesn't matter as long as w1 + w2 = 1.
			w1 = 0.5;
			w2 = 0.5;
		}
		else
		{
			w1 = (zValue - zLo)/configParams.CA_PHI_GRID_SPACE_Z;
			w2 = 1 - w1;
		}

		return w2*loFaceInterp + w1*hiFaceInterp;
	}
}

//===================================================================
// This function returns the xOffset member attribute
vector<double> Detector::getXOffset()
{
	return xOffset;
}

//===================================================================
// This function returns the yOffset member attribute
vector<double> Detector::getYOffset()
{
	return yOffset;
}

//===================================================================
// This function reads the detector specification file and
// initializes the detector specification data structure.
int Detector::importDetectorSpecs(void* ctsiObjectPtr, int indent, void (*dispFuncPtr)(void*, int), char *detectorSpecFile)
{
	char caption[256], comment[256];
	int lineNumber = 1;
	int numParamsImported = 0;
	
	(*dispFuncPtr)(ctsiObjectPtr, indent);
	cout << "Importing detector specifications:" << endl;
	(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
	cout << "Opening " << detectorSpecFile << "..." << endl;
	ifstream inFile(detectorSpecFile);

	try
	{
		// Read the configuration file
		if (!inFile.good())
		{
			(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
			cerr << "Unable to open detector specification file " << detectorSpecFile << endl;
			throw 0;
		}

		while(!inFile.eof())
		{
			comment[0] = '\0';

			try
			{
				switch (lineNumber)
				{
					// Caption is not used for display as a way to catch when parameter and
					// imported field don't match up.
					case 1:
						// Thickness of detector in cm
						inFile >> caption >> detectorSpecs.L;
						inFile.getline(comment, 256);
						if (!strcmp(comment, "\0"))
							throw 0;
						(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
						cout << "L: " << detectorSpecs.L << comment << endl;
						numParamsImported++;
						break;
					case 2:
						// Width of detector in cm
						inFile >> caption >> detectorSpecs.W;
						inFile.getline(comment, 256);
						if (!strcmp(comment, "\0"))
							throw 0;
						(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
						cout << "W: " << detectorSpecs.W << comment << endl;
						numParamsImported++;
						break;
					case 3:
						// Number of anodes on detector
						inFile >> caption >> detectorSpecs.NUM_ANODES;
						inFile.getline(comment, 256);
						if (!strcmp(comment, "\0"))
							throw 0;
						(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
						cout << "NUM_ANODES: " << detectorSpecs.NUM_ANODES << comment << endl;
						numParamsImported++;
						break;
					case 4:
						// Number of cathodes on detector
						inFile >> caption >> detectorSpecs.NUM_CATHODES;
						inFile.getline(comment, 256);
						if (!strcmp(comment, "\0"))
							throw 0;
						(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
						cout << "NUM_CATHODES: " << detectorSpecs.NUM_CATHODES << comment << endl;
						numParamsImported++;
						break;
					case 5:
						// Anode pitch in mm
						inFile >> caption >> detectorSpecs.ANODE_PITCH;
						inFile.getline(comment, 256);
						if (!strcmp(comment, "\0"))
							throw 0;
						(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
						cout << "ANODE_PITCH: " << detectorSpecs.ANODE_PITCH << comment << endl;
						numParamsImported++;
						break;
					case 6:
						// Anode width in um
						inFile >> caption >> detectorSpecs.ANODE_WIDTH;
						inFile.getline(comment, 256);
						if (!strcmp(comment, "\0"))
							throw 0;
						(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
						cout << "ANODE_WIDTH: " << detectorSpecs.ANODE_WIDTH << comment << endl;
						numParamsImported++;
						break;
					case 7:
						// Cathode pitch in mm
						inFile >> caption >> detectorSpecs.CATHODE_PITCH;
						inFile.getline(comment, 256);
						if (!strcmp(comment, "\0"))
							throw 0;
						(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
						cout << "CATHODE_PITCH: " << detectorSpecs.CATHODE_PITCH << comment << endl;
						numParamsImported++;
						break;
					case 8:
						// Cathode width in um
						inFile >> caption >> detectorSpecs.CATHODE_WIDTH;
						inFile.getline(comment, 256);
						if (!strcmp(comment, "\0"))
							throw 0;
						(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
						cout << "CATHODE_WIDTH: " << detectorSpecs.CATHODE_WIDTH << comment << endl;
						numParamsImported++;
						break;
					case 9:
						// Anode potential relative to the cathode in V/cm
						inFile >> caption >> detectorSpecs.BIAS;
						inFile.getline(comment, 256);
						if (!strcmp(comment, "\0"))
							throw 0;
						(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
						cout << "BIAS: " << detectorSpecs.BIAS << comment << endl;
						numParamsImported++;
						break;
					case 10:
						// W factor of CZT (ionization energy) in eV
						inFile >> caption >> detectorSpecs.CZT_W_FACTOR;
						inFile.getline(comment, 256);
						if (!strcmp(comment, "\0"))
							throw 0;
						(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
						cout << "CZT_W_FACTOR: " << detectorSpecs.CZT_W_FACTOR << comment << endl;
						numParamsImported++;
						break;
					case 11:
						// Electron diffusion in cm^2/s
						inFile >> caption >> detectorSpecs.DIFFUSION_E;
						inFile.getline(comment, 256);
						if (!strcmp(comment, "\0"))
							throw 0;
						(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
						cout << "DIFFUSION_E: " << detectorSpecs.DIFFUSION_E << comment << endl;
						numParamsImported++;
						break;
					case 12:
						// Hole diffusion in cm^2/s
						inFile >> caption >> detectorSpecs.DIFFUSION_H;
						inFile.getline(comment, 256);
						if (!strcmp(comment, "\0"))
							throw 0;
						(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
						cout << "DIFFUSION_H: " << detectorSpecs.DIFFUSION_H << comment << endl;
						numParamsImported++;
						break;
					case 13:
						// Electron mobility in CZT in cm^2/Vs
						inFile >> caption >> detectorSpecs.MU_E;
						inFile.getline(comment, 256);
						if (!strcmp(comment, "\0"))
							throw 0;
						(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
						cout << "MU_E: " << detectorSpecs.MU_E << comment << endl;
						numParamsImported++;
						break;
					case 14:
						// Hole mobility in CZT in cm^2/Vs
						inFile >> caption >> detectorSpecs.MU_H;
						inFile.getline(comment, 256);
						if (!strcmp(comment, "\0"))
							throw 0;
						(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
						cout << "MU_H: " << detectorSpecs.MU_H << comment << endl;
						numParamsImported++;
						break;
					case 15:
						// Electron lifetime in s
						inFile >> caption >> detectorSpecs.TAU_E;
						inFile.getline(comment, 256);
						if (!strcmp(comment, "\0"))
							throw 0;
						(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
						cout << "TAU_E: " << detectorSpecs.TAU_E << comment << endl;
						numParamsImported++;
						break;
					case 16:
						// Hole lifetime in s
						inFile >> caption >> detectorSpecs.TAU_H;
						inFile.getline(comment, 256);
						(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
						cout << "TAU_H: " << detectorSpecs.TAU_H << comment << endl;
						numParamsImported++;
						break;
					case 17:
						// Fano factor used for blurring number of electron-hole pairs generated
						inFile >> caption >> detectorSpecs.FANO;
						inFile.getline(comment, 256);
						if (!strcmp(comment, "\0"))
							throw 0;
						(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
						cout << "FANO: " << detectorSpecs.FANO << comment << endl;
						numParamsImported++;
						break;
					case 18:
						// Relative permitivity
						inFile >> caption >> detectorSpecs.EPSILON_R;
						inFile.getline(comment, 256);
						if (!strcmp(comment, "\0"))
							throw 0;
						(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
						cout << "EPSILON_R: " << detectorSpecs.EPSILON_R << comment << endl;
						numParamsImported++;
						break;
					case 19:
						// Absolute temperature in K
						inFile >> caption >> detectorSpecs.TEMP;
						inFile.getline(comment, 256);
						if (!strcmp(comment, "\0"))
							throw 0;
						(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
						cout << caption << ": " << detectorSpecs.TEMP << comment << endl;
						numParamsImported++;
						break;
					default:
						inFile >> caption;
				}
			}
			catch (int)
			{
				cout << endl;
				(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
				cerr << numParamsImported << " of " << NUM_DET_SPEC_PARAMS_TO_IMPORT << " required detector specification parameters imported, import aborted." << endl;
				throw 0;
			}
			lineNumber++;
		}
		inFile.close();

		if (numParamsImported != NUM_DET_SPEC_PARAMS_TO_IMPORT)
		{
			cout << endl;
			(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
			cerr << numParamsImported << " of " << NUM_DET_SPEC_PARAMS_TO_IMPORT << " required detector specification parameters imported, import aborted." << endl;
			throw 0;
		}

		detectorSpecs.LAMBDA_E = detectorSpecs.MU_E*detectorSpecs.BIAS*detectorSpecs.TAU_E; // Electron mean drift length in cm
		detectorSpecs.LAMBDA_H = detectorSpecs.MU_H*detectorSpecs.BIAS*detectorSpecs.TAU_H; // Hole mean drift length in cm
	}
	catch (int)
	{
		cout << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "Contents of the detector specification file must have the following format (comments are not optional):" << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "L VALUE                        // Thickness of detector in cm" << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "W VALUE                        // Width of detector in cm" << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "NUM_ANODES VALUE               // Number of anodes on detector" << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "NUM_CATHODES VALUE             // Number of cathodes on detector" << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "ANODE_PITCH VALUE              // Anode pitch in mm" << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "ANODE_WIDTH VALUE              // Anode width in um" << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "CATHODE_PITCH VALUE            // Cathode pitch in mm" << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "CATHODE_WIDTH VALUE            // Cathode width in um" << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "BIAS VALUE                     // Anode potential relative to the cathode in V/cm" << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "CZT_W_FACTOR VALUE             // W factor (ionization energy) in eV" << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "DIFFUSION_E VALUE              // Electron diffusion" << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "DIFFUSION_H VALUE              // Hole diffusion" << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "MU_E VALUE                     // Electron mobility in cm^2/Vs" << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "MU_H VALUE                     // Hole mobility in cm^2/Vs" << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "TAU_E VALUE                    // Electron lifetime in s" << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "TAU_H VALUE                    // Hole lifetime in s" << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "FANO VALUE                     // Fano factor used for blurring number of electron-hole pairs generated" << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "EPSILON_R VALUE                // Relative permitivity" << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "TEMP VALUE                     // Absolute temperature in K" << endl;
		cout << endl;
		return 1;
	}

	calcOffset();

	(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
	cout << "Done." << endl;

	return 0;
}

//===================================================================
// This function reads in the 2D electric field in ASCII format and
// stores the field data in a vecField2D structure.
int Detector::importEField(void* ctsiObjectPtr, int indent, void (*dispFuncPtr)(void*, int), char *eFieldFile, struct configParamStruct &configParams)
{
	vector<double> tempXVector, tempZVector;
	int numXPos = 0, numZPos = 0, temp = 0;
	double oldXPosVal = 0, xPosVal = 0, zPosVal = 0, xCompVal = 0, zCompVal = 0;
	ifstream inFile(eFieldFile);
	bool warned = false;

	// Flush structure
	eField.xPos.clear();
	eField.zPos.clear();
	eField.xComp.clear();
	eField.zComp.clear();
	tempXVector.clear();
	tempZVector.clear();

	(*dispFuncPtr)(ctsiObjectPtr, indent);
	cout << "Importing electric field matrix:" << endl;
	(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
	cout << "Opening " << eFieldFile << "..." << endl;

	// Check file status
	if (!inFile.good())
	{
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cerr << "Unable to open electric field matrix file " << eFieldFile << endl;
		return 1;
	}

	inFile >> numXPos >> numZPos >> temp >> temp;

	// Store first line of data
	inFile >> xPosVal >> zPosVal >> xCompVal >> zCompVal;
	tempXVector.push_back(xCompVal);
	tempZVector.push_back(zCompVal);
	eField.xPos.push_back(xPosVal);
	oldXPosVal = xPosVal;
	eField.zPos.push_back(zPosVal);

	// Fill the z position vector only once
	while (inFile >> xPosVal >> zPosVal >> xCompVal >> zCompVal)
	{
		if (xPosVal != oldXPosVal)
			break;

		tempXVector.push_back(xCompVal);
		tempZVector.push_back(zCompVal);
		eField.zPos.push_back(zPosVal);
	}

	// Catch inconsistent grid spacing specification
	int vecLength = (int) eField.zPos.size();
	if (abs((eField.zPos[vecLength - 1] - eField.zPos[vecLength - 2]) - configParams.E_FIELD_GRID_SPACE_Z) > 1e-8)
	{
		cout << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cout << "Warning: actual and specified z grid spacing may be inconsistent!" << endl;
		cout << endl;
	}

	do
	{
		// Add temporary vectors to data structure if the x position has changed
		if (xPosVal != oldXPosVal)
		{
			// Catch inconsistent grid spacing specification
			if (!warned)
			{
				if (abs((xPosVal - oldXPosVal) - configParams.E_FIELD_GRID_SPACE_X) > 1e-8)
				{
					warned = true;
					cout << endl;
					(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
					cout << "Warning: actual and specified x grid spacing may be inconsistent!" << endl;
					cout << endl;
				}
			}

			eField.xComp.push_back(tempXVector);
			eField.zComp.push_back(tempZVector);
			eField.xPos.push_back(xPosVal);
			oldXPosVal = xPosVal;
			tempXVector.clear();
			tempZVector.clear();
		}

		// If the x position is the same, then just add to temporary vector
		tempXVector.push_back(xCompVal);
		tempZVector.push_back(zCompVal);
	}
	while (inFile >> xPosVal >> zPosVal >> xCompVal >> zCompVal);

	// Push the last tempXVector and tempZVector's
	eField.xComp.push_back(tempXVector);
	eField.zComp.push_back(tempZVector);

	inFile.close();

	(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
	cout << "Done." << endl;

	return 0;
}

//===================================================================
// This function reads in the 3D weighting potential in ASCII format
// and stores the scalar field data in a scalarField3D structure.
int Detector::importPhi(void* ctsiObjectPtr, int indent, void (*dispFuncPtr)(void*, int), char *weightPotentialFile, struct scalarField3D &phiW3,
						double &phiMinX, double &phiMaxX, double &phiMinY, double &phiMaxY, double &phiMinZ, double &phiMaxZ)
{
	vector<double> tempZVector;
	vector<vector<double> > tempYVector;
	int numXPos = 0, numYPos = 0, numZPos = 0;
	double temp = 0, phiMax = 0;
	ifstream inFile(weightPotentialFile);

	// Flush structure
	phiW3.xPos.clear();
	phiW3.yPos.clear();
	phiW3.zPos.clear();
	phiW3.scalar.clear();
	tempYVector.clear();
	tempZVector.clear();

	phiMinX = numeric_limits<double>::max();
	phiMaxX = numeric_limits<double>::min();
	phiMinY = numeric_limits<double>::max();
	phiMaxY = numeric_limits<double>::min();
	phiMinZ = numeric_limits<double>::max();
	phiMaxZ = numeric_limits<double>::min();

	(*dispFuncPtr)(ctsiObjectPtr, indent);
	cout << "Importing weighting potential matrix:" << endl;
	(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
	cout << "Opening " << weightPotentialFile << "..." << endl;

	// Check file status
	if (!inFile.good())
	{
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cerr << "Unable to open weighting potential matrix file " << weightPotentialFile << endl;
		return 1;
	}

	try
	{
		// Read dimension of matrix
		if (inFile.eof())
			throw 0;
		inFile >> numXPos;
		
		if (inFile.eof())
			throw 0;
		inFile >> numYPos;
		
		if (inFile.eof())
			throw 0;
		inFile >> numZPos;

		// Read in all x grid points
		for (int i = 0; i < numXPos; i++)
		{
			if (inFile.eof())
				throw 0;
			inFile >> temp;
			phiW3.xPos.push_back(temp);

			if (temp < phiMinX)
				phiMinX = temp;
			if (temp > phiMaxX)
				phiMaxX = temp;
		}

		// Read in all y grid points
		for (int i = 0; i < numYPos; i++)
		{
			if (inFile.eof())
				throw 0;
			inFile >> temp;
			phiW3.yPos.push_back(temp);

			if (temp < phiMinY)
				phiMinY = temp;
			if (temp > phiMaxY)
				phiMaxY = temp;
		}

		// Read in all z grid points
		for (int i = 0; i < numZPos; i++)
		{
			if (inFile.eof())
				throw 0;
			inFile >> temp;
			phiW3.zPos.push_back(temp);

			if (temp < phiMinZ)
				phiMinZ = temp;
			if (temp > phiMaxZ)
				phiMaxZ = temp;
		}

		// Beware that the order of the raster read has to match the order in which the
		// matrix was written to file in Matlab.
		// Loop over all x positions (across electrode pitch)
		for (int xIndex = 0; xIndex < numXPos; xIndex++)
		{
			for (int k = 0; k < (int) tempYVector.size(); k++)
				tempYVector[k].clear();
			tempYVector.clear();

			// Loop over all y positions (along electrode length)
			for (int yIndex = 0; yIndex < numYPos; yIndex++)
			{
				tempZVector.clear();
				// Loop over all z positions along the detector thickness
				for (int zIndex = 0; zIndex < numZPos; zIndex++)
				{
					if (inFile.eof())
						throw 0;
					inFile >> temp;
					tempZVector.push_back(temp);
					if (abs(temp) > phiMax)
						phiMax = temp;
				}
				tempYVector.push_back(tempZVector);
			}
			phiW3.scalar.push_back(tempYVector);
		}

		// Normalize so the largets component has value 1
		for (int xIndex = 0; xIndex < numXPos; xIndex++)
		{
			// Loop over all y positions (along electrode length)
			for (int yIndex = 0; yIndex < numYPos; yIndex++)
			{
				// Loop over all z positions along the detector thickness
				for (int zIndex = 0; zIndex < numZPos; zIndex++)
					phiW3.scalar[xIndex][yIndex][zIndex] = phiW3.scalar[xIndex][yIndex][zIndex]/phiMax;
			}
		}

	}
	catch (int)
	{
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cerr << "An error occurred while importing from " << weightPotentialFile << ", import aborted." << endl;
		return 1;
	}

	inFile.close();

	(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
	cout << "Matrix size is " << numXPos << " x " << numYPos << " x " << numZPos << " (x, y, z)" << endl;
	(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
	cout << "Done." << endl;

	return 0;
}

//===================================================================
// This function imports the anode weighting potential.
int Detector::importPhiAnode(void* ctsiObjectPtr, int indent, void (*dispFuncPtr)(void*, int), char *weightPotentialFile)
{
	return importPhi(ctsiObjectPtr, indent, dispFuncPtr, weightPotentialFile, phiW3Anode, anPhiMinX, anPhiMaxX, anPhiMinY, anPhiMaxY, anPhiMinZ, anPhiMaxZ);
}

//===================================================================
// This function imports the cathode weighting potential.
int Detector::importPhiCathode(void* ctsiObjectPtr, int indent, void (*dispFuncPtr)(void*, int), char *weightPotentialFile)
{
	return importPhi(ctsiObjectPtr, indent, dispFuncPtr, weightPotentialFile, phiW3Cathode, caPhiMinX, caPhiMaxX, caPhiMinY, caPhiMaxY, caPhiMinZ, caPhiMaxZ);
}

//===================================================================
// This function prints the electric field matrix to console.
void Detector::printEFieldToConsole()
{
	unsigned int i = 0, j = 0;
	vector<double>::iterator it_vectorDouble1, it_vectorDouble2;

	cout << endl;
	cout << "Printing E field matrix to console:" << endl;
	
	for (it_vectorDouble1 = eField.xPos.begin(); it_vectorDouble1 != eField.xPos.end(); it_vectorDouble1++)
	{
		j = 0;
		for (it_vectorDouble2 = eField.zPos.begin(); it_vectorDouble2 != eField.zPos.end(); it_vectorDouble2++)
		{
			cout.precision(9);
			cout << fixed << "x pos: " << (*it_vectorDouble1) << " z pos: " << (*it_vectorDouble2);
			cout << scientific << " x comp: " << eField.xComp[i][j] << " z comp: " << eField.zComp[i][j] << endl;
			j++;
		}
		i++;
	}
}

//===================================================================
// This function prints the selected entry in the weighting potential
// matrix to console.
void Detector::querryPhiMat()
{
	int toExit = 0, choice = 0, xIndex = 0, yIndex = 0, zIndex = 0, maxXIndex = 0, maxYIndex = 0, maxZIndex = 0;

	while(!toExit)
	{
		// User prompt
		cout << "Enter 0 for anode weighting potential matrix or 1 for cathode weighting potential matrix: ";
		cin >> choice;

		if (!choice)
		{
			maxXIndex = int(phiW3Anode.xPos.size());
			maxYIndex = int(phiW3Anode.yPos.size());
			maxZIndex = int(phiW3Anode.zPos.size());
		}
		else
		{
			maxXIndex = int(phiW3Cathode.xPos.size());
			maxYIndex = int(phiW3Cathode.yPos.size());
			maxZIndex = int(phiW3Cathode.zPos.size());
		}

		cout << "Enter x index: ";
		cin >> xIndex;
		cout << "Enter y index: ";
		cin >> yIndex;
		cout << "Enter z index: ";
		cin >> zIndex;

		if ((0 <= xIndex) & (xIndex <= maxXIndex) & (0 <= yIndex) & (yIndex <= maxYIndex) & (0 <= zIndex) & (zIndex <= maxZIndex))
		{
			if (!choice)
				cout << "Anode phi value: " << phiW3Anode.scalar[xIndex][yIndex][zIndex] << endl;
			else
				cout << "Cathode phi value: " << phiW3Cathode.scalar[xIndex][yIndex][zIndex] << endl;
		}
		else
			cout << "Invalid index values." << endl;
		cout << "Would you like to exit? (Enter 0 for no and 1 for yes): ";
		cin >> toExit;
	}
}
