#ifndef SIMULATION_H
#define SIMULATION_H

#include "Detector.h"
#include "MyRandom.h"
#include "phiLogStruct.h"
#include "trailLogStruct.h"
#include <list>
#include <time.h>
#include <vector>

using std::list;
using std::vector;

class CTSI;
class Event;
class MyRandom;
class Preamplification;

//#define Q 1							// Debug
#define Q 1.60217646e-19				// Electron charge in C
#define EPSILON_0 8.85418782e-12        // Permitivity of free space in F/m
#define K 1.3806503e-23                 // Boltzmann's constant in J/K
#define SQRT3_INVERSE 0.57735026918963  // 1/sqrt(3)

class Simulation
{
	public:
		Simulation();
		~Simulation();
		int				engine(void* ctsiObjectPtr, int indent, void (*dispFuncPtr)(void*, int), void (*outputFuncPtr)(void*, int, char*, int &,
																													   vector<int> &, vector<double> &, vector<double> &, vector<double> &, vector<double> &,
																													   vector<int> &, vector<double> &, vector<double> &, vector<double> &, vector<double> &,
																													   vector<int> &, vector<int> &,
																													   vector<double> &, vector<vector<double> > &, vector<vector<double> > &,
																													   vector<vector<struct trailLog> > &, vector<vector<vector<struct phiLog> > > &, vector<vector<vector<struct phiLog> > > &,
																													   vector<vector<int> > &, vector<vector<int> > &,
																													   vector<vector<double> > &, vector<vector<double> > &,
																													   vector<vector<double> > &, vector<vector<double> > &,
																													   vector<double> &, vector<vector<double> > &, vector<vector<double> > &,
																													   vector<vector<double> > &, vector<vector<double> > &,
																													   bool),
							   struct configParamStruct &configParams, Detector *detector, Preamplification *preamplifier, list<Event> &eventArray, char* runMode);

	protected:
		// Methods
		double						avgVec(vector<double> myVec);
		void						calcChargeChange(struct trailLog &chargeTrailsIn, int &intrxnIndex, int &currSimStep, double &dS, int chargeType);
		void						calcCloudSize(int &numChargeElems, int &intrxnIndex, int &currSimStep);
		void						calcNewChargePos(struct trailLog &chargeTrailsIn, double &dX, double &dY, double &dZ, int chargeType);
		void						clearPhiLog();
		void						clearTrailLog();
		bool						eIsCloseToAnode();
		void						flushVectors();
		double						getChargeElemDisp(struct configParamStruct &configParams, Detector *detector, struct trailLog const &chargeTrailsIn, double &dX, double &dY, double &dZ, int const chargeType);
		double						getDeltaSigmaSq(double sigmaIn, double chargeInCloud, int chargeType);
		void						getTrailPhi(struct configParamStruct &configParams, Detector *detector, int &maxLen, vector<vector<int> > &trailSizeE, vector<vector<int> > &trailSizeH, bool caution);
		void						getTriggers(double &anodeThresh, double &cathodeThresh);
		void						getDetectorOutput(Detector *detector, int &numChargeElems, int &maxLen, vector<vector<int> > &trailSizeE, vector<vector<int> > &trailSizeH);
		void						getPreampOutput(struct configParamStruct &configParams, Preamplification *preamplifier);
		double						getSigma(vector<struct trailLog> &interxnTrails, vector<int> &reachedEnd, int &numChargeElems, int index, int chargeType);
		void						getTimeIncr();
		void						initSim4Event(struct configParamStruct &configParams, Event &eventIn);
		void						initSim4Interxn(Event &eventIn, int intrxnIndex, int numChargeElems, int maxSimSteps);
		vector<double>				interp1(vector<double> &xIrreg, vector<double> &yIrreg, vector<double> &xReg);
		float						listGet(list<float> myList, int index);
		void						setDetectorSpecs(struct detectorSpecStruct &detectorSpecsIn);
		void						setSeed(int seed);
		void						terminateTrail(struct trailLog &chargeTrailsIn, int &intrxnIndex, int &chargeElem, int &currSimStep, int chargeType);

		// Attributes
													// ===== Assign-once items =====
		struct detectorSpecStruct	detectorSpecs;	// Structure for detector specification parameters.
		MyRandom*					randGen;		// Pointer to random number generator.
		struct phiLog				chargePhis;		// Dummy data object storing the phi profile of one charge element.
		struct trailLog				chargeTrails;	// Dummy data object storing the charge and position profile of one charge element.
		vector<int>					zeroVec;		// Dummy vector used in initialization.

													// ===== Core data structures =====
		vector<vector<struct trailLog> >			eventTrails;	// Trail log data object
		vector<vector<vector<struct phiLog> > >		eventAnPhi;		// Anode phi Log data object
		vector<vector<vector<struct phiLog> > >		eventCaPhi;		// Cathode phi Log data object

													// ===== Output data structures =====
		vector<int>					anCollectorOneHot;	// Vector for keeping track of which anode has collected charge.
		vector<int>					caCollectorOneHot;	// Vector for keeping track of which anode has collected charge.
		vector<int>					anCollectorID;		// Vector of 1 positions in anCollectorOneHot.
		vector<int>					caCollectorID;		// Vector of 1 positions in caCollectorOneHot.
		vector<double>				timeVec;			// Time as a function of simulation step.		
		vector<vector<double> >		anVsTime;			// Total induced anode signal as a function of simulation step, refreshed for each event.
		vector<vector<double> >		caVsTime;			// Total induced cathode signal as a function of simulation step, refreshed for each event.
		vector<vector<double> >		anVsTimeReg;		// anVsTime with anode values at regular sampling intervals.
		vector<vector<double> >		caVsTimeReg;		// caVsTime with cathode values at regular sampling intervals.
		vector<vector<double> >		noisyAnVsTimeReg;	// Noisy version of anVsTime with anode values at regular sampling intervals.
		vector<vector<double> >		noisyCaVsTimeReg;	// Noisy version of caVsTime with cathode values at regular sampling intervals.
		vector<double>				timeVecPreamp;		// timeVec with time values at regular sampling intervals.
		vector<vector<double> >		anVsTimePreamp;		// Preamplified version of anVsTime.
		vector<vector<double> >		caVsTimePreamp;		// Preamplified version of caVsTime.
		vector<vector<double> >		noisyAnVsTimePreamp;// Preamplified version of noisyAnVsTime.
		vector<vector<double> >		noisyCaVsTimePreamp;// Preamplified version of noisyCaVsTime.
		vector<int>					anodeID;			// Vector of IDs of anodes that have triggered.
		vector<int>					cathodeID;			// Vector of IDs of cathodes that have triggered.
		vector<double>				anEnergy;			// Vector of signal amplitudes on triggered anodes.
		vector<double>				caEnergy;			// Vector of signal amplitudes on triggered cathodes.
		vector<double>				noisyAnEnergy;		// Vector of noisy signal amplitudes on triggered anodes.
		vector<double>				noisyCaEnergy;		// Vector of noisy signal amplitudes on triggered cathodes.
		vector<double>				anTriggerTime;		// Vector storing the time at which anode triggers occurred.
		vector<double>				caTriggerTime;		// Vector storing the time at which cathode triggers occurred.
		vector<double>				noisyAnTriggerTime;	// Vector storing the time at which noisy anode triggered.
		vector<double>				noisyCaTriggerTime;	// Vector storing the time at which noisy cathode triggered.

													// ===== Time resolution selection =====
		int							actNumChargeElem; // Actual number of charge elements can differ from the minimum specified
		double						dt;				// Current time increment.
		double						targetDt;		// Target time increment.
		double						ddt;			// Change to apply to dt;
		vector<vector<int> >		eCloseToAnode;	// Vector of vector of flags signalling if electrons are beyond the small
													// pixel distance to switch to a finer dt. 1 denotes true and 0 denotes
													// false, the flag is reset to 0 once an electron has reached the anode
													// plane.

													// ===== Cloud size calculation =====
		double						sigmaInit;		// Initial sigma of charge carrier cloud size in cm.
		double						deltaSigmaSqE;	// The sigma-square change in electron cloud size.
		double						deltaSigmaSqH;	// The sigma-square change in hole cloud size.
		vector<vector<double> >		sigmaSqE;		// Log of cloud size in terms of electron charge element variance.
													// The last entry is the expected theoretical variance, all other
													// entries are sample variances. This log is not actually used.
		vector<vector<double> >		sigmaSqH;		// Log of cloud size in terms of hole charge element variance.
													// The last entry is the expected theoretical variance, all other
													// entries are sample variances. This log is not actually used.
		vector<double>				sigmaSqInit;	// Place holder for sigmaInit^2.
													
													// Terminating condition
		vector<vector<int> >		reachedEndE;	// A matrix of flags indicating if each electron element has been collected. 
		vector<int>					reachedEndESum;	// Each element is the total number of electron elements already collected.
													// for the i-th interaction in an event.
		vector<vector<int> >		reachedEndH;	// A matrix of flags indicating if each hole element has been collected.
		vector<int>					reachedEndHSum;	// Each element is the total number of hole elements already collected.
													// for the i-th interaction in an event.
		vector<double>				chargeInCloudE; // Total remaining mobile electron charge in cloud in current simulation step.
		vector<double>				chargeInCloudH; // Total remaining mobile hole charge in cloud in current simulation step.
		int							numReachedEndE;	// Number of electrons already at anode plane at current simulation step.
		int							numReachedEndH;	// Number of holes already at cathode plane at current simulation step.
		int							maxReachedEndE;	// Total number of electron elements there are in one event, across all its interactions.
		int							maxReachedEnd;	// Total number of charge elements there are in one event, across all its interactions.
};

#endif // SIMULATION_H
