#include "Event.h"
#include "MyRandom.h"
#include "Preamplification.h"
#include "Simulation.h"
#include <cstring>
#include <fstream>

using namespace std;

#define UPDATE_INCREMENT 10

//===================================================================
// Constructor
Simulation::Simulation()
{
	dt = 1e-9;					// Changing dt value for different circumstances in simulation
	targetDt = dt;
	ddt = 0;
	randGen = new MyRandom();
	sigmaInit = 1e-3;
	sigmaSqInit.push_back(sigmaInit*sigmaInit);

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

	chargePhis.EPhi.clear();
	chargePhis.HPhi.clear();

	chargeTrails.xPosE.clear();
	chargeTrails.yPosE.clear();
	chargeTrails.zPosE.clear();
	chargeTrails.qETrapped.clear();
	chargeTrails.qEMobile.clear();
	chargeTrails.xPosH.clear();
	chargeTrails.yPosH.clear();
	chargeTrails.zPosH.clear();
	chargeTrails.qHTrapped.clear();
	chargeTrails.qHMobile.clear();

	clearTrailLog();
	flushVectors();

	numReachedEndE = 0;
	numReachedEndH = 0;
	maxReachedEndE = 0;
	maxReachedEnd = 0;

	deltaSigmaSqE = 0;
	deltaSigmaSqH = 0;

	actNumChargeElem = 0;
}

Simulation::~Simulation()
{
	delete randGen;
}

//===================================================================
// This function returns the average of a vector of doubles
double Simulation::avgVec(vector<double> myVec)
{
	int n = int (myVec.size());
	double sum = 0;
	for(int i = 0; i < n; i++)
		sum += myVec[i];
	return sum/n;
}

//===================================================================
// This function updates the mobile and trapped charge values
void Simulation::calcChargeChange(struct trailLog &chargeTrailsIn, int &intrxnIndex, int &currSimStep, double &dS, int chargeType)
{
	double mobileCharge = 0;
	double trappedCharge = 0;

	// Electron
	if (chargeType < 1)
	{
		// Calculate the amount of mobile charge remaining in current simulation step.
		mobileCharge = chargeTrailsIn.qEMobile.back()*exp(-dS/detectorSpecs.LAMBDA_E);
		trappedCharge = chargeTrailsIn.qEMobile.back() - mobileCharge;

		chargeTrailsIn.qEMobile.push_back(mobileCharge);
		chargeTrailsIn.qETrapped.back() = trappedCharge;
		chargeTrailsIn.qETrapped.push_back(0);
		
		// Update mobile electron charge in cloud used to calculate
		// cloud expansion. Subtract off trapped charges.
		chargeInCloudE[intrxnIndex] -= trappedCharge;
	}
	// Holes
	else
	{
		// Update trapped charges
		mobileCharge = chargeTrailsIn.qHMobile.back()*exp(-dS/detectorSpecs.LAMBDA_H);
		trappedCharge = chargeTrailsIn.qHMobile.back() - mobileCharge;

		chargeTrailsIn.qHMobile.push_back(mobileCharge);
		chargeTrailsIn.qHTrapped.back() = trappedCharge;
		chargeTrailsIn.qHTrapped.push_back(0);

		// Update mobile hole charge in cloud used to calculate
		// cloud expansion. Subtract off trapped charges.
		chargeInCloudH[intrxnIndex] -= trappedCharge;
	}
}

//===================================================================
// This function calculates the sample charge cloud size for the
// current simulation step, and the theoretical charge cloud size for
// the next simulation step.
void Simulation::calcCloudSize(int &numChargeElems, int &intrxnIndex, int &currSimStep)
{
	// Use actual sigma^2 value to determine the next value
	double oldSigmaE;
	if (reachedEndESum[intrxnIndex] != numChargeElems) {
		oldSigmaE = getSigma(eventTrails[intrxnIndex], reachedEndE[intrxnIndex], numChargeElems, currSimStep - 1, -1);
		deltaSigmaSqE = getDeltaSigmaSq(oldSigmaE, chargeInCloudE[intrxnIndex], -1);
		// Remember sigmaSqE is not used anywhere else, this is just for record keeping
		// Okay to comment out the next 2 lines for faster run speed
		// sigmaSqE[intrxnIndex][currSimStep - 1] = oldSigmaE*oldSigmaE;
		// sigmaSqE[intrxnIndex].push_back(sigmaSqE[intrxnIndex].back() + deltaSigmaSqE);
	}

	// Use actual sigma^2 value to determine the next value
	double oldSigmaH;
	if (reachedEndHSum[intrxnIndex] != numChargeElems) {
		oldSigmaH = getSigma(eventTrails[intrxnIndex], reachedEndH[intrxnIndex], numChargeElems, currSimStep - 1, 1);
		deltaSigmaSqH = getDeltaSigmaSq(oldSigmaH, chargeInCloudH[intrxnIndex], 1);
		// Remember sigmaSqH is not used anywhere else, this is just for record keeping
		// Okay to comment out the next 2 lines for faster run speed
		// sigmaSqH[intrxnIndex][currSimStep - 1] = oldSigmaH*oldSigmaH;
		// sigmaSqH[intrxnIndex].push_back(sigmaSqH[intrxnIndex].back() + deltaSigmaSqH);
	}
}

//===================================================================
// This function calculates the charge element positions for the
// current simulation step
void Simulation::calcNewChargePos(struct trailLog &chargeTrailsIn, double &dX, double &dY, double &dZ, int chargeType)
{
	// Electron
	if (chargeType < 1)
	{
		chargeTrailsIn.xPosE.push_back(chargeTrailsIn.xPosE.back() + dX);
		chargeTrailsIn.yPosE.push_back(chargeTrailsIn.yPosE.back() + dY);
		chargeTrailsIn.zPosE.push_back(chargeTrailsIn.zPosE.back() + dZ);
	}
	// Holes
	else
	{
		chargeTrailsIn.xPosH.push_back(chargeTrailsIn.xPosH.back() + dX);
		chargeTrailsIn.yPosH.push_back(chargeTrailsIn.yPosH.back() + dY);
		chargeTrailsIn.zPosH.push_back(chargeTrailsIn.zPosH.back() + dZ);
	}
}

//===================================================================
// This functions clears all the fields of the phiLog structure.
void Simulation::clearPhiLog()
{
	chargePhis.EPhi.clear();
	chargePhis.HPhi.clear();
}

//===================================================================
// This functions clears all the fields of the trailLog structure.
void Simulation::clearTrailLog()
{
	chargeTrails.qETrapped.clear();
	chargeTrails.qEMobile.clear();
	chargeTrails.qHTrapped.clear();
	chargeTrails.qHMobile.clear();
	chargeTrails.xPosE.clear();
	chargeTrails.xPosH.clear();
	chargeTrails.yPosE.clear();
	chargeTrails.yPosH.clear();
	chargeTrails.zPosE.clear();
	chargeTrails.zPosH.clear();
}

//===================================================================
// This function returns true if any electron charge element is close
// to the anode plane.
bool Simulation::eIsCloseToAnode()
{
	// Loop over all interactions
	for (int i = 0; i < (int) eCloseToAnode.size(); i++)
		// Loop over all charge elements
		for(int j = 0; j < (int) eCloseToAnode[i].size(); j++)
			if (eCloseToAnode[i][j] == 1)
				return true;

	return false;
}

//===================================================================
// This function performs the actual physics simulation of charge
// transport and signal induction (CTSI) inside the detector for the
// events in the event array.
int Simulation::engine(void*		ctsiObjectPtr,
					   int			indent,
					   void			(*dispFuncPtr)(void*, int),
					   void			(*outputFuncPtr)(void*, int, char*, int &,
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
					   struct configParamStruct		&configParams,
					   Detector*	detector,
					   Preamplification	*preamplifier,
					   list<Event>	&eventArray,
					   char*	runMode)
{
	(*dispFuncPtr)(ctsiObjectPtr, indent);
	cout << "Simulating charge transport and signal induction:" << endl;

	// Set random seed for the engine
	setSeed(configParams.RANDOM_SEED);

	// Copy the detector's specifications into the data member so
	// we don't have to refer to the detector each time we use the
	// specification values.
	setDetectorSpecs(detector->getDetectorSpecs());
	int maxLen = 0;
	vector<vector<int> > trailSizeE, trailSizeH;

	// Vectors to store anode and cathode values
	vector<double> An, Ca;

	// List iterator
	int eventCount = 0;
	int skippedCount = 0;
	list<Event>::iterator it_listEvent;

	// Setting triggering threshold
	double anodeThresh = configParams.ANODE_TRIG_THRESH;
	double cathodeThresh = configParams.CATHODE_TRIG_THRESH;

	int numIntrxn = 0;					// Number of interactions in an event.
	double tempTime = 0, timePerEvent = 0;
	bool terminatedEarly = false;
	int numEvents = (int) eventArray.size();
	int currSimStep = 1;
	int lastEventCount = eventCount;
	char modeStr[256];
	modeStr[0] = 'i';
	modeStr[1] = '\0';

	time_t secSinceUpdate, startTime;
	startTime = time(NULL);
	secSinceUpdate = startTime;

	try
	{
		//=========================================================================================
		// ITERATE OVER ALL EVENTS
		//=========================================================================================
		for (it_listEvent = eventArray.begin(); it_listEvent != eventArray.end(); it_listEvent++)
		{
			// Console update
			if (time(NULL) - secSinceUpdate >= UPDATE_INCREMENT)
			{
				tempTime  = (double) (time(NULL) - startTime);
				timePerEvent = tempTime/((double) eventCount);

				(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
				cout << "Info: " << eventCount - skippedCount << " events processed, " << skippedCount << " events skipped, " << numEvents - eventCount << " events remaining." << endl;
				(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
				cout << "Info: ";
				if (tempTime > 86400) {
					int numDays = (int) tempTime/86400;
					cout << numDays << " days ";
					tempTime -= 86400*numDays;
				}
                if (tempTime > 3600) {
					int numHours = (int) tempTime/3600;
					cout << numHours << " hrs ";
					tempTime -= 3600*numHours;
				}
				if (tempTime > 60) {
					int numMinutes = (int) tempTime/60;
					cout << numMinutes << " minutes ";
					tempTime-= 60*numMinutes;
				}
				cout << tempTime << " seconds elapsed, ";
				
				tempTime = timePerEvent*(numEvents - eventCount);
				if (tempTime > 86400) {
					int numDays = (int) tempTime/86400;
					cout << numDays << " days ";
					tempTime -= 86400*numDays;
				}
                if (tempTime > 3600) {
					int numHours = (int) tempTime/3600;
					cout << numHours << " hrs ";
					tempTime -= 3600*numHours;
				}
				if (tempTime > 60) {
					int numMinutes = (int) tempTime/60;
					cout << numMinutes << " minutes ";
					tempTime-= 60*numMinutes;
				}
				cout << tempTime << " seconds remaining at " << 1/timePerEvent << " events/sec (" << (double) (eventCount - lastEventCount)/(double) UPDATE_INCREMENT << " events/sec based on the last " << eventCount - lastEventCount << " events)." << endl;
                lastEventCount = eventCount;
				secSinceUpdate = time(NULL);
			}

			// Vary the number of charge elements to smear out the effects of charge
			// discretization
			actNumChargeElem = configParams.NUM_CHARGE_ELEMS + (int) (randGen->uniform()*10);

			numIntrxn = (int) (*it_listEvent).getNumIntrxns();

			flushVectors();				// eventTrails refreshed here
			initSim4Event(configParams, (*it_listEvent));
			timeVec.push_back(0);		// Set initial time to 0
			currSimStep = 1;			// Index for tracking point along trajectories
			terminatedEarly = false;

			//=====================================================================================
			// SIMULATE CHARGE TRANSPORT
			// Calculate the trajectories of charge elements
			//=====================================================================================
			try
			{
				// Execute while there are still charge elements in transit
				while(numReachedEndE + numReachedEndH < maxReachedEnd)
				{
					// Check against maximum permitted number of simulation steps
					if (currSimStep > configParams.MAX_NUM_SIM_STEPS)
					{
						(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
						cout << "Info: event " << eventCount << " failed to terminate in " << configParams.MAX_NUM_SIM_STEPS <<" steps, skipped." << endl;
						skippedCount++;
						terminatedEarly = true;
						break;
					}

					// Choose the correct dt value
					getTimeIncr();

					//=================================================================================
					// ITERATE OVER EACH INTERACTION
					//=================================================================================
					int intrxnIndex = 0;
					try
					{
						for(intrxnIndex = 0; intrxnIndex < numIntrxn; intrxnIndex++)
						{
							// Find sample variance
							// This is the only place where deltaSigmaSqE and deltaSigmaSqH are assigned
							calcCloudSize(actNumChargeElem, intrxnIndex, currSimStep);

							//=============================================================================
							// ITERATE OVER EACH CHARGE ELEMENT
							//=============================================================================
							int chargeElem = 0;
							double dX = 0, dY = 0, dZ = 0, dS = 0;
							try
							{
								for (chargeElem = 0; chargeElem < actNumChargeElem; chargeElem++)
								{
									//=====================================================================
									// Simulate mobile electrons
									if (reachedEndE[intrxnIndex][chargeElem] == 0)
									{
										// 1. Calculate charge displacement
										// 2. Update mobile and trapped charge population
										// 3. Update x, y, z positions
										// 4. Turn on signal to change to finer dt when nearing boundary

										/* 1. */ dS = getChargeElemDisp(configParams, detector, eventTrails[intrxnIndex][chargeElem], dX, dY, dZ, -1);
										/* 2. */ calcChargeChange(eventTrails[intrxnIndex][chargeElem], intrxnIndex, currSimStep, dS, -1);
										/* 3. */ calcNewChargePos(eventTrails[intrxnIndex][chargeElem], dX, dY, dZ, -1);
										/* 4. */ if(eventTrails[intrxnIndex][chargeElem].zPosE.back() > configParams.SML_PIXEL_DIST)
													eCloseToAnode[intrxnIndex][chargeElem] = 1;

										// Check if anode plane is reached
										if (eventTrails[intrxnIndex][chargeElem].zPosE.back() < detectorSpecs.L)
										{
											// Correct for boundary cases in x-direction
											// We assume that if the electron charge element hits the x or y boudanries,
											// it would continue drifting toward the anode plane, so z component of
											// position needs no adjustment and the trail needs no terminaiton.
											// 
											// NOTE: since we are doing boundary check after the charge update, the
											// charge update would be inconsistent with the position update. This may or
											// may not be significant depending on the magnitude of the correction.
											double curX = eventTrails[intrxnIndex][chargeElem].xPosE.back();
											double curY = eventTrails[intrxnIndex][chargeElem].yPosE.back();
											if (curX < -detectorSpecs.W/2)
												eventTrails[intrxnIndex][chargeElem].xPosE[currSimStep] = -detectorSpecs.W/2;
											if (curX > detectorSpecs.W/2)
												eventTrails[intrxnIndex][chargeElem].xPosE[currSimStep] = detectorSpecs.W/2;
											if (curY < -detectorSpecs.W/2)
												eventTrails[intrxnIndex][chargeElem].yPosE[currSimStep] = -detectorSpecs.W/2;
											if (curY > detectorSpecs.W/2)
												eventTrails[intrxnIndex][chargeElem].yPosE[currSimStep] = detectorSpecs.W/2;
										}
										else
										{
											// Catch if electrons drifted to the cathode plane
											if (eventTrails[intrxnIndex][chargeElem].zPosE.back() < 0)
												throw 5;
											else
											{
												// Modify variables to reflect termination of charge element trail.
												terminateTrail(eventTrails[intrxnIndex][chargeElem], intrxnIndex, chargeElem, currSimStep, -1);

												// Find out if charge collected by anode, and if so, which anode
												detector->anodeCollection(eventTrails[intrxnIndex][chargeElem], currSimStep, anCollectorOneHot, configParams.CAUTION);
											}
										}
									}

									//=====================================================================
									// Simulate mobile holes
									if (reachedEndH[intrxnIndex][chargeElem] == 0)
									{
										// 1. Calculate charge displacement
										// 2. Update mobile and trapped charge population
										// 3. Update x, y, z positions

										/* 1. */ dS = getChargeElemDisp(configParams, detector, eventTrails[intrxnIndex][chargeElem], dX, dY, dZ, 1);
										/* 2. */ calcChargeChange(eventTrails[intrxnIndex][chargeElem], intrxnIndex, currSimStep, dS, 1);
										/* 3. */ calcNewChargePos(eventTrails[intrxnIndex][chargeElem], dX, dY, dZ, 1);

										// Check if anode plane is reached
										if(eventTrails[intrxnIndex][chargeElem].zPosH.back() >= 0)
										{
											// Correct for boundary cases in x-direction
											// We assume that if the hole charge element hits the x or y boudanries,
											// it would continue drifting toward the cathode plane, so z component of
											// position needs no adjustment and the trail needs no terminaiton.
											double curX = eventTrails[intrxnIndex][chargeElem].xPosH.back();
											double curY = eventTrails[intrxnIndex][chargeElem].yPosH.back();
											if (curX < -detectorSpecs.W/2)
												eventTrails[intrxnIndex][chargeElem].xPosH[currSimStep] = -detectorSpecs.W/2;
											if (curX > detectorSpecs.W/2)
												eventTrails[intrxnIndex][chargeElem].xPosH[currSimStep] = detectorSpecs.W/2;
											if (curY < -detectorSpecs.W/2)
												eventTrails[intrxnIndex][chargeElem].yPosH[currSimStep] = -detectorSpecs.W/2;
											if (curY > detectorSpecs.W/2)
												eventTrails[intrxnIndex][chargeElem].yPosH[currSimStep] = detectorSpecs.W/2;
										}
										else
										{
											// Catch if holes drifted to the anode plane
											if (eventTrails[intrxnIndex][chargeElem].zPosH.back() > detectorSpecs.L)
												throw 5;
											else
											{
												// Modify variables to reflect termination of charge element trail.
												terminateTrail(eventTrails[intrxnIndex][chargeElem], intrxnIndex, chargeElem, currSimStep, 1);

												// Find out if charge collected by cathode, and if so, which cathode
												detector->cathodeCollection(eventTrails[intrxnIndex][chargeElem], currSimStep, caCollectorOneHot, configParams.CAUTION);
											}
										}
									}
								} // End of for loop over charge elements
							} // End of try
							catch (int code)
							{
                                (*dispFuncPtr)(ctsiObjectPtr, indent + 1);
								cerr << "Error: charge element " << chargeElem << "." << endl;

								switch (code)
								{
									case 1:
										cerr << "       calculated triggered anode index out of range." << endl;
										throw 0;
										break;
									case 2:
										(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
										cerr << "       electron collection x pos. is more than half a pitch distance away from assgined anode axis." << endl;
										throw 0;
										break;
									case 3:
										(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
										cerr << "       calculated triggered cathode index out of range." << endl;
										throw 0;
										break;
									case 4:
										(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
										cerr << "       hole collection y pos. is more than half a pitch distance away from assigned cathode axis." << endl;
										throw 0;
										break;
									case 5:
										(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
										cerr << "       found electrons at cathode plane or holes at anode plane." << endl;
										throw 0;
										break;
									default:
										throw 0;
								}
							}
						} // End of for loop over interactions
					}
					catch (...)
					{
						(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
						cerr << "Error: interaction " << intrxnIndex << "." << endl;
						throw 0;
					}

					timeVec.push_back(timeVec.back() + dt);
					currSimStep++;

				} // End of while
			} // End of try
			catch (...)
			{
				(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
				cerr << "Error: simulation step " << currSimStep << "." << endl;
				throw 0;
			}
			
			//=====================================================================================
			// SIMULATE SIGNAL INDUCTION
			//=====================================================================================
			if (!terminatedEarly)
			{
				try
				{
					// Expand collector set to include neighboring electrodes
					// within configurable window for neighbor signal analysis.
					int window = configParams.NEIGHBOR_WINDOW;
					if (window < 0) {
						// -1 = compute ALL electrodes (legacy behavior)
						for (int i = 0; i < (int) anCollectorOneHot.size(); i++)
							anCollectorOneHot[i] = 1;
					} else {
						// Find range of collecting anodes, expand by +/-window
						int anMin = (int) anCollectorOneHot.size();
						int anMax = -1;
						for (int i = 0; i < (int) anCollectorOneHot.size(); i++) {
							if (anCollectorOneHot[i] == 1) {
								if (i < anMin) anMin = i;
								if (i > anMax) anMax = i;
							}
						}
						if (anMax >= 0) {
							int lo = (anMin - window > 0) ? anMin - window : 0;
							int hi = (anMax + window < (int) anCollectorOneHot.size() - 1)
							         ? anMax + window : (int) anCollectorOneHot.size() - 1;
							for (int i = lo; i <= hi; i++)
								anCollectorOneHot[i] = 1;
						}
					}
					// Cathodes: always compute all (only 8 strips, negligible cost)
					for (int i = 0; i < (int) caCollectorOneHot.size(); i++)
						caCollectorOneHot[i] = 1;

					// Retrieve the phi profile and get list-mode output
					getTrailPhi(configParams, detector, maxLen, trailSizeE, trailSizeH, configParams.CAUTION);

					// Get detector output as a function of time
					getDetectorOutput(detector, actNumChargeElem, maxLen, trailSizeE, trailSizeH);

					// Get preamplifier output as a function of time
					getPreampOutput(configParams, preamplifier);

					// Get triggers
					getTriggers(anodeThresh, cathodeThresh);
				}
				catch (int code)
				{
					(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
					switch(code) 
					{
						case 0:
							cerr << "Error: Vectors of the same charge element have different lengths." << endl;
							throw 0;
						case 1:
							cerr << "Error: Charge trail and time vector has different lengths." << endl;
							throw 0;
							break;
						case 2:
							cerr << "Error: electron collection x pos. is more than half a pitch distance away from anode phi axis." << endl;
							throw 0;
							break;
						case 3:
							cerr << "Error: hole collection y pos. is more than half a pitch distance away from cathode phi axis." << endl;
							throw 0;
							break;
						default:
							(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
							cerr << "Error occurred while calculating output." << endl;
							throw 0;
					}
				}

				// Pass simulation results to output function
				// Output only if at least one electrode was triggered or if we're in interactive mode.
				if ((anodeID.size() > 0) || (cathodeID.size() > 0) || !strcmp(runMode, "i"))
					(*outputFuncPtr)(ctsiObjectPtr, indent, runMode, eventCount,
									 anodeID, anEnergy, anTriggerTime, noisyAnEnergy, noisyAnTriggerTime,
									 cathodeID, caEnergy, caTriggerTime, noisyCaEnergy, noisyCaTriggerTime,
									 anCollectorID, caCollectorID,
									 timeVec, anVsTime, caVsTime,
									 eventTrails, eventAnPhi, eventCaPhi,
									 trailSizeE, trailSizeH,
									 anVsTimeReg, caVsTimeReg,
									 noisyAnVsTimeReg, noisyCaVsTimeReg,
									 timeVecPreamp, anVsTimePreamp, caVsTimePreamp,
									 noisyAnVsTimePreamp, noisyCaVsTimePreamp,
									 false);
			}

			eventCount++;

			if (eventCount == numEvents)
			{
				// Call the output function one last time to signal there are no more events
				(*outputFuncPtr)(ctsiObjectPtr, indent, runMode, eventCount,
									anodeID, anEnergy, anTriggerTime, noisyAnEnergy, noisyAnTriggerTime,
									cathodeID, caEnergy, caTriggerTime, noisyCaEnergy, noisyCaTriggerTime,
									anCollectorID, caCollectorID,
									timeVec, anVsTime, caVsTime,
									eventTrails, eventAnPhi, eventCaPhi,
									trailSizeE, trailSizeH,
									anVsTimeReg, caVsTimeReg,
									noisyAnVsTimeReg, noisyCaVsTimeReg,
									timeVecPreamp, anVsTimePreamp, caVsTimePreamp,
									noisyAnVsTimePreamp, noisyCaVsTimePreamp,
									true);
			}
		} // end of for loop through events
	}
	catch (int)
	{
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cerr << "Error: event " << eventCount << endl;
		(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
		cerr << "Simulation aborted." << endl;
		throw 0;
	}

	(*dispFuncPtr)(ctsiObjectPtr, indent + 1);
	cout << "Done." << endl;

	return 0;
}

//===================================================================
// This function clears vector and structure member variables.
void Simulation::flushVectors()
{
	sigmaSqE.clear();
	sigmaSqH.clear();

	eventTrails.clear();
	eventAnPhi.clear();
	eventCaPhi.clear();

	anCollectorOneHot.clear();
	caCollectorOneHot.clear();
	anCollectorID.clear();
	caCollectorID.clear();

	// clearPhiLog();
	// clearTrailLog();

	chargeInCloudE.clear();
	chargeInCloudH.clear();

	reachedEndE.clear();
	reachedEndH.clear();

	reachedEndESum.clear();
	reachedEndHSum.clear();

	eCloseToAnode.clear();

	timeVec.clear();
	anVsTime.clear();
	caVsTime.clear();

	timeVecPreamp.clear();
	anVsTimeReg.clear();
	caVsTimeReg.clear();
	noisyAnVsTimeReg.clear();
	noisyCaVsTimeReg.clear();
	anVsTimePreamp.clear();
	caVsTimePreamp.clear();
	noisyAnVsTimePreamp.clear();
	noisyCaVsTimePreamp.clear();

	anodeID.clear();
	cathodeID.clear();
	anEnergy.clear();
	caEnergy.clear();
	noisyAnEnergy.clear();
	noisyCaEnergy.clear();
	anTriggerTime.clear();
	caTriggerTime.clear();
	noisyAnTriggerTime.clear();
	noisyCaTriggerTime.clear();

	zeroVec.clear();
}

//===================================================================
// This function calculates the amount of displacement of an charge
// element given the electric field and diffusion parameter.
double Simulation::getChargeElemDisp(struct configParamStruct &configParams, Detector *detector, struct trailLog const &chargeTrailsIn, double &dX, double &dY, double &dZ, int const chargeType)
{
	// Electron
	if (chargeType < 1)
	{
		// Interpolate electric field using given values at grid points
		vector<double> eFieldE;
		detector->getEField(configParams, chargeTrailsIn.xPosE.back(), chargeTrailsIn.zPosE.back(), eFieldE);
		double eFieldXE = eFieldE[0];
		double eFieldZE = eFieldE[1];

		// Calculate change in position due to both drift and diffusion
		dX = -detectorSpecs.MU_E*eFieldXE*dt + randGen->normal()*sqrt(deltaSigmaSqE)*SQRT3_INVERSE;
		dY =  randGen->normal()*sqrt(deltaSigmaSqE)*SQRT3_INVERSE;
		dZ = -detectorSpecs.MU_E*eFieldZE*dt + randGen->normal()*sqrt(deltaSigmaSqE)*SQRT3_INVERSE;

		// Debug
		//double randNum = randGen->Gaus()*sqrt(deltaSigmaSqE);
		//cout << randNum << endl;
		//dX = -detectorSpecs.MU_E*eFieldXE*dt + randNum*SQRT3_INVERSE;

		//randNum = randGen->Gaus()*sqrt(deltaSigmaSqE);
		//cout << randNum << endl;
		//dY = randNum*SQRT3_INVERSE;
		//
		//randNum = randGen->Gaus()*sqrt(deltaSigmaSqE);
		//cout << randNum << endl;
		//dZ = -detectorSpecs.MU_E*eFieldZE*dt + randNum*SQRT3_INVERSE;
	}
	// Hole
	else
	{
		// Interpolate electric field using given values at grid points
		vector<double> eFieldH;
		detector->getEField(configParams, chargeTrailsIn.xPosH.back(), chargeTrailsIn.zPosH.back(), eFieldH);
		double eFieldXH = eFieldH[0];
		double eFieldZH = eFieldH[1];

		// Calculate change in position due to both drift and diffusion
		dX = detectorSpecs.MU_H*eFieldXH*dt + randGen->normal()*sqrt(deltaSigmaSqH)*SQRT3_INVERSE;
		dY = randGen->normal()*sqrt(deltaSigmaSqH)*SQRT3_INVERSE;
		dZ = detectorSpecs.MU_H*eFieldZH*dt + randGen->normal()*sqrt(deltaSigmaSqH)*SQRT3_INVERSE;
	}

	return sqrt(dX*dX + dY*dY + dZ*dZ);
}

//===================================================================
// This function calculates the change in sigma^2 of the Gaussian charge
// cloud model from one time step to the next separated by dt seconds. The
// formulation is based on the charge cloud expansion model as described in
//
// Inputs:
// sigmaIn: sigma of the current Gaussian cloud model.
// chargeInCloud: amount of charge in cloud in coulombs.
// dt: the time duration over which the change in sigma is sought.
// chargeType: 1 for holes and -1 for electrons.
//
// M. Benoit and L.A. Hamel, 'Simulation of charge collection prcesses in
// semiconductor CdZnTe gamma-ray detectors', Nuclear Instruments and
// Methods in Physics Research A, 2009.
//
// The function assumes the following detector properties and conditions:
// Material: CZT
// Electron mobility: 1000 cm^2/Vs
// Hole mobility: 50 cm^2/Vs
// Temperature: 293 K (20 degrees C)
// Dielectric constant: 10.9
double Simulation::getDeltaSigmaSq(double sigmaIn, double chargeInCloud, int chargeType)
{
	if(sigmaIn == 0)
		return 0;
	double Dprime;
	if (chargeType > 0)
		Dprime = detectorSpecs.DIFFUSION_H + detectorSpecs.MU_H*abs(chargeInCloud)/(24*pow(PI,1.5)*(EPSILON_0/100)*detectorSpecs.EPSILON_R*sigmaIn);
		// Debug
		// Dprime = detectorSpecs.DIFFUSION_H + detectorSpecs.MU_H*abs(chargeInCloud*1.60217646e-19)/(24*pow(PI,1.5)*(EPSILON_0/100)*detectorSpecs.EPSILON_R*sigmaIn);
	else
		Dprime = detectorSpecs.DIFFUSION_E + detectorSpecs.MU_E*abs(chargeInCloud)/(24*pow(PI,1.5)*(EPSILON_0/100)*detectorSpecs.EPSILON_R*sigmaIn);
		// Debug
		// Dprime = detectorSpecs.DIFFUSION_E + detectorSpecs.MU_E*abs(chargeInCloud*1.60217646e-19)/(24*pow(PI,1.5)*(EPSILON_0/100)*detectorSpecs.EPSILON_R*sigmaIn);

	return 2*Dprime*dt;
}

//===================================================================
// This function calculates the output signal on triggered anodes and
// cathodes as a function of time. The result is stored in anVsTime
// and caVsTime.
void Simulation::getDetectorOutput(Detector *detector,  int &numChargeElems, int &maxLen, vector<vector<int> > &trailSizeE, vector<vector<int> > &trailSizeH)
{
	double signal = 0;
	vector<double> sigVsTime;

	//===================================================================
	// Loop over all triggered anodes
	for (int h = 0; h < (int) anCollectorOneHot.size(); h++)
	{
		if (anCollectorOneHot[h] == 1)
		{
			signal = 0;
			sigVsTime.clear();

			// Loop over all time steps
			for (int i = 0; i < maxLen; i++)
			{
				// Loop over all interactions
				for (int j = 0; j < (int) eventTrails.size(); j++)
				{
					// Loop over all charge elements
					for (int k = 0; k < numChargeElems; k++)
					{
						// Signal component due to electrons
						if (i < trailSizeE[j][k])
						{
							// Cathode plane is at z = 0 mm, anode plane is at z = L mm.
							if (i > 0)
								signal += (eventTrails[j][k].qETrapped[i - 1] - eventTrails[j][k].qEMobile[i - 1])*eventAnPhi[h][j][k].EPhi[i - 1];
							if (i < trailSizeE[j][k] - 1)
								signal += eventTrails[j][k].qEMobile[i]*eventAnPhi[h][j][k].EPhi[i];
							else
								if (i == trailSizeE[j][k] - 1)
									signal += eventTrails[j][k].qETrapped[i]*eventAnPhi[h][j][k].EPhi[i];

						}

						// Signal component due to holes
						if (i < trailSizeH[j][k])
						{
							if (i > 0)
								signal += (eventTrails[j][k].qHTrapped[i - 1] - eventTrails[j][k].qHMobile[i - 1])*eventAnPhi[h][j][k].HPhi[i - 1];
							if (i < trailSizeH[j][k] - 1)
								signal += eventTrails[j][k].qHMobile[i]*eventAnPhi[h][j][k].HPhi[i];
							else
								if (i == trailSizeH[j][k] - 1)
									signal += eventTrails[j][k].qHTrapped[i]*eventAnPhi[h][j][k].HPhi[i];
						}
					} // End of loop over charge elements
				} // End of loop over interactions
				
				sigVsTime.push_back(signal);
			}

			anVsTime.push_back(sigVsTime);
		}
	}

	//===================================================================
	// Loop over all triggered cathodes
	for (int h = 0; h < (int) caCollectorOneHot.size(); h++)
	{
		if (caCollectorOneHot[h] == 1)
		{
			signal = 0;
			sigVsTime.clear();

			// Loop over all time steps
			for (int i = 0; i < maxLen; i++)
			{
				// Loop over all interactions
				for (int j = 0; j < (int) eventTrails.size(); j++)
				{
					// Loop over all charge elements
					for (int k = 0; k < numChargeElems; k++)
					{
						// Signal component due to electrons
						if (i < trailSizeE[j][k])
						{
							// Cathode plane is at z = 0 mm, anode plane is at z = L mm.
							if (i > 0)
								signal += (eventTrails[j][k].qETrapped[i - 1] - eventTrails[j][k].qEMobile[i - 1])*eventCaPhi[h][j][k].EPhi[i - 1];
							if (i < trailSizeE[j][k] - 1)
								signal += eventTrails[j][k].qEMobile[i]*eventCaPhi[h][j][k].EPhi[i];
							else
								if (i == trailSizeE[j][k] - 1)
									signal += eventTrails[j][k].qETrapped[i]*eventCaPhi[h][j][k].EPhi[i];
						}

						// Signal component due to holes
						if (i < trailSizeH[j][k])
						{
							if (i > 0)
								signal += (eventTrails[j][k].qHTrapped[i - 1] - eventTrails[j][k].qHMobile[i - 1])*eventCaPhi[h][j][k].HPhi[i - 1];
							if (i < trailSizeH[j][k] - 1)
								signal += eventTrails[j][k].qHMobile[i]*eventCaPhi[h][j][k].HPhi[i];
							else
								if (i == trailSizeH[j][k] - 1)
									signal += eventTrails[j][k].qHTrapped[i]*eventCaPhi[h][j][k].HPhi[i];
						}
					} // End of loop over charge elements
				} // End of loop over interactions
				
				sigVsTime.push_back(signal);
			}

			caVsTime.push_back(sigVsTime);
		}
	}
}

//===================================================================
// This function simulates the preamplifier output based on the
// detector output
void Simulation::getPreampOutput(struct configParamStruct &configParams, Preamplification *preamplifier)
{
	timeVecPreamp.clear();
	double endTime = timeVec.back();
	double T = preamplifier->getSamplingInterval();
	double stdNoise = configParams.FWHM_PREAMP_NOISE/2.3548;
	int index = 0;
	vector<double> tempVec;

	// Construct time sequence with regular sampling interval
	timeVecPreamp.push_back(0);
	while(timeVecPreamp[index] < endTime)
	{
		timeVecPreamp.push_back(timeVecPreamp[index] + T);
		index++;
	}
	timeVecPreamp.pop_back();

	// Convert sequences to ones with regular sampling interval.
	// Loop over all anodes that collected charge
	anVsTimeReg.clear();
	noisyAnVsTimeReg.clear();
	for (int i = 0; i < (int) anVsTime.size(); i++)
	{
		anVsTimeReg.push_back(interp1(timeVec, anVsTime[i], timeVecPreamp));

		// Generate noisy version of the detector output on regular sampling intervals
		tempVec.clear();
		for (int j = 0; j < (int) anVsTimeReg.back().size(); j++)
			tempVec.push_back(anVsTimeReg.back()[j] + randGen->normal()*stdNoise);
		noisyAnVsTimeReg.push_back(tempVec);
	}

	// Loop over all cathodes that collected charge
	caVsTimeReg.clear();
	noisyCaVsTimeReg.clear();
	for (int i = 0; i < (int) caVsTime.size(); i++)
	{
		caVsTimeReg.push_back(interp1(timeVec, caVsTime[i], timeVecPreamp));

		// Generate noisy version of the detector output on regular sampling intervals
		tempVec.clear();
		for (int j = 0; j < (int) caVsTimeReg.back().size(); j++)
			tempVec.push_back(caVsTimeReg.back()[j] + randGen->normal()*stdNoise);
		noisyCaVsTimeReg.push_back(tempVec);
	}

	// Preamplification
	// Loop over all anodes that collected charge
	anVsTimePreamp.clear();
	noisyAnVsTimePreamp.clear();
	for (int i = 0; i < (int) anVsTimeReg.size(); i++)
	{
		anVsTimePreamp.push_back(preamplifier->preamplify(anVsTimeReg[i]));
		noisyAnVsTimePreamp.push_back(preamplifier->preamplify(noisyAnVsTimeReg[i]));
	}

	// Loop over all cathodes that collected charge
	caVsTimePreamp.clear();
	noisyCaVsTimePreamp.clear();
	for (int i = 0; i < (int) caVsTimeReg.size(); i++)
	{
		caVsTimePreamp.push_back(preamplifier->preamplify(caVsTimeReg[i]));
		noisyCaVsTimePreamp.push_back(preamplifier->preamplify(noisyCaVsTimeReg[i]));
	}
}

//===================================================================
// This function returns the sample standard deviation of the spatial
// distribution of charge elements within a cloud for a given
// simulation step. The return value is calculated from only charge
// elements that have not reached the electrode plane.
double Simulation::getSigma(vector<struct trailLog> &interxnTrails, vector<int> &reachedEnd, int &numChargeElems, int currSimStep, int chargeType)
{
	vector<double> x, y, z;
	double x_avg, y_avg, z_avg;
	int numCounted = 0;

	// Push x, y, z position values into vectors
	for(int chargeElem = 0; chargeElem < numChargeElems; chargeElem++) {
		// Only consider the points that have not reached the boundary edges
		if(reachedEnd[chargeElem] != 1)
		{
			// 1 for holes and -1 for electrons
			if (chargeType == 1) {
				x.push_back(interxnTrails[chargeElem].xPosH[currSimStep]);
				y.push_back(interxnTrails[chargeElem].yPosH[currSimStep]);
				z.push_back(interxnTrails[chargeElem].zPosH[currSimStep]);
			}
			else {
				x.push_back(interxnTrails[chargeElem].xPosE[currSimStep]);
				y.push_back(interxnTrails[chargeElem].yPosE[currSimStep]);
				z.push_back(interxnTrails[chargeElem].zPosE[currSimStep]);
			}
			numCounted++;
		}
	}

	if (numCounted == 1)
		return 0;

	// Calculate average
	x_avg = avgVec(x);
	y_avg = avgVec(y);
	z_avg = avgVec(z);

	// Calculate standard deviation
	double sum = 0;
	for(int i = 0; i < numCounted; i++)
		sum += (x[i]-x_avg)*(x[i]-x_avg)+(y[i]-y_avg)*(y[i]-y_avg)+(z[i]-z_avg)*(z[i]-z_avg);

	// Sample standard deviation
	return sqrt(sum/(numCounted - 1));
}

//===================================================================
// This function selects the simulation time increment based on
// electron and hole positions.
void Simulation::getTimeIncr()
{
	double temp = 0, dtDiff = 0;

	// When all electrons of all interactions have reached the boundary, change to larger dt for holes.
	// If at least one interaction has electrons near the anodes, change to finer dt.
	// After all electrons of interaction have reached the boundary, change back to original dt.
	if(numReachedEndE == maxReachedEndE)
		temp = 5e-8;
	else if(eIsCloseToAnode())
		temp = 1e-10;
	else
		temp = 1e-9;

	// dt target changed
	if (temp != targetDt)
	{
		// Set ddt to a 10th of the difference so we'll reach the target dt
		// in 10 steps.
		dtDiff = temp - dt;
		ddt = dtDiff/10;
		dt += ddt;
		targetDt = temp;
	}
	// Target dt did not change
	else
	{
		// If we are not yet at the target dt, then update dt towards target,
		// otherwise do nothing.
		if (abs(dt - targetDt) > 1e-20)
			dt += ddt;
	}
}

//===================================================================
// Given the log of charge elements positions over time, this
// function returns the list of electrodes that triggered, the time
// of trigger as well as the energy detected on the anode and
// cathodes.
//
void Simulation::getTrailPhi(struct configParamStruct &configParams, Detector *detector, int &maxLen, vector<vector<int> > &trailSizeE, vector<vector<int> > &trailSizeH, bool caution)
{
	trailSizeE.clear();
	trailSizeH.clear();

	int numChargeElems = actNumChargeElem;
	// 0 = analytical phi, 1 = numeric phi
	int source = 0;

	if (!strcmp(configParams.phiSource, "n"))
		source = 1;

	try
	{
		// Find length of the longest trail
		maxLen = detector->getEventTrailsSize(eventTrails, numChargeElems, trailSizeE, trailSizeH, caution);
	}
	catch(int)
	{
		throw 0;
	}
	
	// Check if the trail log and time vectors are the same length
	if (maxLen != (int) timeVec.size())
	{
		// Longest charge trail vector has different length from time vector.
		throw 1;
	}

	//===================================================================
	// Loop over anodes that collected charge
	//===================================================================
	for(int h = 0; h < (int) anCollectorOneHot.size(); h++)
	{
		// Work only on anodes that collected charge
		if (anCollectorOneHot[h] == 1)
		{
			// Loop over all simulation time steps
			for (int i = 0; i < maxLen; i++)
			{
				// Loop over all interactions
				for (int j = 0; j < (int) eventTrails.size(); j++)
				{
					// Loop over all charge elements
					for (int k = 0; k < numChargeElems; k++)
					{
						if (i < trailSizeE[j][k])
						{
							if (source == 0)
								eventAnPhi[h][j][k].EPhi.push_back(detector->getWeightingFn2D(detectorSpecs.ANODE_WIDTH*0.0001, eventTrails[j][k].xPosE[i]-detector->getXOffset()[h], detectorSpecs.L-eventTrails[j][k].zPosE[i]));
							else
                                eventAnPhi[h][j][k].EPhi.push_back(detector->getWeightingFn3D(configParams, eventTrails[j][k].xPosE[i]-detector->getXOffset()[h], eventTrails[j][k].yPosE[i], detectorSpecs.L-eventTrails[j][k].zPosE[i], 0));
						}

						if (i < trailSizeH[j][k])
						{
							if (source == 0)
                                eventAnPhi[h][j][k].HPhi.push_back(detector->getWeightingFn2D(detectorSpecs.ANODE_WIDTH*0.0001, eventTrails[j][k].xPosH[i]-detector->getXOffset()[h], detectorSpecs.L-eventTrails[j][k].zPosH[i]));
							else
                                eventAnPhi[h][j][k].HPhi.push_back(detector->getWeightingFn3D(configParams, eventTrails[j][k].xPosH[i]-detector->getXOffset()[h], eventTrails[j][k].yPosH[i], detectorSpecs.L-eventTrails[j][k].zPosH[i], 0));
						}
					} // End of loop over charge elements
				} // End of loop over interactions
			} // End of loop over simulation steps
		}
	} // End of loop over anodes that collected charge

	//===================================================================
	// Loop over cathodes that collected charge
	//===================================================================
	for(int h = 0; h < (int) caCollectorOneHot.size(); h++)
	{
		// Work only on cathodes that collected charge
		if (caCollectorOneHot[h] == 1)
		{
			// Loop over all simulation time steps
			for (int i = 0; i < maxLen; i++)
			{
				// Loop over all interactions
				for (int j = 0; j < (int) eventTrails.size(); j++)
				{
					// Loop over all charge elements
					for (int k = 0; k < numChargeElems; k++)
					{
						if (i < trailSizeE[j][k])
						{
							if (source == 0)
								eventCaPhi[h][j][k].EPhi.push_back(detector->getWeightingFn2D(detectorSpecs.CATHODE_WIDTH*0.0001, eventTrails[j][k].yPosE[i]-detector->getYOffset()[h], eventTrails[j][k].zPosE[i]));
							else
								eventCaPhi[h][j][k].EPhi.push_back(detector->getWeightingFn3D(configParams, eventTrails[j][k].yPosE[i]-detector->getYOffset()[h], eventTrails[j][k].xPosE[i], eventTrails[j][k].zPosE[i], 1));
						}

						if (i < trailSizeH[j][k])
						{
							if (source == 0)
								eventCaPhi[h][j][k].HPhi.push_back(detector->getWeightingFn2D(detectorSpecs.CATHODE_WIDTH*0.0001, eventTrails[j][k].yPosH[i]-detector->getYOffset()[h], eventTrails[j][k].zPosH[i]));
							else
								eventCaPhi[h][j][k].HPhi.push_back(detector->getWeightingFn3D(configParams, eventTrails[j][k].yPosH[i]-detector->getYOffset()[h], eventTrails[j][k].xPosH[i], eventTrails[j][k].zPosH[i], 1));
						}
					} // End of loop over charge elements
				} // End of loop over interactions
			} // End of loop over simulation steps
		}
	} // End of loop over cathodes that collected charge
}

//===================================================================
// This function indicates which anodes triggered 
void Simulation::getTriggers(double &anodeThresh, double &cathodeThresh)
{
	bool trigger = false;
	double maxSig = 0;

	// Get short list of IDs of anodes that collected charge
	for (int i = 0; i < (int) anCollectorOneHot.size(); i++)
	{
		if (anCollectorOneHot[i] == 1)
			anCollectorID.push_back(i);
	}

	// Get short list of IDs of cathodes that collected charge
	for (int i = 0; i < (int) caCollectorOneHot.size(); i++)
	{
		if (caCollectorOneHot[i] == 1)
			caCollectorID.push_back(i);
	}

	// Check the number of anodes that collected charge
	// and number of traces we have are the same
	if ((anVsTimePreamp.size() != anCollectorID.size()) || (caVsTimePreamp.size() != caCollectorID.size()))
		throw 4;

	// Loop over all anodes that collected charge
	for (int i = 0; i < (int) anVsTimePreamp.size(); i++)
	{
		trigger = false;
		maxSig = 0;

		// Loop over all time steps
		for (int j = 0; j < (int) anVsTimePreamp[i].size(); j++)
		{
			// Take note of time of trigger
			if (!trigger)
			{
				if (anVsTimePreamp[i][j] >= anodeThresh)
				{
					double tempTime = 0;

					anodeID.push_back(anCollectorID[i]);
					if (j == 0)
						anTriggerTime.push_back(0);
					else
					{
						// Linearly interpolate the time trigger
						tempTime = timeVecPreamp[j - 1] + (timeVecPreamp[j] - timeVecPreamp[j - 1])*(anodeThresh - anVsTimePreamp[i][j - 1])/(anVsTimePreamp[i][j] - anVsTimePreamp[i][j - 1]);
						anTriggerTime.push_back(tempTime);
					}
					trigger = true;
				}
			}

			// Detector signal peak
			if (anVsTimePreamp[i][j] > maxSig)
				maxSig = anVsTimePreamp[i][j];
		}

		if (trigger)
			anEnergy.push_back(maxSig);
	}

	// Loop over noisy anodes
	for (int i = 0; i < (int) noisyAnVsTimePreamp.size(); i++)
	{
		trigger = false;
		maxSig = 0;

		// Loop over all time steps
		for (int j = 0; j < (int) noisyAnVsTimePreamp[i].size(); j++)
		{
			// Take note of time of trigger
			if (!trigger)
			{
				if (noisyAnVsTimePreamp[i][j] >= anodeThresh)
				{
					double tempTime = 0;

					if (j == 0)
						noisyAnTriggerTime.push_back(0);
					else
					{
						// Linearly interpolate the time trigger
						tempTime = timeVecPreamp[j - 1] + (timeVecPreamp[j] - timeVecPreamp[j - 1])*(anodeThresh - noisyAnVsTimePreamp[i][j - 1])/(noisyAnVsTimePreamp[i][j] - noisyAnVsTimePreamp[i][j - 1]);
						noisyAnTriggerTime.push_back(tempTime);
					}
					trigger = true;
				}
			}

			// Detector signal peak
			if (noisyAnVsTimePreamp[i][j] > maxSig)
				maxSig = noisyAnVsTimePreamp[i][j];
		}

		if (trigger)
			noisyAnEnergy.push_back(maxSig);
	}

	// Loop over all cathodes that collected charge
	for (int i = 0; i < (int) caVsTimePreamp.size(); i++)
	{
		trigger = false;
		maxSig = 0;

		// Loop over all time steps
		for (int j = 0; j < (int) caVsTimePreamp[i].size(); j++)
		{
			// Take note of time of trigger
			if (!trigger)
			{
				if (caVsTimePreamp[i][j] <= cathodeThresh)
				{
					double tempTime = 0;

					cathodeID.push_back(caCollectorID[i]);
					if (j == 0)
						caTriggerTime.push_back(0);
					else
					{
						// Linearly interpolate the time trigger
						tempTime = timeVecPreamp[j - 1] + (timeVecPreamp[j] - timeVecPreamp[j - 1])*(caVsTimePreamp[i][j - 1] - cathodeThresh)/(caVsTimePreamp[i][j - 1] - caVsTimePreamp[i][j]);
						caTriggerTime.push_back(tempTime);
					}
					trigger = true;
				}
			}

			// Detector signal peak
			if (caVsTimePreamp[i][j] < maxSig)
				maxSig = caVsTimePreamp[i][j];
		}

		if (trigger)
			caEnergy.push_back(maxSig);
	}

	// Loop over noisy cathodes
	for (int i = 0; i < (int) noisyCaVsTimePreamp.size(); i++)
	{
		trigger = false;
		maxSig = 0;

		// Loop over all time steps
		for (int j = 0; j < (int) noisyCaVsTimePreamp[i].size(); j++)
		{
			// Take note of time of trigger
			if (!trigger)
			{
				if (noisyCaVsTimePreamp[i][j] <= cathodeThresh)
				{
					double tempTime = 0;

					if (j == 0)
						noisyCaTriggerTime.push_back(0);
					else
					{
						// Linearly interpolate the time trigger
						tempTime = timeVecPreamp[j - 1] + (timeVecPreamp[j] - timeVecPreamp[j - 1])*(noisyCaVsTimePreamp[i][j - 1] - cathodeThresh)/(noisyCaVsTimePreamp[i][j - 1] - noisyCaVsTimePreamp[i][j]);
						noisyCaTriggerTime.push_back(tempTime);
					}
					trigger = true;
				}
			}

			// Detector signal peak
			if (noisyCaVsTimePreamp[i][j] < maxSig)
				maxSig = noisyCaVsTimePreamp[i][j];
		}

		if (trigger)
			noisyCaEnergy.push_back(maxSig);
	}
}

//===================================================================
// This function initializes simulation variables for one event
void Simulation::initSim4Event(struct configParamStruct	&configParams, Event &eventIn)
{
	vector<struct phiLog> intrxnPhis;
	vector<vector<struct phiLog> > eventPhis;

	// Resest time increment
	dt = 1e-9;
	targetDt = dt;
	ddt = 0;

	// Initialize event-wide variables
	for (int i = 0; i < detectorSpecs.NUM_ANODES; i++)
	{
		// anCollectorOneHot keeps track of which anode has collected charge.
		anCollectorOneHot.push_back(0);
	}

	for (int i = 0; i < detectorSpecs.NUM_CATHODES; i++)
	{
		// caCollectorOneHot keeps track of which anode has collected charge.
		caCollectorOneHot.push_back(0);
	}

	// Loop through all interactions of the event
	for (int i = 0; i < (int) eventIn.getNumIntrxns(); i++)
		initSim4Interxn(eventIn, i, actNumChargeElem, configParams.MAX_NUM_SIM_STEPS);

	// Structures for calculating signal induction on electrodes
	for (int i = 0; i < (int) actNumChargeElem; i++)
		intrxnPhis.push_back(chargePhis);
	for (int i = 0; i < (int) eventIn.getNumIntrxns(); i++)
		eventPhis.push_back(intrxnPhis);
	for (int i = 0; i < detectorSpecs.NUM_ANODES; i++)
		eventAnPhi.push_back(eventPhis);
	for (int i = 0; i < detectorSpecs.NUM_CATHODES; i++)
		eventCaPhi.push_back(eventPhis);
	for (int i = 0; i < (int) eventPhis.size(); i++)
		eventPhis[i].clear();
	eventPhis.clear();

	// Variables to keep track of whether the boundary is reached for all cloud elements
	numReachedEndE = 0;
	numReachedEndH = 0;
	maxReachedEndE = int (actNumChargeElem*eventIn.getNumIntrxns());
	maxReachedEnd = 2*maxReachedEndE;
}

//===================================================================
// This function initializes simulation variables for one interaction
void Simulation::initSim4Interxn(Event &eventIn, int intrxnIndex, int numChargeElems, int maxSimSteps)
{
	// Initialize the zero vector
	for (int i = 0; i < actNumChargeElem; i++)
		zeroVec.push_back(0);

	double tempPos = 0;

	// Local vector
	vector<struct trailLog>	intrxnTrails;

	// Calculate the number of electron-hole pairs generated (energy is in eVs)
	double Nq = listGet(eventIn.getEnergy(), intrxnIndex)/detectorSpecs.CZT_W_FACTOR;

	// Accounting for the Fano factor
	double numEHPairs = floor(Nq + 0.5 + randGen->normal()*sqrt(detectorSpecs.FANO*Nq));
	
	// Charge per cloud element
	double qInCloud = numEHPairs/numChargeElems*Q;

	// Loop through all charge elements of the cloud
	for (int i = 0; i < numChargeElems; i++)
	{
		intrxnTrails.push_back(chargeTrails);

		// Populate the first time step elements
		// Electrons are regarded as positive and holes negative due to the incorporating
		// the sign inversion of the shaper, so that anode signals are positive-going and
		// cathode signals are negative-going like we see on the oscilloscope.
		tempPos = listGet(eventIn.getX(), intrxnIndex) + randGen->normal()*sigmaInit*SQRT3_INVERSE;
		tempPos = (tempPos < -detectorSpecs.W/2) ? -detectorSpecs.W/2 : ((tempPos > detectorSpecs.W/2) ? detectorSpecs.W/2 : tempPos);
		intrxnTrails[i].xPosE.push_back(tempPos);
		tempPos = listGet(eventIn.getY(), intrxnIndex) + randGen->normal()*sigmaInit*SQRT3_INVERSE;
		tempPos = (tempPos < -detectorSpecs.W/2) ? -detectorSpecs.W/2 : ((tempPos > detectorSpecs.W/2) ? detectorSpecs.W/2 : tempPos);
		intrxnTrails[i].yPosE.push_back(tempPos);
		tempPos = listGet(eventIn.getZ(), intrxnIndex) + randGen->normal()*sigmaInit*SQRT3_INVERSE;
		tempPos = (tempPos < 0) ? 0 : ((tempPos > detectorSpecs.L) ? detectorSpecs.L : tempPos);
		intrxnTrails[i].zPosE.push_back(tempPos);
		intrxnTrails[i].qETrapped.push_back(0);
		intrxnTrails[i].qEMobile.push_back(qInCloud);

		tempPos = listGet(eventIn.getX(), intrxnIndex) + randGen->normal()*sigmaInit*SQRT3_INVERSE;
		tempPos = (tempPos < -detectorSpecs.W/2) ? -detectorSpecs.W/2 : ((tempPos > detectorSpecs.W/2) ? detectorSpecs.W/2 : tempPos);
		intrxnTrails[i].xPosH.push_back(tempPos);
		tempPos = listGet(eventIn.getY(), intrxnIndex) + randGen->normal()*sigmaInit*SQRT3_INVERSE;
		tempPos = (tempPos < -detectorSpecs.W/2) ? -detectorSpecs.W/2 : ((tempPos > detectorSpecs.W/2) ? detectorSpecs.W/2 : tempPos);
		intrxnTrails[i].yPosH.push_back(tempPos);
		tempPos = listGet(eventIn.getZ(), intrxnIndex) + randGen->normal()*sigmaInit*SQRT3_INVERSE;
		tempPos = (tempPos < 0) ? 0 : ((tempPos > detectorSpecs.L) ? detectorSpecs.L : tempPos);
		intrxnTrails[i].zPosH.push_back(tempPos);
		intrxnTrails[i].qHTrapped.push_back(0);
		intrxnTrails[i].qHMobile.push_back(-qInCloud);
	}

	eventTrails.push_back(intrxnTrails);

	sigmaSqE.push_back(sigmaSqInit);
	sigmaSqH.push_back(sigmaSqInit);

	chargeInCloudE.push_back(Q*numEHPairs);
	chargeInCloudH.push_back(-Q*numEHPairs);

	reachedEndE.push_back(zeroVec);
	reachedEndH.push_back(zeroVec);
	reachedEndESum.push_back(0);
	reachedEndHSum.push_back(0);
	eCloseToAnode.push_back(zeroVec);
}

//===================================================================
// This implements MATLAB's interp1 function:
// The function interpolates to find yReg, the values of the
// underlying function yIrreg at the points in the vector xReg. The
// vector xIrreg specifies the points at which the data yIrreg is
// given.
vector<double> Simulation::interp1(vector<double> &xIrreg, vector<double> &yIrreg, vector<double> &xReg)
{
	vector<double> yReg;
	double xLo = 0, xHi = 0, yLo = 0, yHi = 0;
	int  index = 0;
	
	// Do nothing if xReg is empty
	if ((int) xReg.size() != 0)
	{
		yReg.push_back(yIrreg[0]);
		// Loop till the end of time sequence with regular time interval.
		for (int i = 1; i < (int) xReg.size(); i++)
		{
			// Find point bounding current xReg point
			while (xIrreg[index] <= xReg[i])
				index++;
			
			// Perform linear interpolation
			xLo = xIrreg[index - 1];
			xHi = xIrreg[index];

			yLo = yIrreg[index - 1];
			yHi = yIrreg[index];

			yReg.push_back(yLo + (yHi - yLo)*(xReg[i] - xLo)/(xHi - xLo));
		}
	}

	return yReg;
}

//===================================================================
// This function gets the value from a list at a given index
float Simulation::listGet(list<float> myList, int index)
{
	list<float>::iterator it_myList;
	int count = 0;
	for(it_myList = myList.begin(); it_myList != myList.end(); it_myList++) {
		if(count == index)
			return (*it_myList);
		count++;
	}
	cout << "Error: Invalid index into list get" << endl;
	return 0;
}

//===================================================================
// This function sets the detector specification structure
void Simulation::setDetectorSpecs(struct detectorSpecStruct &detectorSpecsIn)
{
	detectorSpecs.L = detectorSpecsIn.L;
	detectorSpecs.W = detectorSpecsIn.W;
	detectorSpecs.NUM_ANODES = detectorSpecsIn.NUM_ANODES;
	detectorSpecs.NUM_CATHODES = detectorSpecsIn.NUM_CATHODES;
	detectorSpecs.ANODE_PITCH = detectorSpecsIn.ANODE_PITCH;
	detectorSpecs.ANODE_WIDTH = detectorSpecsIn.ANODE_WIDTH;
	detectorSpecs.CATHODE_PITCH = detectorSpecsIn.CATHODE_PITCH;
	detectorSpecs.CATHODE_WIDTH = detectorSpecsIn.CATHODE_WIDTH;
	detectorSpecs.BIAS = detectorSpecsIn.BIAS;
	detectorSpecs.CZT_W_FACTOR = detectorSpecsIn.CZT_W_FACTOR;
	detectorSpecs.DIFFUSION_E = detectorSpecsIn.DIFFUSION_E;
	detectorSpecs.DIFFUSION_H = detectorSpecsIn.DIFFUSION_H;
	detectorSpecs.MU_E = detectorSpecsIn.MU_E;
	detectorSpecs.MU_H = detectorSpecsIn.MU_H;
	detectorSpecs.TAU_E = detectorSpecsIn.TAU_E;
	detectorSpecs.TAU_H = detectorSpecsIn.TAU_H;
	detectorSpecs.FANO = detectorSpecsIn.FANO;
	detectorSpecs.EPSILON_R = detectorSpecsIn.EPSILON_R;
	detectorSpecs.TEMP = detectorSpecsIn.TEMP;
	detectorSpecs.LAMBDA_E = detectorSpecsIn.LAMBDA_E;
	detectorSpecs.LAMBDA_H = detectorSpecsIn.LAMBDA_H;
}

//===================================================================
// This function sets the seed for the Gaussian random number
// generator
void Simulation::setSeed(int seed)
{
	randGen->setSeed(seed);
}

//===================================================================
// This function performs variable assignments to terminate the trail
// of a charge element when it reaches the electrode plane.
void Simulation::terminateTrail(struct trailLog &chargeTrailsIn, int &intrxnIndex, int &chargeElem, int &currSimStep, int chargeType)
{
	double prevXPos = 0, prevYPos = 0, prevZPos = 0, currXPos = 0, currYPos = 0, currZPos = 0;

	// Electron
	if (chargeType < 1)
	{
		// Update flags and sums
		eCloseToAnode[intrxnIndex][chargeElem] = 0;
		reachedEndE[intrxnIndex][chargeElem] = 1;
		reachedEndESum[intrxnIndex]++;
		numReachedEndE++;

		// Interpolate for the boundary value
		prevXPos = chargeTrailsIn.xPosE[currSimStep - 1];
		prevYPos = chargeTrailsIn.yPosE[currSimStep - 1];
		prevZPos = chargeTrailsIn.zPosE[currSimStep - 1];
		currXPos = chargeTrailsIn.xPosE[currSimStep];
		currYPos = chargeTrailsIn.yPosE[currSimStep];
		currZPos = chargeTrailsIn.zPosE[currSimStep];

		chargeTrailsIn.xPosE[currSimStep] = prevXPos + (currXPos - prevXPos)*(prevZPos - detectorSpecs.L)/(prevZPos - currZPos);
		if (chargeTrailsIn.xPosE[currSimStep] < -detectorSpecs.W/2)
			chargeTrailsIn.xPosE[currSimStep] = -detectorSpecs.W/2;
		if (chargeTrailsIn.xPosE[currSimStep] > detectorSpecs.W/2)
			chargeTrailsIn.xPosE[currSimStep] = detectorSpecs.W/2;
		chargeTrailsIn.yPosE[currSimStep] = prevYPos + (currYPos - prevYPos)*(prevZPos - detectorSpecs.L)/(prevZPos - currZPos);
		if (chargeTrailsIn.yPosE[currSimStep] < -detectorSpecs.W/2)
			chargeTrailsIn.yPosE[currSimStep] = -detectorSpecs.W/2;
		if (chargeTrailsIn.yPosE[currSimStep] > detectorSpecs.W/2)
			chargeTrailsIn.yPosE[currSimStep] = detectorSpecs.W/2;
		chargeTrailsIn.zPosE[currSimStep] = detectorSpecs.L;
		chargeTrailsIn.qETrapped.back() += chargeTrailsIn.qEMobile.back();
		chargeTrailsIn.qEMobile.back() = 0;
	}
	// Hole
	else
	{
		reachedEndH[intrxnIndex][chargeElem] = 1;
		reachedEndHSum[intrxnIndex]++;
		numReachedEndH++;

		prevXPos = chargeTrailsIn.xPosH[currSimStep - 1];
		prevYPos = chargeTrailsIn.yPosH[currSimStep - 1];
		prevZPos = chargeTrailsIn.zPosH[currSimStep - 1];
		currXPos = chargeTrailsIn.xPosH[currSimStep];
		currYPos = chargeTrailsIn.yPosH[currSimStep];
		currZPos = chargeTrailsIn.zPosH[currSimStep];

		chargeTrailsIn.xPosH[currSimStep] = prevXPos + (currXPos - prevXPos)*prevZPos/(prevZPos - currZPos);
		if (chargeTrailsIn.xPosH[currSimStep] < -detectorSpecs.W/2)
			chargeTrailsIn.xPosH[currSimStep] = -detectorSpecs.W/2;
		if (chargeTrailsIn.xPosH[currSimStep] > detectorSpecs.W/2)
			chargeTrailsIn.xPosH[currSimStep] = detectorSpecs.W/2;
		chargeTrailsIn.yPosH[currSimStep] = prevYPos + (currYPos - prevYPos)*prevZPos/(prevZPos - currZPos);
		if (chargeTrailsIn.yPosH[currSimStep] < -detectorSpecs.W/2)
			chargeTrailsIn.yPosH[currSimStep] = -detectorSpecs.W/2;
		if (chargeTrailsIn.yPosH[currSimStep] > detectorSpecs.W/2)
			chargeTrailsIn.yPosH[currSimStep] = detectorSpecs.W/2;
		chargeTrailsIn.zPosH[currSimStep] = 0;
		chargeTrailsIn.qHTrapped.back() += chargeTrailsIn.qHMobile.back();
		chargeTrailsIn.qHMobile.back() = 0;
	}
}
