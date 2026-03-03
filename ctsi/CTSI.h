#ifndef CTSI_H
#define CTSI_H

#include "Event.h"
#include "configParamsStruct.h"
#include "phiLogStruct.h"
#include <fstream>
#include <iostream>
#include <list>
#include <vector>

using std::list;
using std::ofstream;
using std::vector;

class Detector;
class Preamplification;
class Simulation;

class CTSI
{
	public:
		CTSI();
		~CTSI();
		void						printHelp();
		struct configParamStruct	getConfigParams();
		Detector*					getDetector();
		int							importConfigParams(int indent);
		int							importDetectorSpec(int indent);
		int							importEField(int indent);
		int							importPhi(int indent);
		int							importPreampSpec(int indent);
		int							importEvents(int indent);
		int							simulateCtsi(int indent, char* runMode);
		void						dispIndent(int indent);
		list<Event> &				getEventArray();
		void						output(int indent, char* runMode, int &eventCount,
										   vector<int> &anodeID, vector<double> &anEnergy, vector<double> &anTriggerTime, vector<double> &noisyAnEnergy, vector<double> &noisyAnTriggerTime,
										   vector<int> &cathodeID, vector<double> &caEnergy, vector<double> &caTriggerTime, vector<double> &noisyCaEnergy, vector<double> &noisyCaTriggerTime,
										   vector<int> &anCollectorID, vector<int> &caCollectorID,
										   vector<double> &timeVec, vector<vector<double> > &anVsTime, vector<vector<double> > &caVsTime,
										   vector<vector<struct trailLog> > &eventTrailsIn, vector<vector<vector<struct phiLog> > > &eventAnPhiIn, vector<vector<vector<struct phiLog> > > &eventCaPhiIn,
										   vector<vector<int> > &trailSizeE, vector<vector<int> > &trailSizeH,
										   vector<vector<double> > &anVsTimeReg, vector<vector<double> > &caVsTimeReg,
										   vector<vector<double> > &noisyAnVsTimeReg, vector<vector<double> > &noisyCaVsTimeReg,
										   vector<double> &timeVecPreamp, vector<vector<double> > &anVsTimePreamp, vector<vector<double> > &caVsTimePreamp,
										   vector<vector<double> > &noisyAnVsTimePreamp, vector<vector<double> > &noisyCaVsTimePreamp,
										   bool noMoreEvents);

		static void					dispIndentCallBack(void* ctsiObjectPtr, int indent);
		static void					outputCallBack(void* ctsiObjectPtr, int indent, char* runMode, int &eventCount,
										   vector<int> &anodeID, vector<double> &anEnergy, vector<double> &anTriggerTime, vector<double> &noisyAnEnergy, vector<double> &noisyAnTriggerTime,
										   vector<int> &cathodeID, vector<double> &caEnergy, vector<double> &caTriggerTime, vector<double> &noisyCaEnergy, vector<double> &noisyCaTriggerTime,
										   vector<int> &anCollectorID, vector<int> &caCollectorID,
										   vector<double> &timeVec, vector<vector<double> > &anVsTime, vector<vector<double> > &caVsTime,
										   vector<vector<struct trailLog> > &eventTrailsIn, vector<vector<vector<struct phiLog> > > &eventAnPhiIn, vector<vector<vector<struct phiLog> > > &eventCaPhiIn,
										   vector<vector<int> > &trailSizeE, vector<vector<int> > &trailSizeH,
										   vector<vector<double> > &anVsTimeReg, vector<vector<double> > &caVsTimeReg,
										   vector<vector<double> > &noisyAnVsTimeReg, vector<vector<double> > &noisyCaVsTimeReg,
										   vector<double> &timeVecPreamp, vector<vector<double> > &anVsTimePreamp, vector<vector<double> > &caVsTimePreamp,
										   vector<vector<double> > &noisyAnVsTimePreamp, vector<vector<double> > &noisyCaVsTimePreamp,
										   bool noMoreEvents);

	protected:		
		struct configParamStruct		configParams;	// Structure storing configuration parameters
		Detector*						detector;
		Simulation*						simulator;
		Preamplification*				preamplifier;	
		list<Event>						eventArray;
		ofstream						listmodeFile;
		ofstream						dumpFile;
		bool							fileOpened;
		bool							dumpFileOpened;
		int								eventID;
};

#endif // CTSI_H
