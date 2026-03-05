#include "CTSI.h"
#include "Detector.h"
#include "Event.h"
#include "MyRandom.h"
#include "Preamplification.h"
#include "Simulation.h"

#include <cstring>

using namespace std;

#define NUM_CONFIG_PARAMS_TO_IMPORT 32

CTSI::CTSI()
{
	configParams.eFieldFile[0] = '\0';
	configParams.phiSource[0] = '\0';
	configParams.anodeWeightPotentialFile[0] = '\0';
	configParams.cathodeWeightPotentialFile[0] = '\0';
	configParams.eventInputFile[0] = '\0';
	configParams.detectorSpecFile[0] = '\0';
	configParams.preampSpecFile[0] = '\0';
	configParams.listmodeOutputName[0] = '\0';
	configParams.RANDOM_SEED = 0;
	configParams.NUM_CHARGE_ELEMS = 1;
	configParams.MAX_NUM_SIM_STEPS = 1000;
	configParams.SML_PIXEL_DIST = 0;
	configParams.ANODE_TRIG_THRESH = 0;
	configParams.CATHODE_TRIG_THRESH = 0;
	configParams.FWHM_PREAMP_NOISE = 0;
	configParams.FWHM_PREAMP_NOISE = 0;
	configParams.E_FIELD_GRID_SPACE_X = 0;
	configParams.E_FIELD_GRID_SPACE_Y = 0;
	configParams.E_FIELD_GRID_SPACE_Z = 0;
	configParams.AN_PHI_GRID_SPACE_X = 0;
	configParams.AN_PHI_GRID_SPACE_Y = 0;
	configParams.AN_PHI_GRID_SPACE_Z = 0;
	configParams.CA_PHI_GRID_SPACE_X = 0;
	configParams.CA_PHI_GRID_SPACE_Y = 0;
	configParams.CA_PHI_GRID_SPACE_Z = 0;
	configParams.EVENT_POS_SCALE_X = 1;
	configParams.EVENT_POS_SCALE_Y = 1;
	configParams.EVENT_POS_SCALE_Z = 1;
	configParams.EVENT_POS_OFFSET_X = 0;
	configParams.EVENT_POS_OFFSET_Y = 0;
	configParams.EVENT_POS_OFFSET_Z = 0;
	configParams.OUTPUT_MODE[0] = 'm';
	configParams.OUTPUT_MODE[1] = '\0';
	configParams.CAUTION = true;

	detector = new Detector();
	simulator = new Simulation();
	preamplifier = new Preamplification();

	eventArray.clear();
	fileOpened = false;
	dumpFileOpened = false;

	eventID = 0;
}

CTSI::~CTSI()
{
	delete detector;
	delete simulator;
	delete preamplifier;

	if (fileOpened)
		listmodeFile.close();

	if (dumpFileOpened)
		dumpFile.close();
}

//===================================================================
// This function adds indent.
void CTSI::dispIndent(int indent)
{
	for (int i = 0; i < indent; i++)
		cout << "     ";
}

//===================================================================
// This function is a wrapper function for the dispIndent call back.
void CTSI::dispIndentCallBack(void* ctsiObjectPtr, int indent)
{
    // Make a "this" pointer in this static function by explicitly
	// casting to a pointer to CTSI.
    CTSI* dummyThisPtr = (CTSI*) ctsiObjectPtr;

	dummyThisPtr->dispIndent(indent);
}

//===================================================================
// This function returns the structure of configuration parameters.
struct configParamStruct CTSI::getConfigParams()
{
	return configParams;
}

//===================================================================
// This function returns pointer to the detector member object.
Detector* CTSI::getDetector()
{
	return detector;

}
//===================================================================
// This function returns the list of photon interaction events.
list<Event>& CTSI::getEventArray()
{
	return eventArray;
}

//===================================================================
// This function reads the configuration file and returns the
// configuration parameters in a data structure.
int CTSI::importConfigParams(int indent)
{
	char caption[256];
	char inStr[256], refStr[256];
	int lineNumber = 1, numParamsImported = 0;

	dispIndent(indent);
	cout << "Importing configuration parameters: " << endl;
	dispIndent(indent);
	cout << "Opening ctsi.cfg..." << endl;
	ifstream inFile("config/ctsi.cfg");

	try
	{
		// Read the configuration file
		if (!inFile.good())
		{
			dispIndent(indent);
			cerr << "Unable to open configuration file ctsi.cfg." << endl;
			throw 0;
		}

		dispIndent(indent + 1);
		cout << "Opening ctsi.cfg..." << endl;
		while(!inFile.eof())
		{
			caption[0] = '\0';

			try
			{
				switch (lineNumber)
				{
					case 1:
						// Name of file containing electric field grid
						inFile >> caption >> configParams.eFieldFile;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.eFieldFile << endl;
						numParamsImported++;
						break;
					case 2:
						// How the weighting potential is to be calculated
						inFile >> caption >> configParams.phiSource;
						if (!strcmp(caption, "\0"))
							throw 0;
						if (strcmp(configParams.phiSource, "a") && strcmp(configParams.phiSource, "n"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.phiSource << endl;
						numParamsImported++;
						break;
					case 3:
						// Name of file containing anode weighting potential matrix
						inFile >> caption >> configParams.anodeWeightPotentialFile;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.anodeWeightPotentialFile << endl;
						numParamsImported++;
						break;
					case 4:
						// Name of file containing cathode weighting potential matrix
						inFile >> caption >> configParams.cathodeWeightPotentialFile;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.cathodeWeightPotentialFile << endl;
						numParamsImported++;
						break;
					case 5:
						// Name of file containing event list
						inFile >> caption >> configParams.eventInputFile;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.eventInputFile << endl;
						numParamsImported++;
						break;
					case 6:
						// Name of file containing detector specifications
						inFile >> caption >> configParams.detectorSpecFile;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.detectorSpecFile << endl;
						numParamsImported++;
						break;
					case 7:
						// Name of file containing preamplifier specifications
						inFile >> caption >> configParams.preampSpecFile;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.preampSpecFile << endl;
						numParamsImported++;
						break;
					case 8:
						// Number of list-mode output file
						inFile >> caption >> configParams.listmodeOutputName;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.listmodeOutputName << endl;
						numParamsImported++;
						break;
					case 9:
						// Random seed for the simulator
						inFile >> caption >> configParams.RANDOM_SEED;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.RANDOM_SEED << endl;
						numParamsImported++;
						break;
					case 10:
						// Number of charge elements to use in cloud
						inFile >> caption >> configParams.NUM_CHARGE_ELEMS;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.NUM_CHARGE_ELEMS << endl;
						numParamsImported++;
						break;
					case 11:
						// Maximum number of simulation steps to take before skipping
						// the event (giving up).
						inFile >> caption >> configParams.MAX_NUM_SIM_STEPS;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.MAX_NUM_SIM_STEPS << endl;
						numParamsImported++;
						break;
					case 12:
						// Distance from the cathode plane beyond which the small
						// pixel effect is assumed to be no longer effective and a 
						// smaller simulation time step should be used for electrons.
						inFile >> caption >> configParams.SML_PIXEL_DIST;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.SML_PIXEL_DIST << endl;
						numParamsImported++;
						break;
					case 13:
						// Trigger threshold of anodes
						inFile >> caption >> configParams.ANODE_TRIG_THRESH;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.ANODE_TRIG_THRESH << endl;
						numParamsImported++;
						break;
					case 14:
						// Trigger threshold of cathodes
						inFile >> caption >> configParams.CATHODE_TRIG_THRESH;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.CATHODE_TRIG_THRESH << endl;
						numParamsImported++;
						break;
					case 15:
						// FWHM noise at the output of the preamplifier
						inFile >> caption >> configParams.FWHM_PREAMP_NOISE;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.FWHM_PREAMP_NOISE << endl;
						numParamsImported++;
						break;
					case 16:
						// E field matrix grid spacing in the x direction
						inFile >> caption >> configParams.E_FIELD_GRID_SPACE_X;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.E_FIELD_GRID_SPACE_X << endl;
						numParamsImported++;
						break;
					case 17:
						// E field matrix grid spacing in the y direction
						inFile >> caption >> configParams.E_FIELD_GRID_SPACE_Y;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.E_FIELD_GRID_SPACE_Y << endl;
						numParamsImported++;
						break;
					case 18:
						// E field matrix grid spacing in the z direction
						inFile >> caption >> configParams.E_FIELD_GRID_SPACE_Z;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.E_FIELD_GRID_SPACE_Z << endl;
						numParamsImported++;
						break;
					case 19:
						// // Anode phi matrix grid spacing in the x direction
						inFile >> caption >> configParams.AN_PHI_GRID_SPACE_X;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.AN_PHI_GRID_SPACE_X << endl;
						numParamsImported++;
						break;
					case 20:
						// // Anode phi matrix grid spacing in the y direction
						inFile >> caption >> configParams.AN_PHI_GRID_SPACE_Y;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.AN_PHI_GRID_SPACE_Y << endl;
						numParamsImported++;
						break;
					case 21:
						// // Anode phi matrix grid spacing in the z direction
						inFile >> caption >> configParams.AN_PHI_GRID_SPACE_Z;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.AN_PHI_GRID_SPACE_Z << endl;
						numParamsImported++;
						break;
					case 22:
						// // Cathode phi matrix grid spacing in the x direction
						inFile >> caption >> configParams.CA_PHI_GRID_SPACE_X;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.CA_PHI_GRID_SPACE_X << endl;
						numParamsImported++;
						break;
					case 23:
						// // Cathode phi matrix grid spacing in the y direction
						inFile >> caption >> configParams.CA_PHI_GRID_SPACE_Y;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.CA_PHI_GRID_SPACE_Y << endl;
						numParamsImported++;
						break;
					case 24:
						// // Cathode phi matrix grid spacing in the z direction
						inFile >> caption >> configParams.CA_PHI_GRID_SPACE_Z;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.CA_PHI_GRID_SPACE_Z << endl;
						numParamsImported++;
						break;
					case 25:
						// Amount to scale the x position of events
						inFile >> caption >> configParams.EVENT_POS_SCALE_X;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.EVENT_POS_SCALE_X << endl;
						numParamsImported++;
						break;
					case 26:
						// Amount to scale the y position of events
						inFile >> caption >> configParams.EVENT_POS_SCALE_Y;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.EVENT_POS_SCALE_Y << endl;
						numParamsImported++;
						break;
					case 27:
						// Amount to scale the z position of events
						inFile >> caption >> configParams.EVENT_POS_SCALE_Z;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.EVENT_POS_SCALE_Z << endl;
						numParamsImported++;
						break;
					case 28:
						// Amount to offset the x position of events
						inFile >> caption >> configParams.EVENT_POS_OFFSET_X;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.EVENT_POS_OFFSET_X << endl;
						numParamsImported++;
						break;
					case 29:
						// Amount to offset the y position of events
						inFile >> caption >> configParams.EVENT_POS_OFFSET_Y;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.EVENT_POS_OFFSET_Y << endl;
						numParamsImported++;
						break;
					case 30:
						// Amount to offset the z position of events
						inFile >> caption >> configParams.EVENT_POS_OFFSET_Z;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.EVENT_POS_OFFSET_Z << endl;
						numParamsImported++;
						break;
					case 31:
						// Format of the output data
						inFile >> caption >> configParams.OUTPUT_MODE;
						if (!strcmp(caption, "\0"))
							throw 0;
						dispIndent(indent + 1);
						cout << caption << ": " << configParams.OUTPUT_MODE << endl;
						numParamsImported++;
						break;
					case 32:
						// Flag for data consistency check
						inFile >> caption >> inStr;
						if (!strcmp(caption, "\0"))
							throw 0;
						
						strcpy(refStr, "true");
						if (!strcmp(refStr, inStr))
						{
							configParams.CAUTION = true;
							dispIndent(indent + 1);
							cout << caption << ": true" << endl;
							numParamsImported++;
						}
						else
						{
							strcpy(refStr, "false");
							if (!strcmp(refStr, inStr))
							{
								configParams.CAUTION = false;
								dispIndent(indent + 1);
								cout << caption << ": false" << endl;
								numParamsImported++;
							}
							else
							{
								dispIndent(indent + 1);
								cerr << "Error: unrecognized flag value for data consistency check." << endl;
								throw 0;
							}
							break;
						}
					default:
						inFile >> caption;
				}
			}
			catch (int)
			{
				cout << endl;
				dispIndent(indent + 1);
				cerr << numParamsImported << " of " << NUM_CONFIG_PARAMS_TO_IMPORT << " required configuration parameters imported, import aborted." << endl;
				throw 0;
			}
			lineNumber++;
		}
		inFile.close();

		if (numParamsImported != NUM_CONFIG_PARAMS_TO_IMPORT)
		{
			cout << endl;
			dispIndent(indent + 1);
			cerr << numParamsImported << " of " << NUM_CONFIG_PARAMS_TO_IMPORT << " required configuration parameters imported, import aborted." << endl;
			throw 0;
		}
	}
	catch (int)
	{
		cout << endl;
		dispIndent(indent + 1);
		cout << "Contents of ctsi.cfg must have the following format (PATH optional):" << endl;
		dispIndent(indent + 1);
		cout << "1_E_field                     PATH\\E_FIELD_FILENAME" << endl;
		dispIndent(indent + 1);
		cout << "2_Weighting_potential_source  VALUE (n = numeric, a = analytic)" << endl;
		dispIndent(indent + 1);
		cout << "3_Anode_weighting_potential   PATH\\ANODE_PHI_FILENAME" << endl;
		dispIndent(indent + 1);
		cout << "4_Cathode_weighting_potential PATH\\CATHODE_PHI_FILENAME" << endl;
		dispIndent(indent + 1);
		cout << "5_GRAY_output                 PATH\\GRAY_OUTPUT_FILENAME" << endl;
		dispIndent(indent + 1);
		cout << "6_Detector_specification      PATH\\DET_SPEC_FILENAME" << endl;
		dispIndent(indent + 1);
		cout << "7_Preamplifier_specification  PATH\\PREAMP_SPEC_FILENAME" << endl;
		dispIndent(indent + 1);
		cout << "8_List-mode_output_name       PATH\\LISTMODE_OUTPUT_FILENAME" << endl;
		dispIndent(indent + 1);
		cout << "9_Random_seed                 VALUE" << endl;
		dispIndent(indent + 1);
		cout << "10_Number_of_charge_elements  VALUE" << endl;
		dispIndent(indent + 1);
		cout << "11_Max_num_sim_steps          VALUE" << endl;
		dispIndent(indent + 1);
		cout << "12_Small_pixel_distance       VALUE" << endl;
		dispIndent(indent + 1);
		cout << "13_Anode_trigger_Threshold    VALUE" << endl;
		dispIndent(indent + 1);
		cout << "14_Cathode_trigger_Threshold  VALUE" << endl;
		dispIndent(indent + 1);
		cout << "15_FWHM_preamplifier_noise    VALUE" << endl;
		dispIndent(indent + 1);
		cout << "16_E_field_x_grid_space       VALUE" << endl;
		dispIndent(indent + 1);
		cout << "17_E_field_y_grid_space       VALUE" << endl;
		dispIndent(indent + 1);
		cout << "18_E_field_z_grid_space       VALUE" << endl;
		dispIndent(indent + 1);
		cout << "19_Anode_phi_x_grid_space     VALUE" << endl;
		dispIndent(indent + 1);
		cout << "20_Anode_phi_y_grid_space     VALUE" << endl;
		dispIndent(indent + 1);
		cout << "21_Anode_phi_z_grid_space     VALUE" << endl;
		dispIndent(indent + 1);
		cout << "22_Cathode_phi_x_grid_space   VALUE" << endl;
		dispIndent(indent + 1);
		cout << "23_Cathode_phi_y_grid_space   VALUE" << endl;
		dispIndent(indent + 1);
		cout << "24_Cathode_phi_z_grid_space   VALUE" << endl;
		dispIndent(indent + 1);
		cout << "25_Event_x_pos_scale_factor   VALUE" << endl;
		dispIndent(indent + 1);
		cout << "26_Event_y_pos_scale_factor   VALUE" << endl;
		dispIndent(indent + 1);
		cout << "27_Event_z_pos_scale_factor   VALUE" << endl;
		dispIndent(indent + 1);
		cout << "28_Event_x_pos_offset         VALUE" << endl;
		dispIndent(indent + 1);
		cout << "29_Event_y_pos_offset         VALUE" << endl;
		dispIndent(indent + 1);
		cout << "30_Event_z_pos_offset         VALUE" << endl;
		dispIndent(indent + 1);
		cout << "31_Output_mode                VALUE (m = matrix, r = pseudo rena)" << endl;
		dispIndent(indent + 1);
		cout << "32_Caution                    VALUE (true or false)" << endl;
		cout << endl;
		return 1;
	}

	dispIndent(indent + 1);
	cout << "Done." << endl;

	return 0;
}

//===================================================================
// This function imports the detector specifications through the
// detector member.
int CTSI::importDetectorSpec(int indent)
{
	return detector->importDetectorSpecs((void*) this, 0, CTSI::dispIndentCallBack, configParams.detectorSpecFile);
}

//===================================================================
// This function imports the electric field through the detector
// member.
int CTSI::importEField(int indent)
{
	int returnVal = detector->importEField((void*) this, indent, CTSI::dispIndentCallBack, configParams.eFieldFile, configParams);

	// Debug
	// detector.printEFieldToConsole();

	return returnVal;
}

//===================================================================
// This function reads in events from the GRAY txt output and returns
// the events in a STL vector.
int CTSI::importEvents(int indent)
{
	// Declare variables
	int   eventType, eventID, dir, one, detectorID, currentEventID = -1, currentDir = -1, numEvents = 0;
	float time, energy, x, y, z;
	Event tempEvent;
	eventArray.clear();

	ifstream inFile(configParams.eventInputFile);

	dispIndent(indent);
	cout << "Importing events:" << endl;
	dispIndent(indent + 1);
	cout << "Opening " << configParams.eventInputFile << "..." << endl;

	// Check file status
	if (!inFile.good())
	{
		dispIndent(indent + 1);
		cerr << "Unable to open event file " << configParams.eventInputFile << endl;
		return 1;
	}

	// Read in the first field. If it is the end of file this will fail.
	while (inFile >> dec >> eventType)
	{
		// Read file
		inFile >> dec >> eventID;
		inFile >> dec >> dir;
		inFile >> time;
		inFile >> energy;
		inFile >> x;
		inFile >> y;
		inFile >> z;
		inFile >> dec >> one;
		inFile >> dec >> detectorID;
		
		// Check if it's a new event
		if ((currentEventID != eventID) || (currentDir != dir))
		{
			// Update event ID and photon direction
			currentEventID = eventID;
			currentDir = dir;

			// Add last event to array
			if (numEvents != 0)
			{
				eventArray.push_back(tempEvent);
			}

			// Flush the event
			tempEvent.flushEvent();
			tempEvent.setEventID(eventID);
			numEvents++;
		}
		
		// Append to current event
		// On y position:
		// The -5.95 is to center the detector to the origin. Source was at the origin,
		// the detector was assumed to be 39 mm x 39 mm x 5 mm, with its edge 4 cm
		// from the source, so its center is 4 cm + 39mm/2 = 5.95 cm from the origin.
		// On z position:
		// Add L/2 to shift z-values so the bottom is at z = 0
		tempEvent.addHit(time,
						 energy*(float)1e6,
						 (float) ((x*configParams.EVENT_POS_SCALE_X) + configParams.EVENT_POS_OFFSET_X),
						 (float) ((y*configParams.EVENT_POS_SCALE_Y) + configParams.EVENT_POS_OFFSET_Y),
						 (float) ((z*configParams.EVENT_POS_SCALE_Z) + configParams.EVENT_POS_OFFSET_Z),
						 detectorID);
	}

	eventArray.push_back(tempEvent); // Add last event

	inFile.close();

	dispIndent(indent + 1);
	cout << (unsigned int) eventArray.size() << " events." << endl;
	dispIndent(indent + 1);
	cout << "Done." << endl;

	return 0;
}

//===================================================================
// This function imports the anode and cathode weighting potentials
// through the detector member.
int CTSI::importPhi(int indent)
{
	if (detector->importPhiAnode((void*) this, indent, CTSI::dispIndentCallBack, configParams.anodeWeightPotentialFile))
		return 1;

	if (detector->importPhiCathode((void*) this, indent, CTSI::dispIndentCallBack, configParams.cathodeWeightPotentialFile))
		return 1;

	// Debug
	// detector->querryPhiMat();

	return 0;
}

//===================================================================
// This function imports the coefficients of the discrete filter
// that models the preamplifier.
int CTSI::importPreampSpec(int indent)
{
	return preamplifier->importPreamp((void*) this, 0, CTSI::dispIndentCallBack, configParams.preampSpecFile);
}

//===================================================================
// This function outputs simulation result data in the format
// specified by the user.
void CTSI::output(int indent, char* runMode, int &eventCount,
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
				  bool noMoreEvents)
{
	int tempSize = 0;
	char outOption[2], runOption[2];
	runOption[0] = 'i';
	runOption[1] = '\0';

	if (!noMoreEvents)
	{
		// Create output file if not opened already
		if (!fileOpened)
		{
			// TO DO: writing to output in batches may be faster.

			// Open list-mode output file
			listmodeFile.open(configParams.listmodeOutputName);
			fileOpened = true;

			outOption[0] = 'm';
			outOption[1] = '\0';
			if (!strcmp(configParams.OUTPUT_MODE, outOption))
			{
				cout << endl;
				dispIndent(indent + 1);
				cout << "The output file is in matrix format. The columns are:" << endl;
				dispIndent(indent + 1);
				cout << "1. Event ID" << endl;
				dispIndent(indent + 1);
				cout << "2. Original event ID" << endl;
				dispIndent(indent + 1);
				cout << "3. Electrode type (0 = anode, 1 = cathode)" << endl;
				dispIndent(indent + 1);
				cout << "4. Electrode ID" << endl;
				dispIndent(indent + 1);
				cout << "5. Pulse amplitude" << endl;
				dispIndent(indent + 1);
				cout << "6. Trigger time" << endl;
				dispIndent(indent + 1);
				cout << "7. Pulse amplitude with noise" << endl;
				dispIndent(indent + 1);
				cout << "8. Trigger time with noise" << endl;
				cout << endl;
			}
			else
			{
				outOption[0] = 'r';
				outOption[1] = '\0';
				if (!strcmp(configParams.OUTPUT_MODE, outOption))
				{
					listmodeFile << "Event\tChannel\tCharge Induced\tTimestamp\tType" << endl;
					cout << endl;
					dispIndent(indent + 1);
					cout << "The output file is in pseudo RENA-3 ascii output format." << endl;
					cout << endl;
				}
				else
				{
					cout << endl;
					dispIndent(indent + 1);
					cerr << "Error: " << configParams.OUTPUT_MODE << " is unrecognized as output format specification." << endl;
					throw 0;
				}
			}
		}

		//=====================================================================================
		// WRITE TO OUTPUT
		//=====================================================================================
		// Output in matrix format
		outOption[0] = 'm';
		outOption[1] = '\0';
		if (!strcmp(configParams.OUTPUT_MODE, outOption))
		{
			for (int i = 0; i < (int) anEnergy.size(); i++)
				listmodeFile << eventID << "\t" << eventCount << "\t0\t" << anodeID[i] << "\t" << anEnergy[i] << "\t" << anTriggerTime[i] << "\t" << noisyAnEnergy[i] << "\t" << noisyAnTriggerTime[i] << endl;
			for (int i = 0; i < (int) caEnergy.size(); i++)
				listmodeFile << eventID << "\t" << eventCount << "\t1\t" << cathodeID[i] << "\t" << caEnergy[i] << "\t" << caTriggerTime[i] << "\t" << noisyCaEnergy[i] << "\t" << noisyCaTriggerTime[i] << endl;
		}

		// Output as pseudo RENA output
		outOption[0] = 'r';
		outOption[1] = '\0';
		if (!strcmp(configParams.OUTPUT_MODE, outOption))
		{
			// Log anode signal
			for (int i = 0; i < (int) anEnergy.size(); i++)
				listmodeFile << eventCount << "\t" << anodeID[i] << "\t" << anEnergy[i] << "\t" << anTriggerTime[i] << "\tA" << endl;

			// Log cathode signal, but assume one cathode for now
			for (int i = 0; i < (int) caEnergy.size(); i++)
				listmodeFile << eventCount << "\t" << cathodeID[i] << "\t" << -caEnergy[i] << "\t" << caTriggerTime[i] << "\tC" << endl;
		}

		eventID++;

		//===========================================================================
		// Event dump for interctive mode, this is independent from the other output,
		// which is specified by configParams.OUTPUT_MODE.
		//===========================================================================
		if (!strcmp(runMode, runOption))
		{
			// Open pseudo RENA output file
			dumpFile.open("output/interactiveOut.txt");
			dumpFileOpened = true;

			int numFieldsInTrail = 10;

			cout << endl;
			dispIndent(indent + 1);
			cout << "Event data written to output/interactiveOut.txt, use tools/matlab/interactiveOut.m for parsing." << endl;

			//=======================================================================
			// Write event details to file
			//=======================================================================

			// vector<double> &timeVec
			tempSize = (int) timeVec.size();
			dumpFile << tempSize << endl;
			for (int i = 0; i < tempSize; i++)
				dumpFile << timeVec[i] << endl;

			// vector<int> anCollectorID
			tempSize = (int) anCollectorID.size();
			dumpFile << tempSize << endl;
			for (int i = 0; i < tempSize; i++)
				dumpFile << anCollectorID[i] << endl;
			
			// vector<int> caCollectorID
			tempSize = (int) caCollectorID.size();
			dumpFile << tempSize << endl;
			for (int i = 0; i < tempSize; i++)
				dumpFile << caCollectorID[i] << endl;

			// vector<int> &anodeID
			tempSize = (int) anodeID.size();
			dumpFile << tempSize << endl;
			for (int i = 0; i < tempSize; i++)
				dumpFile << anodeID[i] << endl;
			// vector<double> &anEnergy
			for (int i = 0; i < tempSize; i++)
				dumpFile << anEnergy[i] << endl;
			// vector<double> &anTriggerTime
			for (int i = 0; i < tempSize; i++)
				dumpFile << anTriggerTime[i] << endl;
			// vector<double> &noisyAnEnergy
			for (int i = 0; i < tempSize; i++)
				dumpFile << noisyAnEnergy[i] << endl;
			// vector<double> &noisyAnTriggerTime
			for (int i = 0; i < tempSize; i++)
				dumpFile << noisyAnTriggerTime[i] << endl;

			// vector<vector<double> > &anVsTime
			// Loop over all triggered anodes
			tempSize = (int) anVsTime.size();
			dumpFile << tempSize << endl;
			for (int i = 0; i < (int) anVsTime.size(); i++)
				for (int j = 0; j < (int) timeVec.size(); j++)
					dumpFile << anVsTime[i][j] << endl;

			// vector<int> &cathodeID
			tempSize = (int) cathodeID.size();
			dumpFile << tempSize << endl;
			for (int i = 0; i < tempSize; i++)
				dumpFile << cathodeID[i] << endl;
			// vector<double> &caEnergy
			for (int i = 0; i < tempSize; i++)
				dumpFile << caEnergy[i] << endl;
			// vector<double> &caTriggerTime
			for (int i = 0; i < tempSize; i++)
				dumpFile << caTriggerTime[i] << endl;
			// vector<double> &noisyCaEnergy
			for (int i = 0; i < tempSize; i++)
				dumpFile << noisyCaEnergy[i] << endl;
			// vector<double> &noisyCaTriggerTime
			for (int i = 0; i < tempSize; i++)
				dumpFile << noisyCaTriggerTime[i] << endl;

			// vector<vector<double> > &caVsTime
			// Loop over all triggered cathodes
			tempSize = (int) caVsTime.size();
			dumpFile << tempSize << endl;
			for (int i = 0; i < (int) caVsTime.size(); i++)
				for (int j = 0; j < (int) timeVec.size(); j++)
					dumpFile << caVsTime[i][j] << endl;

			// vector<vector<struct trailLog> > &eventTrailsIn
			dumpFile << numFieldsInTrail << endl;
			// Loop over interactions
			tempSize = (int) eventTrailsIn.size();
			dumpFile << tempSize << endl;
			for (int i = 0; i < tempSize; i++)
			{
				// Loop over charge elements
				dumpFile << (int) eventTrailsIn[i].size() << endl;
				for (int j = 0; j < (int) eventTrailsIn[i].size(); j++)
				{
					// Loop over all time steps for all fields of the trailLog struct
					// xPosE
					dumpFile << (unsigned int) eventTrailsIn[i][j].xPosE.size() << endl;
					for (int h = 0; h < (int) eventTrailsIn[i][j].xPosE.size(); h++)
						dumpFile << eventTrailsIn[i][j].xPosE[h] << endl;
					// yPosE
					dumpFile << (unsigned int) eventTrailsIn[i][j].yPosE.size() << endl;
					for (int h = 0; h < (int) eventTrailsIn[i][j].yPosE.size(); h++)
						dumpFile << eventTrailsIn[i][j].yPosE[h] << endl;
					// zPosE
					dumpFile << (unsigned int) eventTrailsIn[i][j].zPosE.size() << endl;
					for (int h = 0; h < (int) eventTrailsIn[i][j].zPosE.size(); h++)
						dumpFile << eventTrailsIn[i][j].zPosE[h] << endl;
					// qETrapped
					dumpFile << (unsigned int) eventTrailsIn[i][j].qETrapped.size() << endl;
					for (int h = 0; h < (int) eventTrailsIn[i][j].qETrapped.size(); h++)
						dumpFile << eventTrailsIn[i][j].qETrapped[h] << endl;
					// qEMobile
					dumpFile << (unsigned int) eventTrailsIn[i][j].qEMobile.size() << endl;
					for (int h = 0; h < (int) eventTrailsIn[i][j].qEMobile.size(); h++)
						dumpFile << eventTrailsIn[i][j].qEMobile[h] << endl;
					// xPosH
					dumpFile << (unsigned int) eventTrailsIn[i][j].xPosH.size() << endl;
					for (int h = 0; h < (int) eventTrailsIn[i][j].xPosH.size(); h++)
						dumpFile << eventTrailsIn[i][j].xPosH[h] << endl;
					// yPosH
					dumpFile << (unsigned int) eventTrailsIn[i][j].yPosH.size() << endl;
					for (int h = 0; h < (int) eventTrailsIn[i][j].yPosH.size(); h++)
						dumpFile << eventTrailsIn[i][j].yPosH[h] << endl;
					// zPosH
					dumpFile << (unsigned int) eventTrailsIn[i][j].zPosH.size() << endl;
					for (int h = 0; h < (int) eventTrailsIn[i][j].zPosH.size(); h++)
						dumpFile << eventTrailsIn[i][j].zPosH[h] << endl;
					// qHTrapped
					dumpFile << (unsigned int) eventTrailsIn[i][j].qHTrapped.size() << endl;
					for (int h = 0; h < (int) eventTrailsIn[i][j].qHTrapped.size(); h++)
						dumpFile << eventTrailsIn[i][j].qHTrapped[h] << endl;
					// qHMobile
					dumpFile << (unsigned int) eventTrailsIn[i][j].qHMobile.size() << endl;
					for (int h = 0; h < (int) eventTrailsIn[i][j].qHMobile.size(); h++)
						dumpFile << eventTrailsIn[i][j].qHMobile[h] << endl;
				}
			}

			// vector<vector<vector<struct phiLog> > > &eventAnPhi
			// Loop over anodes that collected charge
			tempSize = (int) eventAnPhiIn.size();
			dumpFile << tempSize << endl;
			for (int i = 0; i < tempSize; i++)
			{
				// Loop over interactions
				dumpFile << (int) eventAnPhiIn[i].size() << endl;
				for (int j = 0; j < (int) eventAnPhiIn[i].size(); j++)
				{
					// Loop over charge elements
					dumpFile << (int) eventAnPhiIn[i][j].size() << endl;
					for (int h = 0; h < (int) eventAnPhiIn[i][j].size(); h++)
					{
						// Loop over all time steps for all fields of the phiLog struct
						// EPhi
						dumpFile << (unsigned int) eventAnPhiIn[i][j][h].EPhi.size() << endl;
						for (int k = 0; k < (int) eventAnPhiIn[i][j][h].EPhi.size(); k++)
							dumpFile << eventAnPhiIn[i][j][h].EPhi[k] << endl;

						// HPhi
						dumpFile << (unsigned int) eventAnPhiIn[i][j][h].HPhi.size() << endl;
						for (int k = 0; k < (int) eventAnPhiIn[i][j][h].HPhi.size(); k++)
							dumpFile << eventAnPhiIn[i][j][h].HPhi[k] << endl;
					}
				}
			}

			// vector<vector<vector<struct phiLog> > > &eventCaPhiIn
			// Loop over cathodes that collected charge
			tempSize = (int) eventCaPhiIn.size();
			dumpFile << tempSize << endl;
			for (int i = 0; i < tempSize; i++)
			{
				// Loop over interactions
				dumpFile << (int) eventCaPhiIn[i].size() << endl;
				for (int j = 0; j < (int) eventCaPhiIn[i].size(); j++)
				{
					// Loop over charge elements
					dumpFile << (int) eventCaPhiIn[i][j].size() << endl;
					for (int h = 0; h < (int) eventCaPhiIn[i][j].size(); h++)
					{
						// Loop over all time steps for all fields of the phiLog struct
						// EPhi
						dumpFile << (unsigned int) eventCaPhiIn[i][j][h].EPhi.size() << endl;
						for (int k = 0; k < (int) eventCaPhiIn[i][j][h].EPhi.size(); k++)
							dumpFile << eventCaPhiIn[i][j][h].EPhi[k] << endl;

						// HPhi
						dumpFile << (unsigned int) eventCaPhiIn[i][j][h].HPhi.size() << endl;
						for (int k = 0; k < (int) eventCaPhiIn[i][j][h].HPhi.size(); k++)
							dumpFile << eventCaPhiIn[i][j][h].HPhi[k] << endl;
					}
				}
			}

			// vector<vector<int> > &trailSizeE
			// Loop over interactions
			tempSize = (int) trailSizeE.size();
			dumpFile << tempSize <<endl;
			for (int i = 0; i < tempSize; i++)
			{
				// Loop over charge elements
				dumpFile << (int) trailSizeE[i].size() << endl;
				for (int j = 0; j < (int) trailSizeE[i].size(); j++)
					dumpFile << trailSizeE[i][j] << endl;
			}
			
			// vector<vector<int> > &trailSizeH
			tempSize = (int) trailSizeH.size();
			dumpFile << tempSize <<endl;
			for (int i = 0; i < tempSize; i++)
			{
				// Loop over charge elements
				dumpFile << (int) trailSizeH[i].size() << endl;
				for (int j = 0; j < (int) trailSizeH[i].size(); j++)
					dumpFile << trailSizeH[i][j] << endl;
			}

			// vector<double> &timeVecPreamp
			tempSize = (int) timeVecPreamp.size();
			dumpFile << tempSize << endl;
			for (int i = 0; i < tempSize; i++)
				dumpFile << timeVecPreamp[i] << endl;

			// vector<vector<double> > &anVsTimeReg
			// Loop over all triggered anodes
			for (int i = 0; i < (int) anVsTimeReg.size(); i++)
				for (int j = 0; j < tempSize; j++)
					dumpFile << anVsTimeReg[i][j] << endl;

			// vector<vector<double> > &caVsTimeReg
			// Loop over all triggered cathodes
			for (int i = 0; i < (int) caVsTimeReg.size(); i++)
				for (int j = 0; j < tempSize; j++)
					dumpFile << caVsTimeReg[i][j] << endl;

			// vector<vector<double> > &noisyAnVsTimeReg
			for (int i = 0; i < (int) noisyAnVsTimeReg.size(); i++)
				for (int j = 0; j < tempSize; j++)
					dumpFile << noisyAnVsTimeReg[i][j] << endl;

			// vector<vector<double> > &noisyCaVsTimeReg
			for (int i = 0; i < (int) noisyCaVsTimeReg.size(); i++)
				for (int j = 0; j < tempSize; j++)
					dumpFile << noisyCaVsTimeReg[i][j] << endl;

			// vector<vector<double> > &anVsTimePreamp
			for (int i = 0; i < (int) anVsTimePreamp.size(); i++)
				for (int j = 0; j < tempSize; j++)
					dumpFile << anVsTimePreamp[i][j] << endl;
			
			// vector<vector<double> > &caVsTimePreamp
			for (int i = 0; i < (int) caVsTimePreamp.size(); i++)
				for (int j = 0; j < tempSize; j++)
					dumpFile << caVsTimePreamp[i][j] << endl;

			// vector<vector<double> > &noisyAnVsTimePreamp
			for (int i = 0; i < (int) noisyAnVsTimePreamp.size(); i++)
				for (int j = 0; j < tempSize; j++)
					dumpFile << noisyAnVsTimePreamp[i][j] << endl;

			// vector<vector<double> > &noisyCaVsTimePreamp
			for (int i = 0; i < (int) noisyCaVsTimePreamp.size(); i++)
				for (int j = 0; j < tempSize; j++)
					dumpFile << noisyCaVsTimePreamp[i][j] << endl;

			dumpFile.close();
		}
	}
	else
	{
		// Close file
		listmodeFile.close();
		fileOpened = false;
	}

	//	if(drawClouds == 1)
	//		myDraw(time, caVsTime, 1, 4);

	//	// Store induced cathode and anode pulses for event
	//	vector<double>::iterator it_an_ind;
	//	for(it_an_ind = anEnergy.begin(); it_an_ind != anEnergy.end(); it_an_ind++) {
	//		An.push_back(*it_an_ind);
	//		Ca.push_back(caVsTimeFixed);
	//	}
	//	
	//	// Draw electron/hole cloud trajectories
	//	if(drawClouds == 1) {
	//		myDraw_cld(eventTrailsIn[0][0].xPosE, eventTrailsIn[0][0].zPosE, 1, 2);
	//		myDraw_cld(eventTrailsIn[0][0].xPosH, eventTrailsIn[0][0].zPosH, 0, 4);
	//		for(int i = 1; i < configParams.NUM_CHARGE_ELEMS; i++) {
	//			myDraw_cld(eventTrailsIn[0][i].xPosE, eventTrailsIn[0][i].zPosE, 0, 2);
	//			myDraw_cld(eventTrailsIn[0][i].xPosH, eventTrailsIn[0][i].zPosH, 0, 4);
	//		}
	//		for(int i = 1; i < numIntrxn; i++) {
	//			for(int j = 0; j < configParams.NUM_CHARGE_ELEMS; j++) {
	//				myDraw_cld(eventTrailsIn[i][j].xPosE, eventTrailsIn[i][j].zPosE, 0, 2);
	//				myDraw_cld(eventTrailsIn[i][j].xPosH, eventTrailsIn[i][j].zPosH, 0, 4);
	//			}
	//		}
	//		// Save data to file
	//		ofstream myfile;
	//		myfile.open("clouds_out.txt");
	//		myfile << "x position \t y position \t z position \t charge in cloud" << endl;
	//		for(int intrxnIndex = 0; intrxnIndex < numIntrxn; intrxnIndex++) {
	//			myfile << "interaction " << intrxnIndex << endl;
	//			for(int cld_num = 0; cld_num < configParams.NUM_CHARGE_ELEMS; cld_num++) {
	//				myfile << "cloud "<< cld_num << " electrons" << endl;
	//				for(int i = 0; i < int(eventTrailsIn[intrxnIndex][cld_num].xPosE.size()); i++) {
	//					myfile << eventTrailsIn[intrxnIndex][cld_num].xPosE[i] << "\t" << eventTrailsIn[intrxnIndex][cld_num].yPosE[i] << "\t" ;
	//					myfile << eventTrailsIn[intrxnIndex][cld_num].zPosE[i] << "\t" << eventTrailsIn[intrxnIndex][cld_num].qETrapped[i] << endl;
	//				}
	//				myfile << "cloud "<< cld_num << " holes" << endl;
	//				for(int i = 0; i < int(eventTrailsIn[intrxnIndex][cld_num].xPosH.size()); i++) {
	//					myfile << eventTrailsIn[intrxnIndex][cld_num].xPosH[i] << "\t" << eventTrailsIn[intrxnIndex][cld_num].yPosH[i] << "\t" ;
	//					myfile << eventTrailsIn[intrxnIndex][cld_num].zPosH[i] << "\t" << eventTrailsIn[intrxnIndex][cld_num].qHTrapped[i] << endl;
	//				}
	//			}
	//		}
	//		myfile.close();
	//	}

	//// Make anode vs. cathode scatter plot
	//// Scale anode and cathode values
	//if(drawClouds != 1) {
	//	double scalFac1 = -1/maxVec(Ca)*511;
	//	double scalFac2 = 1/maxVec(An)*511;
	//	//double scalFac1 = -511/(511e3/CZT_W_FACTOR*Q);
	//	//double scalFac2 = 511/(511e3/CZT_W_FACTOR*Q);
	//	// cout << "scaling factor: " << 511/(511e3/CZT_W_FACTOR*Q) << endl;
	//	// cout << -1/maxVec(Ca)*511 << " " << 1/maxVec(An)*511 << endl;
	//	for(int i = 0; i < int(An.size()); i++) {
	//		Ca[i] = Ca[i]*scalFac1;
	//		if(An[i] < 0)
	//			An[i] = 0;
	//		else
	//			An[i] = An[i]*scalFac2;
	//	}
	//	myDraw(Ca, An, 1, -1);
	//}
}

//===================================================================
// This function is a wrapper function for the output call back.
void CTSI::outputCallBack(void* ctsiObjectPtr, int indent, char* runMode, int &eventCount,
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
						  bool noMoreEvents)
{
	// Make a "this" pointer in this static function by explicitly
	// casting to a pointer to CTSI.
    CTSI* dummyThisPtr = (CTSI*) ctsiObjectPtr;

	dummyThisPtr->output(indent, runMode, eventCount,
						 anodeID, anEnergy, anTriggerTime, noisyAnEnergy, noisyAnTriggerTime,
						 cathodeID, caEnergy, caTriggerTime, noisyCaEnergy, noisyCaTriggerTime,
						 anCollectorID, caCollectorID,
						 timeVec, anVsTime, caVsTime,
						 eventTrailsIn, eventAnPhiIn, eventCaPhiIn,
						 trailSizeE, trailSizeH,
						 anVsTimeReg, caVsTimeReg,
						 noisyAnVsTimeReg, noisyCaVsTimeReg,
						 timeVecPreamp, anVsTimePreamp, caVsTimePreamp,
						 noisyAnVsTimePreamp, noisyCaVsTimePreamp,
						 noMoreEvents);
}

//===================================================================
// This function prints the manual text to the console.
void CTSI::printHelp()
{
	// The help section
	cout << "|==========================================================================" << endl;
	cout << "| Function call syntax:                                                   |" << endl;
	cout << "|                                                                         |" << endl;
	cout << "| ctsi.exe r                                                              |" << endl;
	cout << "| Calculates the induced anode and cathode pulses for every interaction   |" << endl;
	cout << "| and outputs the result in list-mode format in an ascii output file.     |" << endl;
	cout << "|                                                                         |" << endl;
	cout << "| ctsi.exe i                                                              |" << endl;
	cout << "| In the interactive mode, the program accepts user input for x, y, z     |" << endl;
	cout << "| pos, and energy via the console; event data is then output. Event data  |" << endl;
	cout << "| includes time series data of charge cloud element position, trapped and |" << endl;
	cout << "| mobile charge, corresponding weighting potential values, electrode time |" << endl;
	cout << "| series output and preamplifier time series output. Data is written to   |" << endl;
	cout << "| interactiveOut.txt, which can be parsed and imported into MATLAB by     |" << endl;
	cout << "| interactiveOut.m.                                                       |" << endl;
	cout << "|                                                                         |" << endl;
	cout << "| Provide an empty ctsi.cfg file and run ctsi to invoke help on the       |" << endl;
	cout << "| correct configuration file content and format.                          |" << endl;
	cout << "|                                                                         |" << endl;
	cout << "| By Yi Gu, 2011 January 12                                               |" << endl;
	cout << "|==========================================================================" << endl;
}

//===================================================================
// This function simulates charge transport and signal induction.
int CTSI::simulateCtsi(int indent, char* runMode)
{
	return simulator->engine((void*) this, indent, CTSI::dispIndentCallBack, CTSI::outputCallBack, configParams, detector, preamplifier, eventArray, runMode);
}

////===================================================================
//// This function returns the maximum absolute value of a vector of doubles
//double CTSI::maxVec(vector<double> &myVec)
//{
//	vector<double>::iterator it_vec;
//	it_vec = myVec.begin();
//	double max = abs(*it_vec);
//	it_vec++;
//	for (; it_vec != myVec.end(); it_vec++)
//	{
//		if(abs(*it_vec) > max)
//			max = abs(*it_vec);
//	}
//
//	return max;
//}

////===================================================================
//// This function makes a scatterplot, given two vectors
//int myDraw(vector<double> &x, vector<double> &y, int first, int j)
//{
//	int n = int(x.size());
//	double *x_arr = (double *)malloc(n*sizeof(double));
//	double *y_arr = (double *)malloc(n*sizeof(double));
//
//	// Convert the vectors to arrays
//	for(int i = 0; i < n; i++) {
//		x_arr[i] = x[i];
//		y_arr[i] = y[i];
//	}
//
//	TGraph *gr1 = new TGraph(n, x_arr, y_arr);
//	gr1->SetMarkerSize(.2);
//
//	// For scatter plot
//	if(j==-1) {
//		TCanvas *c1 = new TCanvas();
//		gr1->Draw("A*");
//		return 0;
//	}
//
//	if(first==1) {
//		TCanvas *c1 = new TCanvas();
//		gr1->SetLineColor(j);
//		gr1->SetMarkerColor(j);
//		gr1->Draw("AL*");
//		/*gr1->SetMinimum(0);
//		gr1->SetMaximum(0.5);*/
//		// c1->RangeAxis(-1.95, 0, 1.95, 0.5);
//		// c1->Update();
//	}
//	else {
//		gr1->SetLineColor(j);
//		gr1->SetMarkerColor(j);
//		gr1->Draw("L*");
//	}
//
//	return 0;
//}
//
////===================================================================
//// This function makes a scatterplot, given two vectors
//int myDraw_cld(vector<double> &x, vector<double> &y, int first, int j)
//{
//	int n = int(x.size());
//	double *x_arr = (double *)malloc(n*sizeof(double));
//	double *y_arr = (double *)malloc(n*sizeof(double));
//
//	// Convert the vectors to arrays
//	for(int i = 0; i < n; i++) {
//		x_arr[i] = x[i];
//		y_arr[i] = y[i];
//	}
//
//	TGraph *gr1 = new TGraph(n, x_arr, y_arr);
//	gr1->SetMarkerSize(.2);
//
//	// For scatter plot
//	if(j==-1) {
//		TCanvas *c1 = new TCanvas();
//		gr1->Draw("A*");
//		return 0;
//	}
//
//	if(first==1) {
//		TCanvas *c1 = new TCanvas();
//		gr1->SetLineColor(j);
//		gr1->SetMarkerColor(j);
//		gr1->Draw("AL*");
//		gr1->SetMinimum(0);
//		gr1->SetMaximum(0.5);
//		// c1->RangeAxis(-1.95, 0, 1.95, 0.5);
//		// c1->Update();
//	}
//	else {
//		gr1->SetLineColor(j);
//		gr1->SetMarkerColor(j);
//		gr1->Draw("L*");
//	}
//
//	return 0;
//}
