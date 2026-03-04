//=========================================================================================
// ctsi.cpp : Defines the entry point for the console application.

#include "CTSI.h"
#include "Detector.h"
#include "Event.h"
//#include <TGraph.h>
//#include <TCanvas.h>
//#include <TApplication.h>

using namespace std;

int main(int argc, char* argv[])
{
	// The simulator object
	CTSI ctsiApp;

	//=========================================================================================
	// THE HELP SECTION
	//=========================================================================================
	if (argc == 1)
		ctsiApp.printHelp();

	//=========================================================================================
	// IMPLEMENTATION SECTION
	//=========================================================================================
	char runOption[2];

	if (argc > 1)
	{
		// Import configuration parameters
		if (ctsiApp.importConfigParams(0))
			return 1;

		// Import detector specifications
		if (ctsiApp.importDetectorSpec(0))
			return 1;

		// Import preamplifier specifications
		if (ctsiApp.importPreampSpec(0))
			return 1;

		// Import electric field data
		if (ctsiApp.importEField(0))
			return 1;

		// User specifies numeric weighting potential profile
		if (!strcmp(ctsiApp.getConfigParams().phiSource, "n"))
		{
			// Import weighting potentials
			if (ctsiApp.importPhi(0))
				return 1;
		}

		// Import events
		if (ctsiApp.importEvents(0))
			return 1;

		//=====================================================================================
		// LIST MODE DATA PROCESSING MODE
		//=====================================================================================
		runOption[0] = 'r';
		runOption[1] = '\0';
		if (!strcmp(argv[1], runOption))
		{
			// ROOT:
			//TApplication myApp("ctsi", &argc, argv);

			// Simulate charge transport and signal induction
			ctsiApp.simulateCtsi(0, argv[1]);

			// ROOT:
			//myApp.Run();
		}
		else
		{
			//=====================================================================================
			// INTERACTIVE MODE
			// TO DO: sanitize user input
			//=====================================================================================
			runOption[0] = 'i';
			runOption[1] = '\0';
			if (!strcmp(argv[1], runOption))
			{
				// ROOT:
				// TApplication myApp("ctsi", &argc, argv);
				Event tempEvent;
				float x, y, z, E;
				char keepRunning[256];
				char yes[2], no[2];
				yes[0] = 'y';
				yes[1] = '\0';
				no[0] = 'n';
				no[1] = '\0';

				// Single loop accepting user input for x, z, E values (make infinite?)
				while (true)
				{
					cout << "Enter the values for the event" << endl;
					while(true) {
						cout << "x position (in cm, between " << -ctsiApp.getDetector()->getDetectorSpecs().W/2 << " and " << ctsiApp.getDetector()->getDetectorSpecs().W/2 << ") = ";
						cin >> x;
						if (x >= -ctsiApp.getDetector()->getDetectorSpecs().W/2 && x <= ctsiApp.getDetector()->getDetectorSpecs().W/2)
							break;
					}
					while(true) {
						cout << "y position (in cm, between " << -ctsiApp.getDetector()->getDetectorSpecs().W/2 << " and " << ctsiApp.getDetector()->getDetectorSpecs().W/2 << ") = ";
						cin >> y;
						if (y >= -ctsiApp.getDetector()->getDetectorSpecs().W/2 && y <= ctsiApp.getDetector()->getDetectorSpecs().W/2)
							break;
					}
					while(true) {
						cout << "z position (in cm, between 0 and " << ctsiApp.getDetector()->getDetectorSpecs().L << ") = ";
						cin >> z;
						if(z >= 0 && z <= 0.5)
							break;
					}
					while(true) {
						cout << "Energy (in MeV) = ";
						cin >> E;
						if(E > 0)
							break;
					}

					// Construct a single event
					ctsiApp.getEventArray().clear();
					tempEvent.flushEvent();
					tempEvent.addHit(0,
									 E*(float)1e6,
									 (float) ((x*ctsiApp.getConfigParams().EVENT_POS_SCALE_X) + ctsiApp.getConfigParams().EVENT_POS_OFFSET_X),
									 (float) ((y*ctsiApp.getConfigParams().EVENT_POS_SCALE_Y) + ctsiApp.getConfigParams().EVENT_POS_OFFSET_Y),
									 (float) ((z*ctsiApp.getConfigParams().EVENT_POS_SCALE_Z) + ctsiApp.getConfigParams().EVENT_POS_OFFSET_Z),
									 0);
					ctsiApp.getEventArray().push_back(tempEvent);

					// Simulate charge transport and signal induction
					ctsiApp.simulateCtsi(0, argv[1]);

					// Check if user wants to continue
					while(true) {
						cout << "Continue? (y/n): ";
						cin >> keepRunning;
						if (!strncmp(keepRunning, yes, 1) || !strncmp(keepRunning, no, 1))
							break;
					}

					if (!strncmp(keepRunning, no, 1))
						break;
				}

				// ROOT:
				// myApp.Run();
			}
			else
			{
				cerr << "Error: unrecognized run mode, must be r (list-mode) or i (interactive)." << endl;
				return 1;
			}
		}
	}

	return 0;
}
