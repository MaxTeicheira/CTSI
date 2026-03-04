#ifndef CONFIGPARAMS_H
#define CONFIGPARAMS_H

//===================================================================
// Structure storing run parameters for ctsi executable.
struct configParamStruct{
	char	eFieldFile[256];
	char	phiSource[256];
	char	anodeWeightPotentialFile[256];
	char	cathodeWeightPotentialFile[256];
	char	eventInputFile[256];
	char	detectorSpecFile[256];
	char	preampSpecFile[256];
	char	listmodeOutputName[256];
	int		RANDOM_SEED;
	int		NUM_CHARGE_ELEMS;
	int		MAX_NUM_SIM_STEPS;
	double	SML_PIXEL_DIST;
	double	ANODE_TRIG_THRESH;
	double	CATHODE_TRIG_THRESH;
	double	FWHM_PREAMP_NOISE;
	double	E_FIELD_GRID_SPACE_X;				// E field matrix grid spacing in the x direction
	double	E_FIELD_GRID_SPACE_Y;				// E field matrix grid spacing in the y direction
	double	E_FIELD_GRID_SPACE_Z;				// E field matrix grid spacing in the z direction
	double	AN_PHI_GRID_SPACE_X;				// Anode phi matrix grid spacing in the x direction
	double	AN_PHI_GRID_SPACE_Y;				// Anode phi matrix grid spacing in the y direction
	double	AN_PHI_GRID_SPACE_Z;				// Anode phi matrix grid spacing in the z direction
	double	CA_PHI_GRID_SPACE_X;				// Cathode phi matrix grid spacing in the x direction
	double	CA_PHI_GRID_SPACE_Y;				// Cathode phi matrix grid spacing in the y direction
	double	CA_PHI_GRID_SPACE_Z;				// Cathode phi matrix grid spacing in the z direction
												// Note that scaling is applied before offsetting
	double	EVENT_POS_SCALE_X;					// Amount to scale the x position of events
	double	EVENT_POS_SCALE_Y;					// Amount to scale the y position of events
	double	EVENT_POS_SCALE_Z;					// Amount to scale the z position of events
	double	EVENT_POS_OFFSET_X;					// Amount to offset the x position of events
	double	EVENT_POS_OFFSET_Y;					// Amount to offset the y position of events
	double	EVENT_POS_OFFSET_Z;					// Amount to offset the z position of events
	char	OUTPUT_MODE[256];					// Format of the output data
	bool	CAUTION;							// If caution is on, then data consistency check is performed in run
};

# endif // CONFIGPARAMS_H
