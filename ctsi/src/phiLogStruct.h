#ifndef PHILOG_H
#define PHILOG_H

#include <vector>

using std::vector;

//===================================================================
// Structure storing the weighting potential along the charge trails.
struct phiLog {
	vector<double> EPhi;
	vector<double> HPhi;
};

#endif  // PHILOG_H
