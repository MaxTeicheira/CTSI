#ifndef TRAILLOG_H
#define TRAILLOG_H

#include <vector>

using std::vector;

//===================================================================
// Structure storing charge trail data.
// xPosE[i]
// yPosE[i]
// zPosE[i]
// xPosH[i]
// yPosH[i]
// zPosH[i]
// stores the position of charge element at the BEGINNING of the ith
// step.
//
// qEMobile[i]
// qHMobile[i]
// stores the mobile charge at (x/y/z)Pos(E/H)[i] at the BEGINNING of
// the ith step.
// 
// qETrapped[i]
// qHTrapped[i]
// stores the trapped charge at (x/y/z)Pos(E/H)[i] at the END of the
// ith step.
// xPosE[i]
// yPosE[i]
// zPosE[i]
// xPosH[i]
// yPosH[i]
// zPosH[i]
// stores the position of charge element at the BEGINNING of the ith
// step.
//
// qEMobile[i]
// qHMobile[i]
// stores the mobile charge at (x/y/z)Pos(E/H)[i] at the BEGINNING of
// the ith step.
// 
// qETrapped[i]
// qHTrapped[i]
// stores the trapped charge at (x/y/z)Pos(E/H)[i] at the END of the
// ith step.
struct trailLog{
	vector<double> xPosE;
	vector<double> yPosE;
	vector<double> zPosE;
	vector<double> qETrapped;	// Amount of electron trapped in each simulation step.
								// Last element in the vector is simply the mobile charge in the
								// current simulation step.
	vector<double> qEMobile;	// Amount of mobile electron in current simulation step.
	vector<double> xPosH;
	vector<double> yPosH;
	vector<double> zPosH;
	vector<double> qHTrapped;
	vector<double> qHMobile;	// Amount of mobile holes as function of simulation step number
};

# endif // TRAILLOG_H
