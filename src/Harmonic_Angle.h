#pragma once

#include "Angle.h"
#include <cmath>

class Harmonic_Angle : public Angle {
public:
	Harmonic_Angle(Bead* bead1, Bead* bead2, Bead* bead3, double kk) 
		: Angle{bead1, bead2, bead3}, k{kk} {}

	double k;

	double energy()
	{
		Eigen::RowVector3d displacement1 = pbead2->r - pbead1->r; 
		Eigen::RowVector3d displacement2 = pbead3->r - pbead2->r; 
		displacement1 /= displacement1.norm();
		displacement2 /= displacement2.norm();
		return k*(1-displacement1.dot(displacement2));
	}
};
