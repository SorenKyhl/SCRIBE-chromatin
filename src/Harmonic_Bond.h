#pragma once

#include "Bead.h"
#include "Bond.h"

class Harmonic_Bond : public Bond {
public:
	Harmonic_Bond(Bead* bead1, Bead* bead2, double kk, double r00) 
		: Bond{bead1, bead2}, k{kk}, r0{r00} {}

	double k;
	double r0;

	double energy()
	{
		Eigen::RowVector3d displacement = pbead2->r - pbead1->r; 
		double r = sqrt(displacement.dot(displacement));
		return k*(r- r0)*(r - r0); 
	}
};
