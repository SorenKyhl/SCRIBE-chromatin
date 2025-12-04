#pragma once

#include "Bead.h"

class Angle {
public:
	Angle(Bead* b1, Bead* b2, Bead* b3)
		: pbead1{b1}, pbead2{b2}, pbead3{b3} {}

	Bead* pbead1;
	Bead* pbead2;
	Bead* pbead3;

	void print() {std::cout << pbead1->id <<" "<< pbead2->id <<" "<< pbead3->id << std::endl;}
	virtual double energy() = 0; 
};
