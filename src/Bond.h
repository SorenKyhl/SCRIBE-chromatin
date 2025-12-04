#pragma once

#include <iostream>
#include "Bead.h"

class Bond {
public:
	Bond(Bead* b1, Bead* b2)
		: pbead1{b1}, pbead2{b2} {}

	Bead* pbead1;
	Bead* pbead2;

	void print() {std::cout << pbead1->id <<" "<< pbead2->id << std::endl;}
	virtual double energy() = 0; 
};
