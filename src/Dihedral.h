#pragma once

#include "Bead.h"

class Dihedral {
public:
    Dihedral(Bead *b1, Bead *b2, Bead *b3, Bead *b4)
        : pbead1{b1}, pbead2{b2}, pbead3{b3}, pbead4{b4} {}

    Bead *pbead1;
    Bead *pbead2;
    Bead *pbead3;
    Bead *pbead4;

    void print() {
        std::cout << pbead1->id << " " << pbead2->id << " " << pbead3->id << " "
                  << pbead4->id << std::endl;
    }
    virtual double energy() = 0;
};
