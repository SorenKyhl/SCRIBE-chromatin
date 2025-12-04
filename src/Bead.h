#pragma once

#include <iostream>
#include <vector>
#include "Eigen/Dense"

class Bead {
public:
    Bead(int i, double x = 0, double y = 0, double z = 0, int ntypes = 1)
        : id{i}, r{x, y, z}, d(ntypes) {}

    Bead() : id{0}, r{0, 0, 0}, d(1) {}

    int id; // unique identifier: uniqueness not strictly enforced.
    Eigen::RowVector3d r;  // position
    Eigen::RowVector3d u;  // orientation
    std::vector<double> d; // bead type assignments

    void print() { std::cout << id << " " << r << std::endl; }
};
