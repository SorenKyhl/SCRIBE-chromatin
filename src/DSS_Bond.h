#pragma once

#include "Bond.h"

// Discrete, Shearable, Stretchable bond (see:)
// Koslover, Spakowitz
// "Discretizing elastic chains for coarse-grained polymer models"
// Soft Matter, 2013
// DOI: 10.1039/C3SM50311A
class DSS_Bond : public Bond {
public:
    DSS_Bond(Bead *bead1, Bead *bead2) : Bond{bead1, bead2} {}

    double delta;
    double ata;
    double gamma;
    double eps_bend;
    double eps_parl;
    double eps_perp;

    double energy() {
        // DELTA = 1
        /*
        double delta = 1;
        double ata =  2.7887;
        double gamma = 0.83281;
        double eps_bend = 1.4668;
        double eps_parl = 34.634;
        double eps_perp = 16.438;
        */

        // DELTA = 0.33
        double delta = 16.5;      // nm
        double ata = 0.152;       // nm-1
        double gamma = 0.938;     // dimless
        double eps_bend = 78.309; // kT nm
        double eps_parl = 2.665;  // kT nm-1
        double eps_perp = 1.942;  // kT nm-1

        double U = 0;

        Eigen::RowVector3d R = pbead2->r - pbead1->r;
        Eigen::RowVector3d Rparl = R.dot(pbead1->u) * pbead1->u;
        Eigen::RowVector3d Rperp = R - Rparl;

        // checks:
        // cout << "-------" << endl;
        // cout << "R is :         " << R << endl;
        // cout << "R|_ + R||   =  " << Rparl + Rperp << endl;
        // cout << "R|_ dot R|| =  " << Rparl.dot(Rperp) << endl;
        // cout << "R|| norm is   : " << Rparl.norm() << endl;
        // cout << "R|| norm is   : " << R.dot(pbead1->u) << endl;

        // U += eps_bend*(u.row(i) - u.row(i-1) - ata*Rperp).squaredNorm();  //
        // bend energy U += eps_parl*pow((R.dot(u.row(i-1)) - delta*gamma), 2);
        // // stretch energy U += eps_perp*Rperp.squaredNorm(); // shear energy

        U += eps_bend * (pbead2->u - pbead1->u - ata * Rperp).squaredNorm();
        U += eps_parl * pow((R.dot(pbead1->u) - delta * gamma), 2);
        U += eps_perp * Rperp.squaredNorm();
        U /= 2 * delta;

        return U;
    }
};
