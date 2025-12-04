#pragma once

#include <iostream>
#include <vector>
#include <unordered_set>
#include "Eigen/Dense"
#include "Bead.h"


class Cell {
public:
	Eigen::RowVector3d r; // corner of cell RELATIVE TO ORIGIN... the grid origin diffuses
	std::unordered_set<Bead*> contains; // beads associated inside this gridpoint
	double vol;		                // volume of cell
	static double beadvol; // volume of a bead in the cell.

	static int ntypes;  // number of bead types
	std::vector<double> typenums = std::vector<double>(ntypes); // always up-to-date
	std::vector<double> phis = std::vector<double>(ntypes); // only up-to-date after energy calculation

	static int diag_binsize;
	static int diag_nbins;
	std::vector<double> diag_phis = std::vector<double>(diag_nbins);
	static bool diagonal_linear;

    static bool double_count_main_diagonal;
	static double phi_solvent_max;
	static double phi_chromatin;
	static double kappa;
	static bool density_cap_on;
	static bool compressibility_on;
	static bool diag_pseudobeads_on;
	static bool dense_diagonal_on;
	static int n_small_bins;
	static int n_big_bins;
	static int small_binsize;
	static int big_binsize;
	static int diag_cutoff;
	static int diag_start;
    static bool diagonal_binning;
    static std::vector<int> diagonal_bin_lookup;

	void print();
	void reset();
	void moveIn(Bead* bead);
	void moveOut(Bead* bead);
	double getDensityCapEnergy();
	double getEnergy(const Eigen::MatrixXd &chis);
	double getConstantEnergy(const double constant_chi);
	double getDiagEnergy(const std::vector<double> diag_chis);
	double getBoundaryEnergy(const double boundary_chi, const double delta);
	double getSmatrixEnergy(const Eigen::MatrixXd &Smatrix);
	double getEmatrixEnergy(const Eigen::MatrixXd &Ematrix);
	double getDmatrixEnergy(const Eigen::MatrixXd &Dmatrix);


	double bonds_to_beads(int bonds, int index);
	static int binDiagonal(int d);


};
