#pragma once

#include <cstdint>
#include <iostream>
#include <vector>
#include <unordered_set>
#include "Eigen/Dense"
#include "Bead.h"


class Cell {
public:
	Eigen::RowVector3d r; // corner of cell RELATIVE TO ORIGIN... the grid origin diffuses
	// Beads currently inside this cell. A flat vector (not a hash set): cells
	// hold few beads (bounded by the density cap), so a contiguous scan for
	// moveOut is faster than hashing, and clear()/push_back avoid the per-node
	// malloc/free churn a hash set incurs on every re-mesh and cell crossing.
	std::vector<Bead*> contains;
	double vol;		                // volume of cell
	static double beadvol; // volume of a bead in the cell.

	// Generation stamp used by Sim to deduplicate cells flagged during a single
	// MC move without a hash set (see Sim::flagCell). Compared against Sim's
	// monotonic flag_generation; a stale value reads as "not flagged".
	uint64_t flag_stamp = 0;

	// Cached results of getEnergy / getDiagEnergy. A cell's energy depends only
	// on its contents (typenums, bead ids) and volume, so it is reused until
	// one of those changes: moveIn/moveOut/reset/volume updates call
	// invalidateEnergy(). Most MC moves are accepted, so the "old" energy of a
	// flagged cell is usually the cached "new" energy of an earlier move.
	// Assumes chis/diag_chis are fixed for the lifetime of the simulation.
	bool energy_valid = false;
	bool diag_energy_valid = false;
	double energy_cache = 0;
	double diag_energy_cache = 0;
	void invalidateEnergy() { energy_valid = false; diag_energy_valid = false; }

	static int ntypes;  // number of bead types
	std::vector<double> typenums = std::vector<double>(ntypes); // always up-to-date
	// S n, where S is the symmetric matrix with S_ij = chi_ij from the upper
	// triangle and S_ii = 2 chi_ii, so the plaid sum over i <= j of
	// chi_ij n_i n_j equals n.(S n)/2. Kept up to date like typenums (each
	// bead carries S d), which makes getEnergy O(ntypes), not O(ntypes^2).
	std::vector<double> chi_n = std::vector<double>(ntypes);

	static int diag_binsize;
	static int diag_nbins;
	std::vector<double> diag_phis = std::vector<double>(diag_nbins); // only up-to-date after updateDiagPhis
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

	// Per-separation tables, indexed by genomic separation |i - j| in
	// [0, nbeads): the diagonal bin (-1 if outside [diag_start, diag_cutoff]),
	// the pair count it adds, and diag_chis[bin] * count (0 if outside).
	// Built once by setInteractions so the pair loops do no binning arithmetic.
	static std::vector<int> diag_bin_of;
	static std::vector<int> diag_nbonds_of;
	static std::vector<double> diag_weight_of;
	static void setBeadInteractions(std::vector<Bead> &beads,
	                                const Eigen::MatrixXd &chis);
	static void setInteractions(int nbeads, const std::vector<double> &diag_chis,
	                            bool diagonal_on);

	void print();
	void reset();
	void moveIn(Bead* bead);
	void moveOut(Bead* bead);
	double getDensityCapEnergy();
	double getEnergy(const Eigen::MatrixXd &chis);
	double getConstantEnergy(const double constant_chi);
	double getDiagEnergy(const std::vector<double> &diag_chis);
	void updateDiagPhis();
	double getBoundaryEnergy(const double boundary_chi, const double delta);
	double getSmatrixEnergy(const Eigen::MatrixXd &Smatrix);
	double getEmatrixEnergy(const Eigen::MatrixXd &Ematrix);
	double getDmatrixEnergy(const Eigen::MatrixXd &Dmatrix);


	double bonds_to_beads(int bonds, int index);
	static int binDiagonal(int d);


};
