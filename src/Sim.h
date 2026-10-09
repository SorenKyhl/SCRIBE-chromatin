#pragma once

#include <chrono>
#include <cstdint>
#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_set>
#include <memory>

#include "Eigen/Dense"
#include "nlohmann/json.hpp"

#include "Bead.h"
#include "DSS_Bond.h"
#include "Harmonic_Bond.h"
#include "Harmonic_Angle.h"
#include "Cell.h"
#include "Grid.h"
#include "Analytics.h"
#include "prof_timer.cpp"
#include "random_mars.h"


class Sim {
public:
	Sim();
	Sim(std::string);
	~Sim();
	
	// members
	std::vector<Bead> beads;
	std::vector<std::unique_ptr<Bond>> bonds; // pointers because Bond class is virtual
	std::vector<std::unique_ptr<Angle>> angles; // pointers because Angle class is virtual
	Grid grid;
	std::unique_ptr<RanMars> rng;  // random number generator

	double chi;
	Eigen::MatrixXd chis;
	std::vector<double> diag_chis;
	double boundary_chi;
	double constant_chi;
	int nspecies; // number of different epigenetic marks
	int nbeads;
	double total_volume;
	double grid_size; // nm
	double bond_length; // nm
	std::string bond_type;
	double k_angle;
	// True only when bead orientation (Bead::u) affects the energy, i.e. for
	// DSS bonds. For gaussian/harmonic bonds the energy depends only on bead
	// positions, so the rotational moves skip all orientation bookkeeping.
	bool orientation_active;
	float dense_diagonal_cutoff;
	float dense_diagonal_loading;
	std::string boundary_type;

	// output files
	bool redirect_stdout = false;
	std::streambuf* cout_stream_buffer; // see: redirectStdout
	std::fstream logfile;
	FILE *xyz_out;
	// stats files stay open for the whole run (openStats), flushed per dump
	FILE *energy_out = nullptr;
	FILE *obs_out = nullptr;
	FILE *diag_obs_out = nullptr;
	FILE *constant_obs_out = nullptr;
	FILE *density_out = nullptr;
	FILE *extra_out = nullptr;
	static FILE *openStats(FILE *&f, const std::string &filename);
	void closeStatsFiles();
	std::string data_out_filename;
	std::string log_filename;
	std::string xyz_out_filename;
	std::string energy_out_filename;
	std::string obs_out_filename;
	std::string diag_obs_out_filename;
	std::string constant_obs_out_filename;
	std::string density_out_filename;
	std::string contact_map_filename;
	std::string extra_out_filename;

	// MC variables
	int decay_length;
	// Mean segment length (beads) of the exponential segment-size
	// distributions. Default nbeads/decay_length in *integer* division (the
	// historical behaviour); each can be overridden from config.
	double exp_decay_trans;
	double exp_decay_crank;
	double exp_decay_pivot;
	double step_disp_percentage = 0.30; // step disp is this fraction of bond length
	double step_trans_percentage = 0.30; // step trans is this fraction of bond length
	double step_disp; // nm
	double step_trans; // nm
	double step_crank = M_PI/6; // radians, full width of the uniform angle interval
	double step_pivot = M_PI/6; // radians, full width of the uniform angle interval
	double step_rot = M_PI/12; // radians
	double step_grid; // based off fraction of delta, see initialize

	Analytics analytics;

	int n_disp;
	int n_trans;
	int n_crank;
	int n_pivot;
	int n_rot;
	// config overrides for the move counts per sweep; -1 = default
	int n_disp_cfg = -1;
	int n_trans_cfg = -1;
	int n_crank_cfg = -1;
	int n_pivot_cfg = -1;

	int acc_disp = 0;
	int acc_trans = 0;
	int acc_crank = 0;
	int acc_pivot = 0;
	int acc_rot = 0;

	int nSteps;// = n_trans + n_crank + n_pivot + n_rot

	// Optional per-move-type statistics for tuning the move set
	// (dump_mixing_stats): attempts, acceptances, beads moved and summed
	// squared bead displacement of accepted moves, and time spent. Written to
	// move_stats.txt at the end of the run, and Rg^2 / end-to-end^2 to
	// structure.traj at every stats dump. Pure bookkeeping: the trajectory is
	// unchanged.
	struct MoveStats {
		long long attempts = 0;
		long long accepts = 0;
		long long beads_moved = 0;
		double sq_disp = 0; // sum over accepted moves of sum_i |dr_i|^2, nm^2
		double seconds = 0;
	};
	bool dump_mixing_stats = false;
	MoveStats stats_disp, stats_trans, stats_crank, stats_pivot, stats_rot, stats_grid;
	FILE *structure_out = nullptr;
	std::string structure_out_filename;
	void saveStructure(int sweep);
	void saveMoveStats();
	int nSweeps;
	int dump_frequency; // dump every x sweeps
	int dump_stats_frequency;
	//int nEquilibSweeps; // = 10*dump_frequency;
	int acc = 0;

	bool bonded_on;
	bool nonbonded_on;
	bool angles_on;

	bool displacement_on;
	bool translation_on;
	bool crankshaft_on;
	bool pivot_on;
	bool rotate_on;

	bool gridmove_on;
	bool diagonal_on;
	bool boundary_attract_on;
	bool constant_chi_on;
	bool plaid_on;
	bool cellcount_on = true;

	bool track_contactmap;
	bool contact_bead_skipping;
    bool conservative_contact_pooling;

	bool load_bead_types;
	bool load_configuration;
	std::string load_configuration_filename;
	std::vector<std::string> bead_type_files;

	// analytics
	bool profiling_on; // = false;
	bool print_trans; // = false;

	unsigned long nbeads_moved = 0;

	// Reusable buffers + generation stamp for flagging the cells touched by a
	// single MC move. Replaces a per-move std::unordered_set<Cell*> /
	// std::unordered_map, whose per-node heap allocations dominated the hot
	// path. flagged_cells_buf holds the unique flagged cells (deduplicated via
	// Cell::flag_stamp); bead_swaps_buf records beads that changed grid cell.
	// Both are cleared (capacity retained) at the start of each move, so no
	// allocation happens after warmup.
	struct BeadSwap { int bead; Cell *old_cell; Cell *new_cell; };
	std::vector<Cell *> flagged_cells_buf;
	std::vector<BeadSwap> bead_swaps_buf;
	uint64_t flag_generation = 0;

	// Start a fresh flagging round for one MC move: bump the generation (so all
	// prior stamps read as stale/unflagged) and clear the reused buffers.
	void beginFlagging() {
		++flag_generation;
		flagged_cells_buf.clear();
		bead_swaps_buf.clear();
	}
	// Record a cell as touched by the current move, deduplicating in O(1).
	void flagCell(Cell *c) {
		if (c->flag_stamp != flag_generation) {
			c->flag_stamp = flag_generation;
			flagged_cells_buf.push_back(c);
		}
	}

	bool set_num_threads;
	int num_threads;

	std::vector<std::vector<int>> contact_map;
	int contact_resolution; //= 500;
	bool dump_density;
	bool visit_tracking;
	bool update_contacts_distance;

	bool smatrix_on;
	bool ematrix_on;
	bool dmatrix_on;

	Eigen::MatrixXd smatrix;
	std::string smatrix_filename;
	Eigen::MatrixXd ematrix;
	std::string ematrix_filename;
	Eigen::MatrixXd dmatrix;
	std::string dmatrix_filename;

	// methods
	void run();
	void xyzToContact();
 
	// contact maps
	void initializeContactmap();
	void updateContacts();
	void updateContactsGridConservative();
	void updateContactsGridNonconservative();
	void updateContactsDistance();

	Eigen::RowVector3d unit_vec(Eigen::RowVector3d b);
	void readInput();
	bool outside_boundary(Eigen::RowVector3d r);
	bool allBeadsInBoundary();
    std::vector<int> generate_diagonal_bin_lookup(std::vector<int> diag_bin_boundaries, int nbeads);
	void setInitialConfiguration();
	void initializeObjects();
	void calculateParameters();
	void volParameters();
	void volParameters_new();
	void loadConfiguration();
	void generateRandomCoil(double bondlength);
	int countLines(std::string filepath);
	void loadBeadTypes();
	void constructBonds();
	void constructAngles();
	void print();
	void checkConsistency();

	// energy calculation
	double getAllBondedEnergy();
	double getBondedEnergy(int first, int last);
	double getNonBondedEnergy(const std::vector<Cell*>& flagged_cells);
	double getJustPlaidEnergy(const std::vector<Cell*>& flagged_cells);
	double getJustDiagEnergy(const std::vector<Cell*>& flagged_cells);
	double getJustBoundaryEnergy(const std::vector<Cell*>& flagged_cells);
	double getTotalEnergy(int first, int last, const std::vector<Cell*>& flagged_cells);

	// Monte Carlo moves
	void MC();
	double randomExp(double mu, double decay);
	void MCmove_displace();
	void MCmove_translate();
	void MCmove_crankshaft();
	void MCmove_rotate();
	void MCmove_pivot(int sweep);
	void MCmove_grid();
	void printAcceptanceRates(int sweep);

	// saving data
	void saveXyz() ;
	void saveEnergy(int sweep);
	void saveObservables(int sweep);
	void saveContacts(int sweep);
	void makeDataAndLogFiles();
	void redirectStdout();
	void returnStdout();
	void makeOutputFiles();
	void setupSmatrix();
	void setupEmatrix();
	void setupDmatrix();
};
