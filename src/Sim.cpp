#include "Sim.h"
#include <cstdio>

// assign variable to json key of the same name
#define READ_JSON(json, var) read_json((json), (var), #var)
template <class T>
void read_json(const nlohmann::json &json, T &var, std::string varname) {
    std::cout << "loading " << varname << " ... ";
    if (!json.contains(varname)) {
        throw std::runtime_error("config.contains(\"" + varname +
                                 "\") has failed");
    }
    var = json[varname];
    std::cout << "loaded: " << std::to_string(var) << std::endl;
}

// logging information is sent to stdout
Sim::Sim() {
    data_out_filename = "data_out";
    redirect_stdout = false;
    makeDataAndLogFiles();
}

// logging information is sent to {data_out_filename}/log.log
Sim::Sim(std::string filename) {
    data_out_filename = filename;
    log_filename = data_out_filename + "/log.log";
    redirect_stdout = true;
    makeDataAndLogFiles();
}

// return stdout to original stream, if redirected
Sim::~Sim() {
    if (redirect_stdout)
        returnStdout();
}

// Run TICG simulation
void Sim::run() {
    readInput();           // load parameters from config.json
    calculateParameters(); // calculates derived physical parameters
    initializeObjects();   // set particle positions and construct bonds
    saveXyz();             // dump initial configuration
    MC();                  // MC simulation
}

// compute the contactmap for a given .xyz file
// specify input .xyz file in config file
void Sim::xyzToContact() {
    readInput();
    beads.resize(nbeads);
    calculateParameters();
    loadConfiguration(); // load .xyz file specified in config file
    grid.initialize(beads);
    initializeContactmap();
    saveContacts(0);
}

// set contact map to zeros
// contact resolution is the number of beads per contact map pixel
void Sim::initializeContactmap() {
    int nbins = nbeads / contact_resolution;
    contact_map.resize(nbins, std::vector<int>(nbins, 0));
}

// calcualte contacts between adjacent beads
// the contactmap keeps a running sum of contacts
// does not write out out contact map to file
void Sim::updateContacts() {
    if (update_contacts_distance) {
        updateContactsDistance();
    } else {
        if (conservative_contact_pooling)
            updateContactsGridConservative();
        else
            updateContactsGridNonconservative();
    }
}

// accumulates contacts for the current configuration in the contact map.
// two beads are considered in contact if they reside in the same grid cell
// self-self contacts are included
// optional behavior:
// - contact_bead_skipping: (e.g.) if contact_resolution is 5, only check every
// 5th bead
//		for the purposes of contact counting
// - visit_tracking: only count a maximum of one contact per hic-matrix pixel
// for a
//		given configuration (only relevant when contact_resolution > 1)
void Sim::updateContactsGridConservative() {
    std::vector<std::vector<bool>> visited(
        contact_map.size(),
        std::vector<bool>(contact_map.size(), false)); // for visit tracking

    int pixel1, pixel2;
    for (Cell *cell : grid.active_cells) {

        // enumerate pairs of beads - (don't double count!!)
        for (auto bead1 = cell->contains.begin(); bead1 != cell->contains.end(); bead1++) {
            for (auto bead2 = bead1; bead2 != cell->contains.end(); bead2++) {
                bool ignore = false;
                pixel1 = (*bead1)->id / contact_resolution;
                pixel2 = (*bead2)->id / contact_resolution;

                // ignore every other bead for the purposes of contact counting
                // (in the case where contact_resolution=2)
                if (contact_bead_skipping) {
                    if ((*bead1)->id % contact_resolution != 0 ||
                        (*bead2)->id % contact_resolution != 0) {
                        ignore = true;
                    }
                }

                if (visited[pixel1][pixel2] == false) {
                    if (!ignore) {
                        contact_map[pixel1][pixel2] += 1; // increment contacts
                        if (pixel1 != pixel2) {
                            // symmetric, but don't double count main diagonal
                            contact_map[pixel2][pixel1] += 1;
                        }
                    }

                    if (visit_tracking)
                        visited[pixel1][pixel2] = true;
                }
            }
        }
    }
}

void Sim::updateContactsGridNonconservative()
{
    std::vector<std::vector<bool>> visited(contact_map.size(), std::vector<bool>(contact_map.size(), false));   
    int pixel1, pixel2;
	for(Cell* cell : grid.active_cells)
	{
		for(Bead* bead1 : cell->contains)
		{
			for(Bead* bead2 : cell->contains)
			{
				bool ignore = false;
				pixel1 = bead1->id/contact_resolution;
				pixel2 = bead2->id/contact_resolution;

				// ignore every other bead for the purposes of contact counting (in the case where contact_resolution=2)
				if (contact_bead_skipping) {
					if (bead1->id % contact_resolution != 0 || bead2->id % contact_resolution != 0) {
						ignore = true;
					}
				}

				if (visited[pixel1][pixel2] == false) {
					if (!ignore) contact_map[pixel1][pixel2] += 1; // increment contacts
					if (visit_tracking) visited[pixel1][pixel2] = true;
				}
			}
        }
    }
}

// alternative method for calculating bead contacts
// two beads are considered in contact if they are closer than a cutoff radius
void Sim::updateContactsDistance() {
    double cutoff = 28.7; // nm
    int pixel1, pixel2;
    for (std::size_t i = 0; i < beads.size(); i++) {
        for (std::size_t j = i; j < beads.size(); j++) {
            Eigen::RowVector3d distance = beads[i].r - beads[j].r;
            if (distance.norm() < cutoff) {
                pixel1 = beads[i].id / contact_resolution;
                pixel2 = beads[j].id / contact_resolution;
                contact_map[pixel1][pixel2] += 1;
                if (pixel1 != pixel2)
                    contact_map[pixel2][pixel1] += 1;
            }
        }
    }
}

// generates random unit vector
Eigen::MatrixXd Sim::unit_vec(Eigen::MatrixXd b) {
    double R1, R2, R3;
    do {
        R1 = (2 * rng->uniform() - 1);
        R2 = (2 * rng->uniform() - 1);
        R3 = R1 * R1 + R2 * R2;
    } while (R3 >= 1);

    b(0) = 2 * sqrtl(1.0 - R3) * R1;
    b(1) = 2 * sqrtl(1.0 - R3) * R2;
    b(2) = 1 - 2.0 * R3;
    return b;
}

// read simulation parameters from config.json file
void Sim::readInput() {
    std::cout << "reading input config.json file ... " << std::endl;

    std::ifstream i("config.json");
    if (!i.good()) {
        throw std::runtime_error(
            "config.json does not exist or cannot be opened");
    }

    nlohmann::json config;
    i >> config;

    assert(config.contains("nbeads"));
    nbeads = config["nbeads"];
    assert(config.contains("plaid_on"));
    plaid_on = config["plaid_on"];
    assert(config.contains("diagonal_on"));
    diagonal_on = config["diagonal_on"];

    if (plaid_on) {
        assert(config.contains("nspecies"));
        nspecies = config["nspecies"];
        assert(config.contains("load_bead_types"));
        load_bead_types = config["load_bead_types"];

        if (load_bead_types) {
            assert(config.contains("bead_type_files"));
            for (auto file : config["bead_type_files"]) {
                bead_type_files.push_back(file);
            }

            if (bead_type_files.size() != nspecies) {
                throw std::logic_error("Number of bead type files: " +
                                       std::to_string(bead_type_files.size()) +
                                       " must equal number of species: " +
                                       std::to_string(nspecies));
            }

            Cell::ntypes = nspecies;
            int chirows = config["chis"].size();
            int chicols = config["chis"][0].size();

            if (chirows != chicols || chirows != nspecies) {
                throw std::logic_error(
                    "Size of chi matrix : (" + std::to_string(chirows) + ", " +
                    std::to_string(chicols) +
                    ") must be an (n x n) matrix of size nspecies: " +
                    std::to_string(nspecies));
            }

            // set up chi matrix
            chis = Eigen::MatrixXd::Zero(nspecies, nspecies);
            assert(config.contains("chis"));
            for (int i = 0; i < nspecies; i++) {
                for (int j = 0; j < nspecies; j++) {
                    chis(i, j) = config["chis"][i][j];
                }
            }

            std::cout << "chis are: " << std::endl;
            std::cout << chis << std::endl;
        } else {
            nspecies = 0;
        }
    } else {
        nspecies = 0;
        load_bead_types = false;
    }

    if (diagonal_on) {
        assert(config.contains("diag_pseudobeads_on"));
        Cell::diag_pseudobeads_on = config["diag_pseudobeads_on"];
        assert(config.contains("diagonal_linear"));
        Cell::diagonal_linear = config["diagonal_linear"];
        assert(config.contains("dense_diagonal_on"));
        Cell::dense_diagonal_on = config["dense_diagonal_on"];
        assert(config.contains("diag_cutoff"));
        Cell::diag_cutoff = config["diag_cutoff"];
        assert(config.contains("diag_start"));
        Cell::diag_start = config["diag_start"];
        assert(config.contains("diagonal_binning"));
        Cell::diagonal_binning = config["diagonal_binning"];

        assert(config.contains("diag_chis"));
        for (auto e : config["diag_chis"]) {
            diag_chis.push_back(e);
        }

        Cell::diag_nbins = diag_chis.size();

        int ndiag_beads =
            Cell::diag_cutoff - Cell::diag_start; // only count beads that will
                                                  // have diagonal interactions
        if (Cell::dense_diagonal_on) {
            if (config.contains("n_small_bins") &&
                config.contains("n_big_bins")) {
                Cell::n_small_bins = config["n_small_bins"];
                Cell::n_big_bins = config["n_big_bins"];
                assert(Cell::n_small_bins + Cell::n_big_bins ==
                       Cell::diag_nbins);
            } else {
                assert(config.contains("dense_diagonal_loading"));
                dense_diagonal_loading = config["dense_diagonal_loading"];

                Cell::n_small_bins =
                    int(dense_diagonal_loading * Cell::diag_nbins);
                Cell::n_big_bins = Cell::diag_nbins - Cell::n_small_bins;
            }

            if (config.contains("small_binsize") &&
                config.contains("big_binsize")) {
                Cell::small_binsize = config["small_binsize"];
                Cell::big_binsize = config["big_binsize"];
            } else {
                assert(config.contains("dense_diagonal_cutoff"));
                dense_diagonal_cutoff = config["dense_diagonal_cutoff"];
                assert(floorf(ndiag_beads * dense_diagonal_cutoff) ==
                       ndiag_beads * dense_diagonal_cutoff);
                int dividing_line = ndiag_beads * dense_diagonal_cutoff;
                Cell::small_binsize = int(dividing_line / Cell::n_small_bins);
                Cell::big_binsize =
                    int((ndiag_beads - dividing_line) / Cell::n_big_bins);
            }
            std::cout << "number of small bins: " << Cell::n_small_bins
                      << ", of size: " << Cell::small_binsize << std::endl;
            std::cout << "number of big bins: " << Cell::n_big_bins
                      << ", of size: " << Cell::big_binsize << std::endl;
            std::cout << "number of diag beads: " << ndiag_beads << std::endl;
            assert(Cell::n_big_bins + Cell::n_small_bins == Cell::diag_nbins);
            assert(Cell::small_binsize * Cell::n_small_bins +
                       Cell::big_binsize * Cell::n_big_bins ==
                   ndiag_beads);
        }
        else if (Cell::diagonal_binning) 
        {
            std::cout << "DIAGONAL BINNING" << std::endl;
            assert(config.contains("diagonal_bin_boundaries"));
            std::vector<int> diagonal_bin_boundaries;
            for(auto e: config["diagonal_bin_boundaries"])
            {
                diagonal_bin_boundaries.push_back(e);
            }

            if (diagonal_bin_boundaries.size() != diag_chis.size())
            {
                throw std::logic_error("number of diag chis must be equal to number of diagonal bin boundaries");
            }

            Cell::diagonal_bin_lookup = generate_diagonal_bin_lookup(diagonal_bin_boundaries, nbeads);
        }
        else {
            std::cout << "number of bins: " << Cell::diag_nbins << std::endl;
            Cell::diag_binsize = ndiag_beads / diag_chis.size();
            std::cout << "binsize " << Cell::diag_binsize << std::endl;
            assert(ndiag_beads % Cell::diag_binsize == 0);
        }
    }

    assert(config.contains("gridmove_on"));
    gridmove_on = config["gridmove_on"];
    std::cout << "grid move is : " << gridmove_on << std::endl;

    // MC move params
    assert(config.contains("decay_length"));
    decay_length = config["decay_length"];
    assert(config.contains("displacement_on"));
    displacement_on = config["displacement_on"];
    assert(config.contains("translation_on"));
    translation_on = config["translation_on"];
    assert(config.contains("crankshaft_on"));
    crankshaft_on = config["crankshaft_on"];
    assert(config.contains("pivot_on"));
    pivot_on = config["pivot_on"];
    assert(config.contains("rotate_on"));
    rotate_on = config["rotate_on"];

    // energy/chi params
    assert(config.contains("boundary_chi"));
    boundary_chi = config["boundary_chi"];
    assert(config.contains("constant_chi_on"));
    constant_chi_on = config["constant_chi_on"];
    assert(config.contains("constant_chi"));
    constant_chi = config["constant_chi"];
    smatrix_filename = "none";
    ematrix_filename = "none";
    dmatrix_filename = "none";
    assert(config.contains("smatrix_on"));
    smatrix_on = config["smatrix_on"];
    if (config.contains("smatrix_filename")) {
        smatrix_filename = config["smatrix_filename"];
    }
    assert(config.contains("ematrix_on"));
    ematrix_on = config["ematrix_on"];
    if (config.contains("ematrix_filename")) {
        ematrix_filename = config["ematrix_filename"];
    }
    assert(config.contains("dmatrix_on"));
    dmatrix_on = config["dmatrix_on"];
    if (config.contains("dmatrix_filename")) {
        dmatrix_filename = config["dmatrix_filename"];
    }
    assert(config.contains("boundary_attract_on"));
    boundary_attract_on = config["boundary_attract_on"];

    // bead/bond params
    assert(config.contains("beadvol"));
    Cell::beadvol = config["beadvol"];
    assert(config.contains("bond_type"));
    bond_type = config["bond_type"];
    assert(config.contains("bond_length"));
    bond_length = config["bond_length"];
    assert(config.contains("bonded_on"));
    bonded_on = config["bonded_on"];
    assert(config.contains("nonbonded_on"));
    nonbonded_on =
        config["nonbonded_on"]; // TODO does this still work as intended?

    // dump params
    assert(config.contains("dump_frequency"));
    dump_frequency = config["dump_frequency"];
    assert(config.contains("dump_stats_frequency"));
    dump_stats_frequency = config["dump_stats_frequency"];
    assert(config.contains("dump_density"));
    dump_density = config["dump_density"];

    assert(config.contains("nSweeps"));
    nSweeps = config["nSweeps"];
    assert(config.contains("load_configuration"));
    load_configuration = config["load_configuration"];
    assert(config.contains("load_configuration_filename"));
    load_configuration_filename = config["load_configuration_filename"];
    assert(config.contains("profiling_on"));
    profiling_on = config["profiling_on"];
    assert(config.contains("print_trans"));
    print_trans = config["print_trans"];
    assert(config.contains("contact_resolution"));
    contact_resolution = config["contact_resolution"];
    assert(config.contains("grid_size"));
    grid_size = config["grid_size"];
    assert(config.contains("track_contactmap"));
    track_contactmap = config["track_contactmap"];
    assert(config.contains("visit_tracking"));
    visit_tracking = config["visit_tracking"];
    assert(config.contains("update_contacts_distance"));
    update_contacts_distance = config["update_contacts_distance"];
    assert(config.contains("phi_solvent_max"));
    Cell::phi_solvent_max = config["phi_solvent_max"];
    assert(config.contains("phi_chromatin"));
    Cell::phi_chromatin = config["phi_chromatin"];
    assert(config.contains("density_cap_on"));
    Cell::density_cap_on = config["density_cap_on"];
    assert(config.contains("compressibility_on"));
    Cell::compressibility_on = config["compressibility_on"];
    assert(config.contains("kappa"));
    Cell::kappa = config["kappa"];

    assert(config.contains("parallel"));
    Grid::parallel = config["parallel"];
    assert(config.contains("beadvol"));
    Cell::beadvol = config["beadvol"];
    assert(config.contains("cell_volumes"));
    Grid::cell_volumes = config["cell_volumes"];
    assert(config.contains("bond_type"));
    bond_type = config["bond_type"];
    assert(config.contains("bond_length"));
    bond_length = config["bond_length"];
    assert(config.contains("contact_bead_skipping"));
    contact_bead_skipping = config["contact_bead_skipping"];
    assert(config.contains("boundary_type"));
    boundary_type = config["boundary_type"];
    assert(config.contains("angles_on"));
    angles_on = config["angles_on"];
    assert(config.contains("k_angle"));
    k_angle = config["k_angle"];

    // parallel config params
    assert(config.contains("parallel"));
    Grid::parallel = config["parallel"];
    if (Grid::parallel) {
        assert(config.contains("set_num_threads"));
        set_num_threads = config["set_num_threads"];
        if (set_num_threads) {
            assert(config.contains("num_threads"));
            num_threads = config["num_threads"];
            // omp_set_num_threads(num_threads);
        }
    }
    // cellcount_on = config["cellcount_on"];

    assert(config.contains("double_count_main_diagonal"));
    Cell::double_count_main_diagonal = config["double_count_main_diagonal"];

    assert(config.contains("conservative_contact_pooling"));
    conservative_contact_pooling = config["conservative_contact_pooling"];

    assert(config.contains("seed"));
    int seed = config["seed"];
    rng = std::make_unique<RanMars>(seed);

    std::cout << "config_file read successfully" << std::endl;
}

void Sim::makeDataAndLogFiles() {
    std::string command = "mkdir " + data_out_filename;
    const int dir_err = system(command.c_str());

    if (dir_err == -1) {
        std::cout << "error making data_out directory" << std::endl;
    }

    if (redirect_stdout) {
        redirectStdout();
    }
}

void Sim::redirectStdout() {
    cout_stream_buffer = std::cout.rdbuf(); // save original std::cout stream
                                            // buffer for later, to un-redirect
    logfile.open(log_filename, std::ios::out);
    auto cout_buf = std::cout.rdbuf(logfile.rdbuf()); // redirect stdout to log
}

// return std::cout to original stream buffer
// typically done after after done redirecting to log file.
// see: deconstructor
void Sim::returnStdout() { std::cout.rdbuf(cout_stream_buffer); }

void Sim::makeOutputFiles() {
    std::cout << " making output files ... ";

    xyz_out_filename = "./" + data_out_filename + "/output.xyz";
    energy_out_filename = "./" + data_out_filename + "/energy.traj";
    obs_out_filename = "./" + data_out_filename + "/observables.traj";
    diag_obs_out_filename = "./" + data_out_filename + "/diag_observables.traj";
    constant_obs_out_filename =
        "./" + data_out_filename + "/constant_observable.traj";
    density_out_filename = "./" + data_out_filename + "/density.traj";
    extra_out_filename = "./" + data_out_filename + "/extra.traj";

    // fopen(char*, char*) function signature takes c strings...
    // std::string overloads + operator and resturns std::string, need to
    // convert back to c_str()
    xyz_out = fopen(xyz_out_filename.c_str(), "w");
    energy_out = fopen((energy_out_filename).c_str(), "w");
    obs_out = fopen((obs_out_filename).c_str(), "w");
    diag_obs_out = fopen((diag_obs_out_filename).c_str(), "w");
    constant_obs_out = fopen((constant_obs_out_filename).c_str(), "w");
    density_out = fopen((density_out_filename).c_str(), "w");
    extra_out = fopen((extra_out_filename).c_str(), "w");

    std::cout << "created successfully" << std::endl;

    std::string command = "cp config.json " + data_out_filename;
    const int result = system(command.c_str());

    for (const std::string file : bead_type_files) {
        command = "cp " + file + " " + data_out_filename;
        const int result = system(command.c_str());
    }
}

std::vector<int> Sim::generate_diagonal_bin_lookup(std::vector<int> diag_bin_boundaries, int nbeads)
{
    //
    if (diag_bin_boundaries[diag_bin_boundaries.size()-1] != nbeads)
    {
        throw std::logic_error("diagonal binning error; the end of the last diagonal bin must be equal to the numbeer of beads");
    }

    std::vector<int> diag_bin_lookup = std::vector<int>(nbeads, 0);

    int curr = 0;
    int bin_id = 0;
    for (int i = 0; i < nbeads; i++)
    {
        if (i >= diag_bin_boundaries[curr])
        {
            bin_id += 1;
            curr += 1;
        }
        diag_bin_lookup[i] = bin_id;
    }
    return diag_bin_lookup;
}


// check if vector r is inside simulation boundary
bool Sim::outside_boundary(Eigen::RowVector3d r) {
    bool is_out = false;

    if (grid.cubic_boundary) {
        is_out = (r.minCoeff() < 0 || r.maxCoeff() > grid.side_length);
    } else if (grid.spherical_boundary) {
        is_out = (r - grid.sphere_center).norm() > grid.radius;
    }

    return is_out;
}

bool Sim::allBeadsInBoundary() {
    for (const Bead &b : beads) {
        if (outside_boundary(b.r)) {
            return false;
        }
    }
    return true;
}

// set initial bead positions, either by loading from an existing
// .xyz file, or otherwise generating a random coil
void Sim::setInitialConfiguration() {
    beads.resize(nbeads);
    load_configuration ? loadConfiguration() : generateRandomCoil(bond_length);
    assert(allBeadsInBoundary());
}

// initializes all objects prior to simulation
void Sim::initializeObjects() {
    std::cout << "initializing simulation objects ... " << std::endl;
    Timer t_init("Initializing");

    // output files
    makeOutputFiles();

    // physical objects
    setInitialConfiguration();
    if (load_bead_types)
        loadBeadTypes();
    constructBonds();
    if (angles_on)
        constructAngles();

    // energy matrices
    if (dmatrix_on) {
        setupDmatrix();
    }
    if (smatrix_on) {
        setupSmatrix();
    }
    if (ematrix_on) {
        setupEmatrix();
    }

    // grid
    grid.initialize(beads);

    // contactmap
    initializeContactmap();

    std::cout << "Simulation objects initialized" << std::endl;
}

void Sim::volParameters_new() {
    assert(Cell::phi_chromatin > 0 && Cell::phi_chromatin < 1);
    double vol_beads = nbeads * Cell::beadvol;
    double vol = vol_beads / Cell::phi_chromatin;

    if (boundary_type == "cube" || boundary_type == "cubic") {
        grid.cubic_boundary = true;
        grid.spherical_boundary = false;
    }

    if (boundary_type == "sphere" || boundary_type == "spherical") {
        grid.spherical_boundary = true;
        grid.cubic_boundary = false;
    }

    if (grid.cubic_boundary) {
        std::cout << "cubic boundary" << std::endl;
        grid.side_length = std::pow(vol, 1.0 / 3.0);

        grid.L =
            std::ceil(grid.side_length /
                      grid.delta); // number of grid cells per side // ROUNDED,
                                   // won't exactly equal a desired volume frac

        std::cout << "grid.side_length is: " << grid.side_length << std::endl;
        std::cout << "grid.L is: " << grid.L << std::endl;
        total_volume = pow(grid.side_length / 1000.0, 3.0); // [um^3]
        std::cout << "simulation volume is: " << total_volume << " um^3"
                  << std::endl;
        std::cout << "volume fraction is: "
                  << nbeads * Cell::beadvol /
                         (total_volume * 1000 * 1000 * 1000)
                  << std::endl;
    }

    if (grid.spherical_boundary) {
        std::cout << "spherical boundary" << std::endl;
        // float total_volume = 3*vol/(4*M_PI);
        total_volume = vol;
        grid.radius = std::pow(3 * vol / (4 * M_PI),
                               1.0 / 3.0); // [nm] radius of simulation volume
        grid.L = std::ceil(2 * grid.radius / grid.delta);
        grid.sphere_center = {grid.radius, grid.radius, grid.radius};

        std::cout << "grid.radius is " << grid.radius << std::endl;
        std::cout << "grid.L is: " << grid.L << std::endl;
        std::cout << "simulation volume is: "
                  << total_volume / 1000 / 1000 / 1000 << " um^3" << std::endl;
        // std::cout << "sphere_center is : " << grid.sphere_center[0] << " nm "
        // << std::endl;
    }
}

// calculate all derived physical parameters relevant to simulation
void Sim::calculateParameters() {
    grid.delta = grid_size;
    std::cout << "grid size is : " << grid.delta << std::endl;

    // size of Monte-Carlo proposal steps for each type:
    step_grid = grid.delta / 10.0;
    step_disp = step_disp_percentage * bond_length;
    step_trans = step_trans_percentage * bond_length;

    std::cout << "bead volume is : " << Cell::beadvol << std::endl;
    volParameters_new();

    // grid.boundary_radius = std::round(grid.radius); // radius in units of
    // grid cells
    //  sphere center needs to be centered on a multiple of grid delta
    // grid.sphere_center = {grid.boundary_radius*grid.delta,
    // grid.boundary_radius*grid.delta, grid.boundary_radius*grid.delta};

    exp_decay = nbeads / decay_length; // size of exponential falloff for MCmove
                                       // second bead choice
    exp_decay_crank = nbeads / decay_length;
    exp_decay_pivot = nbeads / decay_length;

    if (bond_type == "gaussian" && rotate_on) {
        throw std::runtime_error(
            "bead type is gaussian, set rotate_on = false");
    }

    // number of Monte-Carlo proposal steps for each type
    n_disp = displacement_on ? nbeads : 0;
    n_trans = translation_on ? decay_length : 0;
    n_crank = crankshaft_on ? decay_length : 0;
    n_pivot = pivot_on ? decay_length / 10 : 0;
    n_rot = rotate_on ? nbeads : 0;
    nSteps = n_trans + n_crank + n_pivot + n_rot;
}

void Sim::loadConfiguration() {
    // loads x,y,z positions for every particle from <inputfile>.xyz file

    std::cout << "Loading configuration from " << load_configuration_filename
              << std::endl;
    std::ifstream IFILE;
    IFILE.open(load_configuration_filename);
    if (!IFILE.good()) {
        throw std::runtime_error(load_configuration_filename +
                                 " does not exist or could not be opened");
    }
    std::string line;
    getline(IFILE, line); // nbeads line
    std::cout << line << std::endl;

    int init_nbeads = stoi(line);
    std::cout << "checking if nbeads in config.json matches number of beads in "
                 "the first line of <input>.xyz ... "
              << std::endl;
    assert(init_nbeads == nbeads);

    getline(IFILE, line); // comment line
    std::cout << line << std::endl;

    // first bead
    getline(IFILE, line);
    // std::cout << line << std::endl;
    std::stringstream ss;
    ss << line;
    ss >> beads[0].id;
    ss >> beads[0].r(0);
    ss >> beads[0].r(1);
    ss >> beads[0].r(2);

    for (int i = 1; i < nbeads; i++) {
        getline(IFILE, line);
        // std::cout << line << std::endl;
        std::stringstream ss; // new stream so no overflow from last line
        ss << line;

        ss >> beads[i].id;
        ss >> beads[i].r(0);
        ss >> beads[i].r(1);
        ss >> beads[i].r(2);

        beads[i - 1].u = beads[i].r - beads[i - 1].r;
        beads[i - 1].u = beads[i - 1].u.normalized();
    }
    beads[nbeads - 1].u =
        unit_vec(beads[nbeads - 1].u); // random orientation for last bead
}

void Sim::generateRandomCoil(double bondlength) {
    // generates x,y,z positions for all particles according to a random coil

    double center; // center of simulation box
    if (grid.cubic_boundary) {
        std::cout << " cubic centering" << std::endl;
        center = grid.delta * grid.L / 2;
    } else if (grid.spherical_boundary) {
        std::cout << " sphere centering" << std::endl;
        std::cout << "  assigning center" << std::endl;
        center = grid.sphere_center[0];
        std::cout << "  center is " << center << std::endl;
    } else {
        std::cout << " did not center" << std::endl;
    }

    beads[0].r = {center, center, center}; // start in middlle of the box
    beads[0].u = unit_vec(beads[0].u);
    beads[0].id = 0;

    for (int i = 1; i < nbeads; i++) {
        do {
            beads[i].u = unit_vec(beads[i].u);
            beads[i].r =
                beads[i - 1].r +
                bondlength *
                    beads[i].u; // orientations DO NOT point along contour
            beads[i].id = i;
        } while (outside_boundary(beads[i].r));
    }
}

int Sim::countLines(std::string filepath) {
    int count = 0;
    std::string line;

    std::ifstream file(filepath);
    while (getline(file, line))
        count++;
    return count;
}

void Sim::loadBeadTypes() {
    // set up bead types
    int marktype = 0;
    for (std::string bead_type_file : bead_type_files) {
        std::ifstream IFBEADTYPE;
        IFBEADTYPE.open(bead_type_file);

        int nlines = countLines(bead_type_file);
        if (nlines != nbeads) {
            throw std::runtime_error(
                bead_type_file + " (length : " + std::to_string(nlines) +
                ") is not the right size for a simulation with " +
                std::to_string(nbeads) + " particles.");
        }

        if (IFBEADTYPE.good()) {
            for (int i = 0; i < nbeads; i++) {
                // beads[i].d.reserve(nspecies);
                beads[i].d.resize(nspecies);
                IFBEADTYPE >> beads[i].d[marktype];
            }
            marktype++;
            IFBEADTYPE.close();
        } else {
            throw std::runtime_error(bead_type_file +
                                     " does not exist or could not be opened");
        }
    }
    std::cout << " loaded beads, first bead, first mark:" << beads[0].d[0]
              << std::endl;
}

void Sim::constructBonds() {
    bonds.resize(nbeads - 1); // use default constructor
    for (int i = 0; i < nbeads - 1; i++) {
        if (bond_type == "DSS") {
            bonds[i] = std::make_unique<DSS_Bond>(&beads[i], &beads[i + 1]);
        }
        if (bond_type == "gaussian") {
            // gaussian coil spring constant is 3/(2b^2)
            double k = 3 / (2 * bond_length * bond_length);
            bonds[i] =
                std::make_unique<Harmonic_Bond>(&beads[i], &beads[i + 1], k, 0);
        }
    }
    std::cout << " bonds constructed " << std::endl;
}

void Sim::constructAngles() {
    angles.resize(nbeads - 2);
    for (int i = 0; i < nbeads - 2; i++) {
        if (bond_type == "gaussian" && angles_on) {
            angles[i] = std::make_unique<Harmonic_Angle>(
                &beads[i], &beads[i + 1], &beads[i + 2], k_angle);
        }
    }
    std::cout << "angles constructed" << std::endl;
}

void Sim::print() {
    std::cout << "simulation in : ";
    std::cout << "With beads: " << std::endl;
    for (auto &bead : beads)
        bead.print(); // use reference to avoid copies
    std::cout << "And bonds: " << std::endl;
    for (auto &bond : bonds)
        bond->print();
}

double Sim::getAllBondedEnergy() {
    double U = 0;
    for (auto &bond : bonds) {
        U += bond->energy();
    }

    if (angles_on) {
        for (auto &angle : angles) {
            U += angle->energy();
        }
    }

    return U;
}

double Sim::getBondedEnergy(int first, int last) {
    double U = 0;
    // for(Bond* bo : bonds) U += bo->energy();  // inefficient
    if (first > 0)
        U += bonds[first - 1]->energy(); // move affects bond going into first
    if (last < (nbeads - 1))
        U += bonds[last]->energy(); // ... and leaving the second

    if (angles_on) {
        if (first < (nbeads - 1) && first >= 2) {
            U += angles[first - 2]->energy();
            U += angles[first - 1]->energy();
        }
        if (last < (nbeads - 2) && last >= 1) {
            U += angles[last]->energy();
            U += angles[last - 1]->energy();
        }

        // the above is equivalent to the limiting case below:
        // for(Angle* angle : angles) {U += angle->energy();}
    }
    return U;
}

double
Sim::getNonBondedEnergy(const std::unordered_set<Cell *> &flagged_cells) {
    // gets all the nonbonded energy
    auto start = std::chrono::high_resolution_clock::now();

    double U = grid.densityCapEnergy(flagged_cells);
    if (plaid_on) {
        if (smatrix_on) {
            U += grid.SmatrixEnergy(flagged_cells, smatrix);
        } else if (ematrix_on) {
            U += grid.EmatrixEnergy(flagged_cells, ematrix);
        } else {
            U += grid.energy(flagged_cells, chis);
        }
    }

    if (constant_chi > 0) {
        U += grid.constantEnergy(flagged_cells, constant_chi);
    }

    if (diagonal_on) {
        if (dmatrix_on) {
            U += grid.DmatrixEnergy(flagged_cells, dmatrix);
        } else {
            U += grid.diagEnergy(flagged_cells, diag_chis);
        }
    }

    if (boundary_attract_on) {
        U += grid.boundaryEnergy(flagged_cells, boundary_chi);
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    // std::cout << "NonBonded took " << duration.count() << "microseseconds "<<
    // std::endl;
    return U;
}

double
Sim::getJustPlaidEnergy(const std::unordered_set<Cell *> &flagged_cells) {
    // for when dumping energy;
    double U = grid.energy(flagged_cells, chis);
    return U;
}

double
Sim::getJustBoundaryEnergy(const std::unordered_set<Cell *> &flagged_cells) {
    // for when dumping energy;
    double U = grid.boundaryEnergy(flagged_cells, boundary_chi);
    return U;
}

double Sim::getTotalEnergy(int first, int last,
                           const std::unordered_set<Cell *> &flagged_cells) {
    double U = 0;
    if (bonded_on)
        U += getBondedEnergy(first, last);
    if (nonbonded_on)
        U += getNonBondedEnergy(flagged_cells);
    return U;
}

double Sim::randomExp(double mu, double decay) {
    // generates random number distributed according to two-sided exponential
    // distribution centered about mu, with characteristic decay length
    double cdf_y;
    do {
        cdf_y = rng->uniform();
    } while (
        cdf_y <=
        0); // cdf_y cannot exactly equal 0, otherwise inverse cdf is -infty

    if (cdf_y > 0.5) {
        return mu - decay * log(1 - 2 * std::abs(cdf_y - 0.5)); // inverse cdf
    } else {
        return mu + decay * log(1 - 2 * std::abs(cdf_y - 0.5)); // inverse cdf
    }
}

void Sim::MC() {
    std::cout << "Beginning Simulation" << std::endl;

    analytics.startTimer();
    checkConsistency();

    for (int sweep = 1; sweep < nSweeps + 1; sweep++) {
        double nonbonded;

        Timer t_pivot("pivoting", profiling_on);
        for (int j = 0; j < n_pivot; j++) {
            MCmove_pivot(sweep);
        }
        // t_pivot.~Timer();

        Timer t_crankshaft("cranking", profiling_on);
        for (int j = 0; j < n_crank; j++) {
            MCmove_crankshaft();
        }
        // t_crankshaft.~Timer();

        Timer t_translation("translating", profiling_on);
        for (int j = 0; j < n_trans; j++) {
            MCmove_translate();
        }
        // t_translation.~Timer();

        if (gridmove_on)
            MCmove_grid();

        Timer t_displace("displacing", profiling_on);
        for (int j = 0; j < n_disp; j++) {
            MCmove_displace();
        }
        // t_displace.~Timer();

        Timer t_rotation("Rotating", profiling_on);
        for (int j = 0; j < n_rot; j++) {
            MCmove_rotate();
        }
        // t_rotation.~Timer();

        if (sweep % dump_frequency == 0 || sweep == nSweeps) {
            analytics.log(sweep);
            saveXyz();
            saveContacts(sweep);
            printAcceptanceRates(sweep);
        }

        if (sweep % dump_stats_frequency == 0) {
            saveEnergy(sweep);
            updateContacts();
            saveObservables(sweep);
        }
    }

    checkConsistency();
    std::cout << "overall acceptance rate: "
              << (float)acc / (nSweeps * nSteps) * 100.0 << "%" << std::endl;
}

void Sim::printAcceptanceRates(int sweep) {
    std::cout << "acceptance rate: "
              << (float)acc / ((sweep + 1) * nSteps) * 100.0 << "%"
              << std::endl;
    if (displacement_on)
        std::cout << "disp: " << (float)acc_disp / (sweep * n_disp) * 100
                  << "% \t";
    if (translation_on)
        std::cout << "trans: " << (float)acc_trans / (sweep * n_trans) * 100
                  << "% \t";
    if (crankshaft_on)
        std::cout << "crank: " << (float)acc_crank / (sweep * n_crank) * 100
                  << "% \t";
    if (pivot_on)
        std::cout << "pivot: " << (float)acc_pivot / (sweep * n_pivot) * 100
                  << "% \t";
    if (rotate_on)
        std::cout << "rot: " << (float)acc_rot / (sweep * n_rot) * 100
                  << "% \t";
    // std::cout << "cellcount: " << grid.cellCount();
    std::cout << std::endl;
}

void Sim::checkConsistency() { assert(grid.checkCellConsistency(nbeads)); }

void Sim::MCmove_displace() {
    Timer t_displacemove("Displacement move", profiling_on);
    // pick random particle
    int o = floor(beads.size() * rng->uniform());

    // copy old info (don't forget orientation, etc)
    Cell *old_cell = grid.getCell(beads[o]);

    // proposed displacement
    Eigen::RowVector3d displacement;
    displacement = step_disp * unit_vec(displacement);

    Eigen::RowVector3d new_location = beads[o].r + displacement;

    // check if exited the simulation box, if so reject the move
    if (outside_boundary(new_location)) {
        return;
    }

    // update grid
    Cell *new_cell = grid.getCell(new_location);

    std::unordered_set<Cell *> flagged_cells;
    flagged_cells.insert(old_cell);
    flagged_cells.insert(new_cell);

    double Uold = getTotalEnergy(o, o, flagged_cells);

    // move
    beads[o].r = new_location;

    // check if moved grid into new grid cell, update grid
    if (new_cell != old_cell) {
        new_cell->moveIn(&beads[o]);
        old_cell->moveOut(&beads[o]);
    }

    double Unew = getTotalEnergy(o, o, flagged_cells);

    if (rng->uniform() < exp(Uold - Unew)) {
        // move accepted
        acc += 1;
        acc_disp += 1;
        analytics.nbeads_moved += 1;
    } else {
        // move rejected
        beads[o].r -= displacement;
        new_cell->moveOut(&beads[o]);
        old_cell->moveIn(&beads[o]);
    }
}

void Sim::MCmove_translate() {
    // Timer t_trans("Translate", print_trans);
    // Timer t_setup("setup", print_trans);

    // select last bead from two-sided exponential distribution around first
    int first = floor(beads.size() * rng->uniform());
    int last = -1;
    while (last < 0 || last >= nbeads) {
        last = std::round(
            randomExp(first, exp_decay)); // does this obey detailed balance?
    }

    if (last < first) {
        std::swap(first, last);
    } // swap first and last to ensure last > first
    // generate displacement vector with magnitude step_trans
    Eigen::RowVector3d displacement;
    displacement = step_trans * unit_vec(displacement);

    // memory storage objects
    std::unordered_set<Cell *> flagged_cells;
    std::unordered_map<int, std::pair<Cell *, Cell *>>
        bead_swaps; // index of beads that swapped cell locations

    flagged_cells.reserve(last - first);
    bead_swaps.reserve(last - first);

    Cell *old_cell_tmp;
    Cell *new_cell_tmp;
    Eigen::RowVector3d new_loc;

    // t_setup.~Timer();

    // execute move
    // Timer t_flag("Flag cells", print_trans);
    // flag appropriate cells for energy calculation and find beads that swapped
    // cells
    for (int i = first; i <= last; i++) {
        new_loc = beads[i].r + displacement;
        if (outside_boundary(new_loc)) {
            return;
        }

        old_cell_tmp = grid.getCell(beads[i]);
        new_cell_tmp = grid.getCell(new_loc);

        if (new_cell_tmp != old_cell_tmp) {
            bead_swaps[i] = std::make_pair(old_cell_tmp, new_cell_tmp);
            flagged_cells.insert(new_cell_tmp);
            flagged_cells.insert(old_cell_tmp);
        }
    }
    // t_flag.~Timer();

    // Timer t_uold("Uold", print_trans);
    // std::cout << "Beads: " << last-first << " Cells: " <<
    // flagged_cells.size() << std::endl;
    double Uold = getTotalEnergy(first, last, flagged_cells);
    // t_uold.~Timer();

    // Timer t_disp("Displacement", print_trans);
    for (int i = first; i <= last; i++) {
        beads[i].r += displacement;
    }
    // t_disp.~Timer();

    // Timer t_swap("Bead Swaps", print_trans);
    //  update grid <bead index,   <old cell , new cell>>
    // for(std::pair<int, std::pair<Cell*, Cell*>> &x : bead_swaps)
    for (auto const &x : bead_swaps) {
        x.second.first->moveOut(&beads[x.first]);
        x.second.second->moveIn(&beads[x.first]);
    }
    // t_swap.~Timer();

    // Timer t_unew("Unew", print_trans);
    double Unew = getTotalEnergy(first, last, flagged_cells);
    // t_unew.~Timer();

    if (rng->uniform() < exp(Uold - Unew)) {
        // move accepted
        acc += 1;
        acc_trans += 1;
        analytics.nbeads_moved += (last - first);
    } else {
        // move rejected
        for (int i = first; i <= last; i++) {
            beads[i].r -= displacement; // restore particle positions
        }

        if (bead_swaps.size() > 0) {
            // restore old grid populations
            for (auto const &x : bead_swaps) {
                x.second.first->moveIn(&beads[x.first]);
                x.second.second->moveOut(&beads[x.first]);
            }
        }
    }
}

void Sim::MCmove_crankshaft() {
    // index of first bead through index of last bead represent all beads that
    // MOVE in the crankshaft axis of crankshaft motion is between beads
    // (first-1) and (last+1), because those two do not move.

    // select a first bead at random (cannot be first or last bead)
    int first = -1;
    while (first < 1 || first > nbeads - 2) {
        first = floor(beads.size() * rng->uniform());
    }

    // choose second bead from two-sided exponential distribution around first
    int last = -1;
    while (last < 1 || last > nbeads - 2) {
        last = std::round(randomExp(first, exp_decay_crank));
    }

    // swap first and last to ensure last > first
    if (last < first) {
        std::swap(first, last);
    }

    // compute axis of rotation, create quaternion
    Eigen::RowVector3d axis = beads[last + 1].r - beads[first - 1].r;
    double angle =
        step_crank *
        (rng->uniform() - 0.5); // random symmtric angle in cone size step_crank
    Eigen::Quaterniond du;
    du = Eigen::AngleAxisd(
        angle, axis.normalized()); // object representing this rotation

    // memory storage objects
    std::vector<Eigen::RowVector3d> old_positions;
    std::vector<Eigen::RowVector3d> old_orientations;
    std::unordered_set<Cell *> flagged_cells;
    std::unordered_map<int, std::pair<Cell *, Cell *>> bead_swaps;

    old_positions.reserve(last - first);
    old_orientations.reserve(last - first);
    flagged_cells.reserve(last - first);
    bead_swaps.reserve(last - first);

    Cell *old_cell_tmp;
    Cell *new_cell_tmp;

    // execute move
    try {
        double Uold = 0;
        if (bonded_on)
            Uold += getBondedEnergy(first, last);

        for (int i = first; i <= last; i++) {
            // save old configuration
            // --------------------- can this be done more efficiently?
            // ------------------------------------------
            old_positions.push_back(beads[i].r);
            old_orientations.push_back(beads[i].u);

            // step to new configuration, but don't update grid yet (going to
            // check if in bounds first)
            beads[i].r = du * (beads[i].r - beads[first - 1].r) +
                         beads[first - 1].r.transpose();
            beads[i].u = du * beads[i].u;
        }

        // reject if moved out of simulation box, need to restore old bead
        // positions
        for (int i = first; i <= last; i++) {
            if (outside_boundary(beads[i].r)) {
                throw "exited simulation box";
            }
        }

        // flag cells and bead swaps, but do not update the grid
        for (int i = first; i <= last; i++) {
            new_cell_tmp = grid.getCell(beads[i]);
            old_cell_tmp = grid.getCell(old_positions[i - first]);

            if (new_cell_tmp != old_cell_tmp) {
                bead_swaps[i] = std::make_pair(old_cell_tmp, new_cell_tmp);
                flagged_cells.insert(new_cell_tmp);
                flagged_cells.insert(old_cell_tmp);
            }
        }

        // calculate old nonbonded energy based on flagged cells
        if (nonbonded_on)
            Uold += getNonBondedEnergy(flagged_cells);

        // Update grid
        // for(std::pair<int, std::pair<Cell*, Cell*>> &x : bead_swaps)
        for (auto const &x : bead_swaps) {
            x.second.first->moveOut(&beads[x.first]); // out of the old cell
            x.second.second->moveIn(&beads[x.first]); // in to the new cell
        }

        double Unew = getTotalEnergy(first, last, flagged_cells);

        if (rng->uniform() < exp(Uold - Unew)) {
            // std::cout << "Accepted"<< std::endl;
            acc += 1;
            acc_crank += 1;
            analytics.nbeads_moved += (last - first);
        } else {
            // std::cout << "Rejected" << std::endl;
            throw "rejected";
        }
    }
    // REJECTION CASES -- restore old conditions
    catch (const char *msg) {
        // restore particle positions
        for (std::size_t i = 0; i < old_positions.size(); i++) {
            beads[first + i].r = old_positions[i];
            beads[first + i].u = old_orientations[i];
        }

        // restore grid allocations
        if (bead_swaps.size() > 0) {
            // for(std::pair<int, std::pair<Cell*, Cell*>> &x : bead_swaps)
            for (auto const &x : bead_swaps) {
                x.second.first->moveIn(&beads[x.first]); // back in to the old
                x.second.second->moveOut(
                    &beads[x.first]); // back out of the new
            }
        }
    }
}

void Sim::MCmove_rotate() {
    // pick random particle
    int o = floor(beads.size() * rng->uniform());

    // old configuration
    Eigen::RowVector3d u_old = beads[o].u;

    double Uold =
        getBondedEnergy(o, o); // only need bonded energy for roation moves

    // step to new configuration
    double angle =
        step_rot *
        (rng->uniform() - 0.5); // random symmtric angle in cone size step_rot
    Eigen::RowVector3d axis;    // random axis
    axis = unit_vec(axis);
    Eigen::Quaterniond du{
        Eigen::AngleAxisd(angle, axis)}; // object representing this rotation

    beads[o].u = du * beads[o].u;

    double Unew = getBondedEnergy(o, o);

    if (rng->uniform() < exp(Uold - Unew)) {
        // std::cout << "Accepted"<< std::endl;
        acc += 1;
        acc_rot += 1;
    } else {
        // std::cout << "Rejected" << std::endl;
        beads[o].u = u_old;
    }
}

void Sim::MCmove_pivot(int sweep) {
    // index terminology:
    // pivot === the bead being pivoted around, but not itself moved
    // [first,last] === the interval of beads physically moved in the pivot
    // depending on which end the pivot is executed, this is either [0, pivot-1]
    // or [pivot+1, nbeads-1]
    // The only bonds affected are pivot-1 OR pivot;

    // chose one end of the polymer and a pivot bead
    int end = (nbeads - 1) *
              std::round(rng->uniform()); //  either first bead or last bead

    // end = nbeads-1;
    //  pick second bead according to single-sided exponential distribution away
    //  from end
    int length;
    do {
        length = std::abs(std::round(randomExp(
            0, exp_decay_pivot))); // length down the polymer from the end
    } while (length < 1 || length > nbeads - 1);

    int pivot = (end == 0) ? length : (nbeads - 1 - length);

    int first = (pivot < end) ? pivot + 1 : end;
    int last = (pivot < end) ? end : pivot - 1;

    // rotation objects
    double angle =
        step_pivot *
        (rng->uniform() - 0.5); // random symmtric angle in cone size step_pivot
    Eigen::RowVector3d axis;    // random axis
    axis = unit_vec(axis);
    Eigen::Quaterniond du{
        Eigen::AngleAxisd(angle, axis)}; // object representing this rotation

    // memory storage objects
    std::vector<Eigen::RowVector3d> old_positions;
    std::vector<Eigen::RowVector3d> old_orientations;
    std::unordered_set<Cell *> flagged_cells;
    std::unordered_map<int, std::pair<Cell *, Cell *>> bead_swaps;

    Cell *old_cell_tmp;
    Cell *new_cell_tmp;

    // execute move
    try {
        double Uold = 0;
        if (bonded_on)
            Uold += getBondedEnergy(pivot - 1, pivot);

        for (int i = first; i <= last; i++) {
            // save old positions
            old_positions.push_back(beads[i].r);
            old_orientations.push_back(beads[i].u);

            // step to new configuration, but don't update grid yet (going to
            // check if in bounds first)
            beads[i].r =
                du * (beads[i].r - beads[pivot].r) + beads[pivot].r.transpose();
            beads[i].u = du * beads[i].u;
        }

        // reject if moved out of simulation box
        for (int i = first; i <= last; i++) {
            if (outside_boundary(beads[i].r)) {
                throw "exited simulation box";
            }
        }

        // flag cells and bead swaps, but do not update the grid
        for (int i = first; i <= last; i++) {
            new_cell_tmp = grid.getCell(beads[i]);
            old_cell_tmp = grid.getCell(old_positions[i - first]);

            if (new_cell_tmp != old_cell_tmp) {
                bead_swaps[i] = std::make_pair(old_cell_tmp, new_cell_tmp);
                flagged_cells.insert(old_cell_tmp);
                flagged_cells.insert(new_cell_tmp);
            }
        }

        // calculate old nonbonded energy based on flagged cells
        if (nonbonded_on)
            Uold += getNonBondedEnergy(flagged_cells);

        // Update grid
        // for(std::pair<int, std::pair<Cell*, Cell*>> &x : bead_swaps)
        for (auto const &x : bead_swaps) {
            x.second.first->moveOut(&beads[x.first]); // out of the old cell
            x.second.second->moveIn(&beads[x.first]); // in to the new cell
        }

        double Unew = getTotalEnergy(pivot - 1, pivot, flagged_cells);

        if (rng->uniform() < exp(Uold - Unew)) {
            acc += 1;
            acc_pivot += 1;
            analytics.nbeads_moved += (last - first);
        } else {
            throw "rejected";
        }
    }
    // REJECTION CASES -- restore old conditions
    catch (const char *msg) {
        // restore particle positions
        for (std::size_t i = 0; i < old_positions.size(); i++) {
            beads[first + i].r = old_positions[i];
            beads[first + i].u = old_orientations[i];
        }

        // restore bead allocations
        if (bead_swaps.size() > 0) {
            // for(std::pair<int, std::pair<Cell*, Cell*>> &x : bead_swaps)
            for (auto const &x : bead_swaps) {
                x.second.first->moveIn(&beads[x.first]); // back in to the old
                x.second.second->moveOut(
                    &beads[x.first]); // back out of the new
            }
        }
    }
}

void Sim::MCmove_grid() {
    // not really a MC move (metropolis criterion doesn't apply)
    // don't need to choose active cells; they are chosen at the beginning of
    // the simulation to include all cells that could possibly include
    // particles.
    bool flag = true;
    double U;
    while (flag) {
        Eigen::RowVector3d displacement;
        Eigen::RowVector3d old_origin = grid.origin;

        displacement = step_grid * unit_vec(displacement);
        grid.origin += displacement;

        // periodic boundary conditions: move inside volume bounded by (-delta,
        // -delta, -delta) and (0,0,0)
        grid.origin(0) -= std::ceil(grid.origin(0) / grid.delta) * grid.delta;
        grid.origin(1) -= std::ceil(grid.origin(1) / grid.delta) * grid.delta;
        grid.origin(2) -= std::ceil(grid.origin(2) / grid.delta) * grid.delta;

        // remesh beads.
        grid.meshBeads(beads);

        U = getNonBondedEnergy(grid.active_cells);

        // don't accept if move violates density maximum
        if (U < 9999999999) {
            // std::cout << "passed grid move" << std::endl;
            flag = false;
            // grid.setActiveCells();
        } else {
            // std::cout << "failed grid move" << std::endl;
            flag = false;
            grid.origin = old_origin;
            grid.meshBeads(beads); // remesh back with old origin
        }
    }
}

// write bead coordinates and bead types to .xyz file
// format: {id} {x,y,z} {bead_types}
void Sim::saveXyz() {
    xyz_out = fopen(xyz_out_filename.c_str(), "a");
    fprintf(xyz_out, "%d\n", nbeads);
    fprintf(xyz_out, "atoms\n");

    for (Bead bead : beads) {
        fprintf(xyz_out, "%d\t %lf\t %lf\t %lf\t", bead.id, bead.r(0),
                bead.r(1), bead.r(2));

        if (plaid_on) {
            for (int i = 0; i < nspecies; i++) {
                fprintf(xyz_out, "%f\t", bead.d[i]);
            }
        }

        fprintf(xyz_out, "\n");
    }
    fclose(xyz_out);
}

void Sim::saveEnergy(int sweep) {
    double bonded = 0;
    bonded = bonded_on ? getAllBondedEnergy() : 0;
    double plaid = 0;
    if (plaid_on) {
        if (smatrix_on) {
            plaid = grid.SmatrixEnergy(grid.active_cells, smatrix);
        } else if (ematrix_on) {
            plaid = grid.EmatrixEnergy(grid.active_cells, ematrix);
        } else {
            plaid = grid.energy(grid.active_cells, chis);
        }
    }
    double diagonal = 0;
    if (diagonal_on) {
        if (dmatrix_on) {
            diagonal = grid.DmatrixEnergy(grid.active_cells, dmatrix);
        } else {
            diagonal = grid.diagEnergy(grid.active_cells, diag_chis);
        }
    }

    double boundary = 0;
    boundary =
        boundary_attract_on ? getJustBoundaryEnergy(grid.active_cells) : 0;
    energy_out = fopen(energy_out_filename.c_str(), "a");
    fprintf(energy_out, "%d\t %lf\t %lf\t %lf\t %lf\t %lf\n", sweep, bonded,
            plaid, diagonal, boundary, bonded + plaid + diagonal + boundary);
    fclose(energy_out);
}

void Sim::saveObservables(int sweep) {
    // TODO phis are not updated unless energy function is called
    // leads to error if dumping observables after a rejected move;
    // beads are returned to their original state and typenums is updated
    // but cell.phis is not
    double U = grid.energy(grid.active_cells, chis); // to update phis in cells
    if (plaid_on) {
        obs_out = fopen(obs_out_filename.c_str(), "a");
        fprintf(obs_out, "%d", sweep);

        for (int i = 0; i < nspecies; i++) {
            for (int j = i; j < nspecies; j++) {
                double ij_contacts = grid.get_ij_Contacts(i, j);
                fprintf(obs_out, "\t%lf", ij_contacts);
            }
        }

        fprintf(obs_out, "\n");
        fclose(obs_out);
    }

    if (constant_chi_on) {
        obs_out = fopen(constant_obs_out_filename.c_str(), "a");
        fprintf(obs_out, "%d", sweep);

        double contacts = grid.getContacts();
        fprintf(obs_out, "\t%lf", contacts);

        fprintf(obs_out, "\n");
        fclose(obs_out);
    }

    if (diagonal_on || dmatrix_on)
    // if dmatrix_on and (smatrix_on or ematrix_on), diagonal_on will be set to
    // False for computational efficiency
    {
        double Udiag = grid.diagEnergy(
            grid.active_cells, diag_chis); // to update phis_diag? jan 28-2022
        diag_obs_out = fopen(diag_obs_out_filename.c_str(), "a");
        fprintf(diag_obs_out, "%d", sweep);

        std::vector<double> diag_obs(diag_chis.size(), 0.0);
        grid.getDiagObs(diag_obs);

        for (auto &e : diag_obs) {
            fprintf(diag_obs_out, "\t%lf", e);
        }

        fprintf(diag_obs_out, "\n");
        fclose(diag_obs_out);
    }

    if (dump_density) {
        density_out = fopen(density_out_filename.c_str(), "a");
        fprintf(density_out, "%d", sweep);

        double avg_density = 0;
        int i = 0;
        for (Cell *cell : grid.active_cells) {
            i++;
            // fprintf(density_out, " %lf", cell->phis[0]);
            avg_density += cell->phis[0];
        }
        avg_density /= i;
        fprintf(density_out, " %lf\n", avg_density);
        fclose(density_out);
    }

    extra_out = fopen(extra_out_filename.c_str(), "a");
    double phi_c = grid.getChromatinVolfrac();
    double phi_c2 = grid.getChromatinVolfrac2();
    double phi_cD = grid.getChromatinVolfracD();
    fprintf(extra_out, "%.8f %.8f %.8f\n", phi_c, phi_c2, phi_cD);
    fclose(extra_out);
}

// write contact map to file
void Sim::saveContacts(int sweep) {
    if (track_contactmap) {
        // outputs new contact map, doesn't override
        contact_map_filename = "./" + data_out_filename + "/contacts" +
                               std::to_string(sweep) + ".txt";
    } else {
        // overwrites contact file with most current values
        contact_map_filename = "./" + data_out_filename + "/contacts.txt";
    }

    std::ofstream contactsOutFile(contact_map_filename);
    for (const auto &row : contact_map) {
        for (std::size_t i = 0; i < row.size(); i++) {
            if (i == row.size() - 1) {
                contactsOutFile << row[i];
            } else {
                contactsOutFile << row[i] << " ";
            }
        }
        contactsOutFile << "\n";
    }
}

void Sim::setupSmatrix() {
    std::ifstream smatrixfile(smatrix_filename);
    smatrix.resize(nbeads, nbeads);

    if (smatrixfile.good()) {
        std::cout << smatrix_filename << "smatrix_filename is good\n";
        for (int i = 0; i < nbeads; i++) {
            for (int j = 0; j < nbeads; j++) {
                smatrixfile >> smatrix(i, j);
            }
        }
    } else {
        std::cout << "smatrix_filename (" << smatrix_filename
                  << ") does not exist or cannot be opened\n";
        Eigen::MatrixXd psi;

        // need to define psi
        psi = Eigen::MatrixXd::Zero(nbeads, nspecies);
        for (int i = 0; i < nbeads; i++) {
            for (int k = 0; k < nspecies; k++) {
                psi(i, k) = beads[i].d[k];
            }
        }
        // ensure chis are triu
        Eigen::MatrixXd chis_triu;
        chis_triu = chis.triangularView<Eigen::Upper>();

        // try and create smatrix from chi and psi
        smatrix = psi * chis * psi.transpose();
    }

    if (dmatrix_on) {
        Eigen::MatrixXd dmatrix_diag = dmatrix.diagonal().asDiagonal();
        dmatrix += dmatrix_diag;
        smatrix += dmatrix;
        diagonal_on = false;
    }

    std::cout << "loaded Smatrix, first element:" << smatrix(0, 0) << std::endl;
}

void Sim::setupEmatrix() {
    std::ifstream ematrixfile(ematrix_filename);
    ematrix.resize(nbeads, nbeads);

    if (ematrixfile.good()) {
        for (int i = 0; i < nbeads; i++) {
            for (int j = 0; j < nbeads; j++) {
                ematrixfile >> ematrix(i, j);
            }
        }
    } else {
        // first get smatrix
        Sim::setupSmatrix();

        Eigen::MatrixXd left = smatrix + smatrix.transpose();
        Eigen::MatrixXd right = smatrix.diagonal().asDiagonal();

        ematrix = left - right;
    }
    // can delete dmatrix and smatrix
    std::cout << "loaded Ematrix, first element:" << ematrix(0, 0) << std::endl;
}

void Sim::setupDmatrix() {
    std::ifstream dmatrixfile(dmatrix_filename);
    dmatrix.resize(nbeads, nbeads);

    if (dmatrixfile.good()) {
        std::cout << dmatrix_filename << " opened\n";
        for (int i = 0; i < nbeads; i++) {
            for (int j = 0; j < nbeads; j++) {
                dmatrixfile >> dmatrix(i, j);
            }
        }
    } else {
        std::cout << "dmatrix_filename (" << dmatrix_filename
                  << ") does not exist or cannot be opened\n";

        // try and create dmatrix from diag_chis
        for (int i = 0; i < nbeads; i++) {
            for (int j = 0; j < nbeads; j++) {
                int d = std::abs(i - j);
                if ((d <= Cell::diag_cutoff) && (d >= Cell::diag_start)) {
                    d -= Cell::diag_start; // TODO check that this works for
                                           // non-zero diag_start
                    int d_index = Cell::binDiagonal(d);
                    dmatrix(i, j) = diag_chis[d_index];
                } else {
                    dmatrix(i, j) = 0;
                }
            }
        }
    }
    std::cout << "loaded Dmatrix, first element:" << dmatrix(0, 0) << std::endl;
}
