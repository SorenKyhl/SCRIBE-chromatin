#include "Cell.h"
#include "cmath"

// TODO: phase out
bool Cell::double_count_main_diagonal;
double Cell::beadvol;
int Cell::ntypes;
int Cell::diag_nbins;
int Cell::diag_binsize;
bool Cell::diagonal_linear;
double Cell::phi_solvent_max;
double Cell::phi_chromatin;
double Cell::kappa;
bool Cell::density_cap_on;
bool Cell::compressibility_on;
bool Cell::diag_pseudobeads_on;
bool Cell::dense_diagonal_on;
int Cell::n_small_bins;
int Cell::n_big_bins;
int Cell::small_binsize;
int Cell::big_binsize;
int Cell::diag_cutoff;
int Cell::diag_start;
bool Cell::diagonal_binning;
std::vector<int> Cell::diagonal_bin_lookup;
std::vector<int> Cell::diag_bin_of;
std::vector<int> Cell::diag_nbonds_of;
std::vector<double> Cell::diag_weight_of;

void Cell::setInteractions(int nbeads, const std::vector<double> &diag_chis,
                           bool diagonal_on) {
    // must run after all diagonal parameters are read
    diag_bin_of.assign(nbeads, -1);
    diag_nbonds_of.assign(nbeads, 0);
    diag_weight_of.assign(nbeads, 0.0);
    if (!diagonal_on) {
        return; // binning parameters are unset
    }
    for (int sep = 0; sep < nbeads; sep++) {
        if ((sep <= diag_cutoff) && (sep >= diag_start)) {
            int d = sep - diag_start; // TODO check that this works for
                                      // non-zero diag_start
            diag_bin_of[sep] = binDiagonal(d);
            if (double_count_main_diagonal) {
                diag_nbonds_of[sep] = 2; // both main and off diagonal count twice
            } else {
                diag_nbonds_of[sep] = d ? 2 : 1; // count two for all off-diagonal
            }
            diag_weight_of[sep] = diag_chis[diag_bin_of[sep]] * diag_nbonds_of[sep];
        }
    }
}

void Cell::setBeadInteractions(std::vector<Bead> &beads,
                               const Eigen::MatrixXd &chis) {
    // must run after bead types are loaded and before beads enter cells;
    // moveIn/moveOut read chi_d[0, ntypes)
    bool have_chis = chis.rows() == ntypes && chis.cols() == ntypes;
    for (Bead &bead : beads) {
        bead.chi_d.assign(ntypes, 0.0);
        if (!have_chis || (int)bead.d.size() != ntypes) {
            continue;
        }
        for (int i = 0; i < ntypes; i++) {
            double s = 0;
            for (int j = 0; j < ntypes; j++) {
                double chi = (i <= j) ? chis(i, j) : chis(j, i);
                s += (i == j ? 2 * chi : chi) * bead.d[j];
            }
            bead.chi_d[i] = s;
        }
    }
}

void Cell::print() {
    std::cout << r << "     N: " << contains.size() << std::endl;
    for (Bead *bead : contains) {
        bead->print();
    };
};

void Cell::reset() {
    // clears population trackers
    contains.clear();
    invalidateEnergy();
    std::fill(typenums.begin(), typenums.end(), 0); // DO NOT USE .clear()
    std::fill(chi_n.begin(), chi_n.end(), 0);
};

void Cell::moveIn(Bead *bead) {
    // updates local number of each type of bead
    // TODO update populations of distance ids
    contains.push_back(bead);
    invalidateEnergy();
    for (int i = 0; i < ntypes; i++) {
        typenums[i] += bead->d[i];
        chi_n[i] += bead->chi_d[i];
    }
};

void Cell::moveOut(Bead *bead) {
    // updates local number of each type of bead
    // TODO update populations of distance ids
    invalidateEnergy();
    // swap-and-pop: O(k) find, O(1) removal; order is irrelevant (this is a set)
    for (std::size_t i = 0; i < contains.size(); i++) {
        if (contains[i] == bead) {
            contains[i] = contains.back();
            contains.pop_back();
            break;
        }
    }
    for (int i = 0; i < ntypes; i++) {
        typenums[i] -= bead->d[i];
        chi_n[i] -= bead->chi_d[i];
    }
};

double Cell::getDensityCapEnergy() {
    // Density in each cell is capped at phi_solvent_max
    // otherwise, incur a large energy penalty

    float phi_beads = contains.size() * beadvol / vol;
    float phi_solvent = 1 - contains.size() * beadvol / vol;

    double U = 0;
    if (density_cap_on) {
        if (phi_solvent < phi_solvent_max) {
            // high volume fraction occurs when more than 50% of the volume is
            // occupied by beads
            U = 99999999 * phi_beads;
        }
    } else if (compressibility_on) {
        U = (phi_beads - phi_chromatin) * (phi_beads - phi_chromatin) * kappa;
    }
    return U;
};

double Cell::getEnergy(const Eigen::MatrixXd &chis) {
    // U = sum_{i<=j} chi_ij phi_i phi_j vol/beadvol with phi = n beadvol/vol,
    // i.e. beadvol/vol * sum_{i<=j} chi_ij n_i n_j = beadvol/vol * n.(S n)/2,
    // with S n maintained incrementally in chi_n (from the same chis, see
    // setBeadInteractions).
    if (energy_valid) {
        return energy_cache;
    }
    double U = 0;
    for (int i = 0; i < ntypes; i++) {
        U += typenums[i] * chi_n[i];
    }
    U *= 0.5 * beadvol / vol;
    energy_cache = U;
    energy_valid = true;
    return U;
};

double Cell::getConstantEnergy(const double constant_chi) {
    // constant nonbonded interaction between all pairs of beads
    double U = constant_chi * pow(contains.size(), 2) * beadvol / vol;

    return U;
};

double Cell::getSmatrixEnergy(const Eigen::MatrixXd &Smatrix) {
    double U = 0;

    std::vector<int> indices;
    int imax = (int)contains.size();
    for (const auto &elem : contains) {
        indices.push_back(elem->id);
    }

    assert(imax == indices.size());

    for (int i = 0; i < imax; i++) {
        for (int j = 0; j < imax; j++) {

            U += Smatrix(indices[i], indices[j]) * beadvol / vol;
        }
    }
    return U;
}

double Cell::getEmatrixEnergy(const Eigen::MatrixXd &Ematrix) {
    double U = 0;
    std::vector<int> indices;
    int imax = (int)contains.size();
    for (const auto &elem : contains) {
        indices.push_back(elem->id);
    }
    assert(imax == indices.size());
    for (int i = 0; i < imax; i++) {
        for (int j = i; j < imax; j++) {
            U += Ematrix(indices[i], indices[j]) * beadvol / vol;
        }
    }
    return U;
}

double Cell::getDmatrixEnergy(const Eigen::MatrixXd &Dmatrix) {
    double U = 0;

    std::vector<int> indices;
    int imax = (int)contains.size();
    for (const auto &elem : contains) {
        indices.push_back(elem->id);
    }

    assert(imax == indices.size());

    // this also works
    // for (int i=0; i<imax; i++)
    // {
    // 	for(int j=0; j<imax; j++)
    // 	{
    // 		if (i == j)
    // 		{
    // 			U += Dmatrix[indices[i]][indices[j]] * beadvol/vol * 2;
    // 		}
    // 		else
    // 		{
    // 			U += Dmatrix[indices[i]][indices[j]] * beadvol/vol;
    // 		}
    // 	}
    // }
    for (int i = 0; i < imax; i++) {
        for (int j = i; j < imax; j++) {
            U += Dmatrix(indices[i], indices[j]) * beadvol / vol * 2;
        }
    }
    return U;
}

int Cell::binDiagonal(int d) {
    int bin_index = -1;

    if (Cell::dense_diagonal_on) {
        int dense_cutoff = Cell::n_small_bins * Cell::small_binsize;
        // diagonal chis are binned in a dense set (small bins) from d=0 to
        // d=dense_cutoff, then a sparse set (large bins) from d=cutoff to
        // d=diag_cutoff
        if (d > dense_cutoff) {
            bin_index = Cell::n_small_bins +
                        std::floor((d - dense_cutoff) / Cell::big_binsize);
        } else {
            bin_index = std::floor(d / Cell::small_binsize);
        }
    } 
    else if (diagonal_binning){
        return diagonal_bin_lookup[d];
    } else {
        // diagonal chis are linearly spaced from d=0 to d=nbeads
        bin_index = std::floor(d / Cell::diag_binsize);
    }
    return bin_index;
}

double Cell::getDiagEnergy(const std::vector<double> &diag_chis) {
    // U = beadvol/vol * sum_bins diag_chis[b] * count[b], summed directly per
    // pair as beadvol/vol * sum_pairs diag_weight_of[|i - j|] (built from the
    // same diag_chis by setInteractions). Does not fill diag_phis; the
    // observables get those from updateDiagPhis.
    if (diag_energy_valid) {
        return diag_energy_cache;
    }
    std::size_t imax = contains.size();
    double Udiag = 0;
    // pairwise contacts -- include self-self interaction!!
    for (std::size_t i = 0; i < imax; i++) {
        int id_i = contains[i]->id;
        for (std::size_t j = i; j < imax; j++) {
            Udiag += diag_weight_of[std::abs(id_i - contains[j]->id)];
        }
    }
    diag_energy_cache = Udiag * beadvol / vol;
    diag_energy_valid = true;
    return diag_energy_cache;
};

void Cell::updateDiagPhis() {
    // per-bin pair counts, for the diagonal observables
    std::fill(diag_phis.begin(), diag_phis.end(), 0);
    std::size_t imax = contains.size();
    // count pairwise contacts  -- include self-self interaction!!
    for (std::size_t i = 0; i < imax; i++) {
        int id_i = contains[i]->id;
        for (std::size_t j = i; j < imax; j++) {
            int sep = std::abs(id_i - contains[j]->id);
            int bin = diag_bin_of[sep];
            if (bin >= 0) {
                diag_phis[bin] += diag_nbonds_of[sep];
            }
        }
    }
};

double Cell::getBoundaryEnergy(const double boundary_chi, const double delta) {
    // TODO: this is broken if the grid is moving;
    // need to check relative to origin
    double Uboundary = 0;
    for (const auto &bead : contains) {
        if (bead->r(0) < delta) {
            Uboundary += boundary_chi;
        }
    }
    return Uboundary;
};
