

#include <numeric>
#include "Grid.h"
#include "cmath"

bool Grid::parallel;
bool Grid::cell_volumes;

void Grid::generate() {

    origin = {-delta / 2.0, -delta / 2.0, -delta / 2.0};

    // need L+1 grid cells per direction to fully enclose simulation boundary
    std::cout << "Initializing cells: " << std::endl;
    int cells_per_dim = L + 1;
    cells.resize(cells_per_dim);
    for (int i = 0; i < cells_per_dim; i++) {
        cells[i].resize(cells_per_dim);
        for (int j = 0; j < cells_per_dim; j++) {
            cells[i][j].resize(cells_per_dim);
            for (int k = 0; k < cells_per_dim; k++) {

                cells[i][j][k].vol = delta * delta * delta;
                cells[i][j][k].r = {i * delta, j * delta,
                                    k * delta}; // position relative to origin
                                                // cells[i][j][k].print();
            }
        }
    }
    std::cout << "Cells initialized, number: "
              << cells.size() * cells[0].size() * cells[0][0].size()
              << std::endl;
};

void Grid::setActiveCells() {
    std::cout << "Setting active cells" << std::endl;
    for (int i = 0; i <= L; i++) {
        for (int j = 0; j <= L; j++) {
            for (int k = 0; k <= L; k++) {

                if (cubic_boundary) {
                    // all cells are active at all times
                    active_cells.insert(&cells[i][j][k]);
                } else if (spherical_boundary) {
                    // std::cout << "spherical boundary!" << std::endl;
                    //  only cells within sqrt(3)*delta of radius are active, at
                    //  all times but: origin can move around, as much as a full
                    //  grid cell. so actually need 2*sqrt(3) distance of buffer
                    //  to include all possible grid cells
                    Eigen::RowVector3d cell_corner;
                    cell_corner(0) = i * delta;
                    cell_corner(1) = j * delta;
                    cell_corner(2) = k * delta;

                    Eigen::RowVector3d difference = cell_corner - sphere_center;

                    if (difference.norm() < radius + 2 * sqrt(3) * delta) {
                        active_cells.insert(&cells[i][j][k]);
                    }
                }
            }
        }
    }
    std::cout << "Active cells selected; number: " << active_cells.size()
              << std::endl;
};

void Grid::printActiveCells() {
    // prints only the active cells that contains particles
    std::cout << "Printing Active Cells: number:" << active_cells.size()
              << std::endl;
    for (Cell *cell : active_cells) {
        if (cell->contains.size() > 0)
            cell->print();
    }
};

void Grid::meshBeads(std::vector<Bead> &beads) {
    // Inserts all beads into their corresponding grid cells
    // and recomputes cell volumes
    for (std::size_t i = 0; i < cells.size(); i++) {
        for (std::size_t j = 0; j < cells[i].size(); j++) {
            for (std::size_t k = 0; k < cells[i][j].size(); k++) {
                cells[i][j][k].reset(); // initialize cells to contain no beads
            }
        }
    }

    int i, j, k;
    for (Bead &bead : beads) {
        i = floor((bead.r(0) - origin(0)) / delta);
        j = floor((bead.r(1) - origin(1)) / delta);
        k = floor((bead.r(2) - origin(2)) / delta);

        cells[i][j][k].moveIn(&bead);
    }

    if (cell_volumes)
        getCellVolumes();
};

void Grid::initialize(std::vector<Bead> &beads) {
    generate();       // create grid locations
    setActiveCells(); // activate only grid cells containing beads
    meshBeads(beads); // populate grid cells with the appropriate beads
}

void Grid::getCellVolumes() {

    // float total_vol = 0;
    for (Cell *cell : active_cells) {
        Eigen::RowVector3d s{delta, delta, delta}; // cell side lengths

        for (int alpha = 0; alpha <= 2; alpha++) {
            float left_edge = cell->r(alpha);
            float right_edge = cell->r(alpha) + delta;

            // hanging off left side
            if (left_edge < 0 && right_edge > 0) {
                s(alpha) = right_edge;
            }

            // hanging off right side
            if (left_edge < side_length && right_edge > side_length) {
                s(alpha) = side_length - left_edge;
            }

            /*
            // all the way off left
            if ( left_edge < 0 && right_edge < 0 )   {
                    s(alpha) = 0;
            }

            // all the way off left
            if ( left_edge > side_length && right_edge > side_length )   {
                    s(alpha) = 0;
            }
            */
        }

        cell->vol = s(0) * s(1) * s(2);
        // total_vol += s(0)*s(1)*s(2);
    }
    // std::cout << "total vol " << total_vol << std::endl;
    // std::cout << "side length cubed" << side_length*side_length*side_length
    // << std::endl; assert(total_vol - side_length*side_length*side_length >
    // 1e-11);
}

Cell *Grid::getCell(const Bead &bead) {
    // Returns a pointer to the cell in which a bead is located
    int i, j, k;
    i = floor((bead.r(0) - origin(0)) / delta);
    j = floor((bead.r(1) - origin(1)) / delta);
    k = floor((bead.r(2) - origin(2)) / delta);
    return &cells[i][j][k];
};

Cell *Grid::getCell(const Eigen::RowVector3d &r) {
    // Returns a pointer to the cell in which a bead is located
    int i, j, k;
    i = floor((r(0) - origin(0)) / delta);
    j = floor((r(1) - origin(1)) / delta);
    k = floor((r(2) - origin(2)) / delta);
    return &cells[i][j][k];
};

bool Grid::checkCellConsistency(int nbeads) {
    // checks to see if number of beads in the grid is the same as simulation
    // nbeads
    double cellbeads = 0;
    for (Cell *cell : active_cells) {
        cellbeads += cell->contains.size();
    }

    return (cellbeads == nbeads);
};

double Grid::densityCapEnergy(const std::unordered_set<Cell *> &flagged_cells) {
    // energy penalty due to density cap
    double U = 0;
    for (Cell *cell : flagged_cells) {
        U += cell->getDensityCapEnergy();
    }
    return U;
};

double Grid::energy(const std::unordered_set<Cell *> &flagged_cells,
                    const Eigen::MatrixXd &chis) {
    // nonbonded volume interactions
    double U = 0;
    for (Cell *cell : flagged_cells) {
        U += cell->getEnergy(chis);
    }
    return U;
};

double Grid::constantEnergy(const std::unordered_set<Cell *> &flagged_cells,
                            const double constant_chi) {
    // nonbonded volume interactions
    double U = 0;
    for (Cell *cell : flagged_cells) {
        U += cell->getConstantEnergy(constant_chi);
    }
    return U;
};

double Grid::diagEnergy(const std::unordered_set<Cell *> &flagged_cells,
                        const std::vector<double> diag_chis) {
    // nonbonded volume interactions
    double U = 0;

    if (parallel) {
        std::vector<Cell *> flagged_cells_vec(flagged_cells.begin(),
                                              flagged_cells.end());
        for (Cell *cell : flagged_cells_vec) {
            U += cell->getDiagEnergy(diag_chis);
        }
    } else {
        for (Cell *cell : flagged_cells) {
            U += cell->getDiagEnergy(diag_chis);
        }
    }

    return U;
};

double Grid::boundaryEnergy(const std::unordered_set<Cell *> &flagged_cells,
                            const double boundary_chi) {
    // nonbonded volume interactions
    double U = 0;
    for (Cell *cell : flagged_cells) {
        U += cell->getBoundaryEnergy(boundary_chi, delta);
    }
    return U;
};

double Grid::SmatrixEnergy(const std::unordered_set<Cell *> &flagged_cells,
                           const Eigen::MatrixXd &Smatrix) {
    // nonbonded volume interactions
    double U = 0;
    for (Cell *cell : flagged_cells) {
        double smatrixenergy = cell->getSmatrixEnergy(Smatrix);
        U += smatrixenergy;
    }
    return U;
};

double Grid::EmatrixEnergy(const std::unordered_set<Cell *> &flagged_cells,
                           const Eigen::MatrixXd &Ematrix) {
    // nonbonded volume interactions
    double U = 0;
    for (Cell *cell : flagged_cells) {
        double ematrixenergy = cell->getEmatrixEnergy(Ematrix);
        U += ematrixenergy;
    }
    return U;
};

double Grid::DmatrixEnergy(const std::unordered_set<Cell *> &flagged_cells,
                           const Eigen::MatrixXd &Dmatrix) {
    // nonbonded volume interactions
    double U = 0;
    for (Cell *cell : flagged_cells) {
        double dmatrixenergy = cell->getDmatrixEnergy(Dmatrix);
        U += dmatrixenergy;
    }
    return U;
};

double Grid::get_ij_Contacts(int i, int j) {
    // calculates average phi_i phi_j
    double obs = 0;
    for (Cell *cell : active_cells) {
        obs +=
            cell->typenums[i] * cell->typenums[j] * cell->beadvol / cell->vol;
    }

    // obs /= active_cells.size();
    return obs;
};

double Grid::getContacts() {
    // calculates total number of contacts
    double obs = 0;
    for (Cell *cell : active_cells) {
        obs += pow(cell->contains.size(), 2) * cell->beadvol / cell->vol;
    }

    return obs;
};

void Grid::getDiagObs(std::vector<double> &diag_obs) {
    for (Cell *cell : active_cells) {
        for (std::size_t i = 0; i < diag_obs.size(); i++) {
            // diag_obs[i] += cell->diag_phis[i] * cell->diag_phis[i] *
            // cell->vol / cell->beadvol;
            diag_obs[i] += cell->diag_phis[i] * cell->beadvol / cell->vol;
        }
    }
};

double Grid::getChromatinVolfrac2() {
    double obs = 0;
    for (Cell *cell : active_cells) {

        double phi_c = cell->contains.size() * cell->beadvol / cell->vol;
        obs += phi_c * phi_c;
    }

    obs /= active_cells.size();
    return obs;
};

double Grid::getChromatinVolfrac() {
    double obs = 0;
    for (Cell *cell : active_cells) {

        double phi_c = cell->contains.size() * cell->beadvol / cell->vol;
        obs += phi_c;
    }

    obs /= active_cells.size();
    return obs;
};

double Grid::getChromatinVolfracD() {
    double obs = 0;
    for (Cell *cell : active_cells) {

        double phi_c = cell->contains.size() * cell->beadvol / cell->vol;
        obs += (phi_c - Cell::phi_chromatin) * (phi_c - Cell::phi_chromatin);
    }

    obs /= active_cells.size();
    return obs;
};

double Grid::cellCount() {
    int sum = 0;
    int num = 0;
    for (Cell *cell : active_cells) {
        sum += cell->contains.size();
        num++;
    }

    return (double)sum / (double)num;
};
