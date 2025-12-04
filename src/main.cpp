#include <chrono>

#include "Eigen/Dense"
#include "random_mars.h"
#include "nlohmann/json.hpp"
#include "prof_timer.cpp"

#include "Bead.h"
#include "Bond.h"
#include "DSS_Bond.h"
#include "Cell.cpp"
#include "Grid.h"
#include "Sim.h"

int main(int argc, char *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();

    Sim *mySim;

    if (argc == 1) {
        mySim = new Sim();
        // log goes to stdout
    } else {
        std::string data_out_file = argv[1];
        mySim = new Sim(data_out_file);
        // log goes to {data_out_file}/log.log
    }

    mySim->run();

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Took " << duration.count() << "seconds " << std::endl;
    std::cout << "Moved " << mySim->nbeads_moved << " beads " << std::endl;
    delete mySim;
    return 0;
}
