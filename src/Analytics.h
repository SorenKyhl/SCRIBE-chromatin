#pragma once

#include <iostream>
#include <chrono>

// keep track of total time elapsed,
// time of last block,
// total beads moved,
// beads moved last block
class Analytics {
public:
    Analytics() : beads_moved_last_sweep{0}, nbeads_moved{0} {}

    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> laststop;
    std::chrono::seconds totaltime, blocktime;
    int nbeads_moved, beads_moved_last_sweep;

    void startTimer() {
        start = std::chrono::high_resolution_clock::now();
        laststop = std::chrono::high_resolution_clock::now();
    }

    void lapTimer() {
        auto stop = std::chrono::high_resolution_clock::now();
        totaltime =
            std::chrono::duration_cast<std::chrono::seconds>(stop - start);
        blocktime =
            std::chrono::duration_cast<std::chrono::seconds>(stop - laststop);
        laststop = stop;
    }

    void log(int sweep) {
        lapTimer();
        int beads_moved_this_block = nbeads_moved - beads_moved_last_sweep;
        beads_moved_last_sweep = nbeads_moved;

        std::cout << "------ Sweep number " << sweep << " ------\n";
        std::cout << "elapsed: " << totaltime.count() << "sec \t|\t";
        if (totaltime.count() > 0) {
            std::cout << sweep / totaltime.count() << " sweep/sec \n";
        }
        std::cout << "Total beads moved " << nbeads_moved << " beads \t|\t ";
        if (blocktime.count() > 0) {
            std::cout << "Rate (this block)"
                      << beads_moved_this_block / blocktime.count()
                      << " beads/s \n";
        }
    }
};
