#pragma once

#include <iostream>
#include <chrono>

class Timer {
public:
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<float, std::micro> duration;
    std::string title;

    bool destructor_called = false;
    bool clock_on;

    Timer(std::string msg, bool on = true) {
        start = std::chrono::high_resolution_clock::now();
        title = msg;
        clock_on = on;
    }

    ~Timer() {
        if (destructor_called == false && clock_on) {
            destructor_called = true;
            end = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            std::cout << title << " took " << duration.count()
                      << " microseconds " << std::endl;
        }
    }
};
