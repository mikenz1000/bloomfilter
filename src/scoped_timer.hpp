/*
    Simple scope-based timer that prints the amount of time between instantiation and destruction - to be used to time blocks of code
*/
#ifndef __SCOPED_TIMER_H
#define __SCOPED_TIMER_H
#include <chrono>
#include <iostream>
#include <string>

class scoped_timer {
    std::chrono::high_resolution_clock::time_point started;
    std::string description;
    int iterations;
public:
    
    /* constructs and starts the timer */
    scoped_timer(const char * description, int iterations = 0) 
        : started(std::chrono::high_resolution_clock::now()), 
        description(description),
        iterations(iterations)
    {
    }
    
    /* destroys and writes the message */
    ~scoped_timer()
    {
        std::chrono::high_resolution_clock::time_point stopped = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> seconds = std::chrono::duration_cast<std::chrono::duration<double> >(stopped-started);
        std::cout << description << " " << std::fixed << seconds.count() << "s";
//        if (iterations > 0) std::cout << std::scientific << " (" << (seconds.count()/iterations) << " x " << iterations << " its)";
        //std::cout << std::endl;
    }
};

#endif