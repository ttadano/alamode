#include "timer.h"
#include <iostream>
#include <iomanip>

using namespace ALM_NS;

Timer::Timer(ALM *alm): Pointers(alm)
{
#if defined(WIN32) || defined(_WIN32)
    QueryPerformanceCounter(&time_ref);
    QueryPerformanceFrequency(&frequency);
#else
    gettimeofday(&time_ref, NULL);
#endif
}

Timer::~Timer() {}

void Timer::reset()
{
#if defined(WIN32) || defined(_WIN32)
    QueryPerformanceCounter(&time_ref);
#else
    gettimeofday(&time_ref, NULL);
#endif
}

double const Timer::elapsed()
{
#if defined(WIN32) || defined(_WIN32)
    LARGE_INTEGER time_now;
    QueryPerformanceCounter(&time_now);
    return static_cast<double>(time_now.QuadPart - time_ref.QuadPart) / static_cast<double>(frequency.QuadPart);
#else
    timeval stopTime;
    gettimeofday(&time_now, NULL);
    return (time_now.tv_sec - time_ref.tv_sec) + (time_now.tv_usec - time_ref.tv_usec) * 1.0e-6;
#endif
}

void Timer::print_elapsed()
{
    std::cout << std::endl << "Time Elapsed: " << elapsed() << " sec." << std::endl << std::endl;
}