#include <iostream>
#include "parameters.hpp"

int main()
{
    std::cout << "=====================================\n";
    std::cout << " Quasi-1D Nozzle Solver (C++)\n";
    std::cout << "=====================================\n";
    
    std::cout << "Grid points : "
              << Parameters::Imax + 1 << "\n";

    return 0;
}