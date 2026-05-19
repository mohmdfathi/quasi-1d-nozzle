#include <iostream>
// #include "parameters.hpp"
#include "state.hpp"

int main()
{
    std::cout << "=====================================\n";
    std::cout << " Quasi-1D Nozzle Solver (C++)\n";
    std::cout << "=====================================\n";
    
    std::cout << "Grid points : "
              << Parameters::Imax + 1 << "\n";

    FlowState flow;

    flow.initialize(0.5);

    std::cout << "rho[0] = "
              << flow.rho[0]
              << "\n";
    return 0;
}