#include <iostream>
#include "parameters.hpp"
#include "state.hpp"
#include "grid.hpp"

int main()
{
    std::cout << "=====================================\n";
    std::cout << " Quasi-1D Nozzle Solver (C++)\n";
    std::cout << "=====================================\n";
    
    std::cout << "Grid points : "
              << Parameters::Imax + 1 << "\n";

    FlowState flow;
    NozzleGrid grid;

    flow.initialize(0.5);

    std::cout << "rho[0] = "
              << flow.rho[0]
              << "\n";

    
    std::cout << "x[0] = "
              << grid.x[0]
              << "\n";

    return 0;
}