#pragma once

#include <vector>
#include "parameters.hpp"

struct NozzleGrid
{
    // =========================================
    // Geometry
    // =========================================

    std::vector<double> x;
    std::vector<double> A;
    std::vector<double> dAdx;

    // =========================================
    // Constructor
    // =========================================

    NozzleGrid();
};