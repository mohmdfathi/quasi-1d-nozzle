#include <vector>

#include "parameters.hpp"


struct FlowState
{
    // =========================================
    // Primitive variables
    // =========================================

    std::vector<double> rho;
    std::vector<double> u;
    std::vector<double> p;

    // =========================================
    // Constructor
    // =========================================

    FlowState();

    // =========================================
    // Methods
    // =========================================

    void initialize(double mach);

    std::vector<double> sound_speed() const;
};