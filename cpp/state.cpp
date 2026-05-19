#include "state.hpp"
#include <cmath>

// =============================================
// Constructor
// =============================================

FlowState::FlowState()
{
    const int N = Parameters::Imax + 1;

    rho.resize(N);
    u.resize(N);
    p.resize(N);
}

// =============================================
// Initialize flow
// =============================================

void
FlowState::initialize(double mach)
{
    const int N = Parameters::Imax + 1;
    const double k  = Parameters::k;
    const double R  = Parameters::R;
    const double p0 = Parameters::p0;
    const double T0 = Parameters::T0;

    double fac = 1.0 + 0.5 * (k - 1.0) * mach * mach;
    double pressure = p0 / std::pow(fac, k / (k - 1.0));
    double temperature = T0 / fac;
    double density = pressure / (R * temperature);
    double velocity = mach * std::sqrt(k * R * temperature);

    for(int i = 0; i < N; ++i)
    {
        rho[i] = density;
        u[i]   = velocity;
        p[i]   = pressure;
    }
}

// =============================================
// Sound speed
// =============================================

std::vector<double>
FlowState::sound_speed() const
{
    const int N = Parameters::Imax + 1;
    const double k  = Parameters::k;

    std::vector<double> a(N);

    for(int i = 0; i < N; ++i)
    {
        a[i] = std::sqrt( k * p[i] / rho[i] );
    }

    return a;
}