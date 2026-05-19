#pragma once

struct Parameters
{
    // =========================================
    // Solver
    // =========================================
    static constexpr int Imax = 1280;
    static constexpr int iter_max = 500000;

    // =========================================
    // Physics
    // =========================================
    static constexpr double R  = 286.9865;
    static constexpr double k  = 1.4;
    static constexpr double cv = 718.0;

    static constexpr double p0 = 25000.0;
    static constexpr double T0 = 300.0;
    static constexpr double pb = 15000.0;

    // =========================================
    // Numerics
    // =========================================
    static constexpr double CFL = 0.4;
    static constexpr double Cx  = 5.0;

    // =========================================
    // Geometry
    // =========================================
    static constexpr double x0 = -0.25;
    static constexpr double x1 =  1.03;
    static constexpr double dx = (x1 - x0) / Imax;
    static constexpr double A_throat = 0.03150;
};