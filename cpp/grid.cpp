#include "grid.hpp"

// =============================================
// Constructor
// =============================================

NozzleGrid::NozzleGrid()
{
    const int N = Parameters::Imax + 1;

    x.resize(N);
    A.resize(N);
    dAdx.resize(N);

    // -----------------------------------------
    // Grid generation
    // -----------------------------------------

    for(int i = 0; i < N; ++i)
    {
        x[i] = Parameters::x0 + i * Parameters::dx;
    }

    // -----------------------------------------
    // Geometry
    // -----------------------------------------

    for(int i = 0; i < N; ++i)
    {
        if(x[i] <= 0.0)
        {
            dAdx[i] = -0.0318;
        }
        else
        {
            dAdx[i] = 0.0510;
        }

        A[i] = Parameters::A_throat + dAdx[i] * x[i];
    }
}