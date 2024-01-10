#include "eulerSolver1D_main.H"

// Global Variables

double x0 = 0, x1 = 1, T = 0.25, Gamma = 1.4;
int nxCells = 100;
double dx = (x1 - x0)/nxCells;
int nVar = 3;



int main()
{

    mat u0;
    u0.resize(nxCells+2);
    for (int i = 0; i != nxCells+2; i++){
        double x = x0 + (i-0.5)*dx;
        u0[i][0] = ((x < 0.5) ? 1 : 0.125);
        u0[i][1] = ((x < 0.5) ? 0 : 0);
        u0[i][2] = ((x < 0.5) ? 1 : 0.1);
    }

    solver(u0, "SLIC");
    return 0;
}
