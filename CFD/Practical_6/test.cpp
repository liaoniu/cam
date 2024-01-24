#include "riemann_solver.H"

double x0 = 0, x1 = 1, T = 0.25, Gamma = 1.4;
int nxCells = 2;
double dx = (x1 - x0)/nxCells;
int nVar = 3;

int main()
{
    arr uL = {1, 0, 1};
    arr uR = {0.125, 0, 0.1};
    double p_star = compute_p_star(uL, uR);
    cout << p_star << endl;
    return 0;
}