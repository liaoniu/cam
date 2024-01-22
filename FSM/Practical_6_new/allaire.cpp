#include "allaire_main.H"

// Global Variables

// double x0 = 0, x1 = 1, T = 0.00023744, Gamma1 = 1.4, Gamma2 = 4.4, p_infty1 = 0, p_infty2 = 6*1e8;
double x0 = 0, x1 = 1, T = 0.00023744, Gamma1 = 4.4, Gamma2 = 1.4, p_infty1 = 6*1e8, p_infty2 = 0;
// double x0 = 0, x1 = 1, T = 0.25, Gamma1 = 1.4, Gamma2 = 1.4, p_infty1 = 0, p_infty2 = 0;
// double x0 = 0, x1 = 1, T = 0.25, Gamma1 = 1.4, Gamma2 = 1.4, p_infty1 = 0, p_infty2 = 0;

int nxCells = 1000;
double dx = (x1 - x0)/nxCells, C = 0.4;
int nVar = 5;


int main()
{
    mat u0(nxCells+2);
    for (int i = 0; i != nxCells+2; i++)
    {
        double x = x0 + (i-0.5)*dx;
        u0[i][0] = ((x < 0.7) ? 1 - 1e-6 : 1e-6);
        u0[i][1] = ((x < 0.7) ? 1000 : 1000);
        u0[i][2] = ((x < 0.7) ? 50 : 50);
        u0[i][3] = ((x < 0.7) ? 0 : 0);
        u0[i][4] = ((x < 0.7) ? 1e9 : 1e5);
    }


    // for (int i = 0; i != nxCells+2; i++)
    // {
    //     double x = x0 + (i-0.5)*dx;
    //     u0[i][0] = ((x < 0.7) ? 1e-6 : 1 - 1e-6);
    //     u0[i][1] = ((x < 0.7) ? 50 : 50);
    //     u0[i][2] = ((x < 0.7) ? 1000 : 1000);
    //     u0[i][3] = ((x < 0.7) ? 0 : 0);
    //     u0[i][4] = ((x < 0.7) ? 1e9 : 1e5);
    // }

    // for (int i = 0; i != nxCells + 2; i++)
    // {
    //     double x = x0 + (i-0.5)*dx;
    //     u0[i][0] = ((x < 0.5) ? 1e-6 : 1 - 1e-6);
    //     u0[i][1] = ((x < 0.5) ? 0.125 : 0.125);
    //     u0[i][2] = ((x < 0.5) ? 1 : 1);
    //     u0[i][3] = ((x < 0.5) ? 0 : 0);
    //     u0[i][4] = ((x < 0.5) ? 1 : 0.1);
    // }




    solver(u0, "HLLC");

    return 0;
}