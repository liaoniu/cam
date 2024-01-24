#include "riemann_solver.H"

double x0 = 0, x1 = 1, T = 0.25, Gamma = 1.4;
int nxCells = 100;
double dx = (x1 - x0)/nxCells;
int nVar = 3;


void testSolver(mat &u, string method)
{
    double t = 0, dt = 0;
    mat flux(nxCells+1);
    u = primitiveToConserved(u);
    // do
    // {
    for (int i = 0; i != 1; i++)
    {
        
        dt = (calc_dt(u));
        t += dt;
        update(u, flux, dt, method);
        transmissiveBC(u);
    }
    // } while (t < T);
    
    u = conservedToPrimitive(u);
    storeData(u);
}


int main()
{
    mat u(nxCells+2);
    for (int i = 0; i != nxCells+2; i++)
    {
        double x = x0 + (i-0.5)*dx;
        u[i][0] = ((x < 0.5) ? 1 : 0.125);
        u[i][1] = ((x < 0.5) ? 0 : 0);
        u[i][2] = ((x < 0.5) ? 1 : 0.1);
    }

    // for (int i = 0; i != nxCells+2; i++)
    // {
    //     double x = x0 + (i-0.5)*dx;
    //     u[i][0] = ((x < 0.5) ? 0.392333 : 0.125);
    //     u[i][1] = ((x < 0.5) ? 0.982022 : 0);
    //     u[i][2] = ((x < 0.5) ? 0.336438 : 0.1);
    // }

    solver(u, "exact");
    return 0;

}