#include "eulerSolver1D_main.H"


// Global Variables

double x0 = 0, x1 = 1, T = 0.012, Gamma = 1.4;
int nxCells = 1000;
double dx = (x1 - x0)/nxCells;
int nVar = 3;
double xI = 0.5;

int main()
{
    mat u0(nxCells+2);
    mat u1(nxCells+2);
    vector<double> phi0(nxCells+2);
    initializePhi(phi0);
    for (int i = 0; i != nxCells+2; i++)
    {
        double x = x0 + (i-0.5)*dx;
        u0[i][0] = ((x < 0.5) ? 1 : 0.125);
        u0[i][1] = ((x < 0.5) ? 0 : 0);
        u0[i][2] = ((x < 0.5) ? 1 : 0.1);
    }
    
    for (int i = 0; i != nxCells+2; i++)
    {
        double x = x0 + (i-0.5)*dx;
        u1[i][0] = ((x < 0.5) ? 1 : 1);
        u1[i][1] = ((x < 0.5) ? 0 : 0);
        u1[i][2] = ((x < 0.5) ? 1000 : 0.01);
    }

    solver_LevelSet2Var(u0, u1, phi0, "HLLC");
    // solver_LevelSet(u1, phi0, "HLLC");
    return 0;
}


// MAIN SOLVER

void solver_LevelSet(mat &u0, vector<double> &phi0, string method = "FORCE")
{
    double t = 0;
    double dt = 0;
    mat flux;
    mat u = u0;
    flux.resize(nxCells+1);
    u = primitiveToConserved(u);
    vector<double> phi(nxCells+2);
    initializePhi(phi);

    do{
        dt = calc_dt(u);
        t += dt;
        //Transmissive BC
        transmissiveBC(u);
        //Dimensional Split 1
        split(u, flux, dt, method);
        updatePhi(phi, u, dt);
        reinitialization(phi);
        //cout << t << endl;

    }while(t<T);

    u = conservedToPrimitive(u);
    ofstream output1("rho.dat");
    ofstream output2("v.dat");
    ofstream output3("p.dat");
    ofstream output4("phi.dat");

    for (int i = 1; i != nxCells + 1; i++){
        double x = x0 + (i-0.5)*dx;
        output1 << x << " " << u[i][0] << endl;
        output2 << x << " " << u[i][1] << endl;
        output3 << x << " " << u[i][2] << endl;
        output4 << x << " " << phi[i] << endl;
    }

    output1.close();
    output2.close();
    output3.close();
    output4.close();

}


void solver_LevelSet2Var(mat &u0, mat &u1, vector<double> &phi0, string method)
{
    double t = 0;
    double dt = 0;
    mat flux(nxCells+1);
    mat flux1(nxCells+1);
    mat u = u0;
    u = primitiveToConserved(u);
    u1 = primitiveToConserved(u1);
    vector<double> phi(nxCells+2);
    initializePhi(phi);

    do{
        dt = min(calc_dt(u), calc_dt(u1));
        t += min(dt, T-t);
        //Transmissive BC
        transmissiveBC(u);
        transmissiveBC(u1);
        //Dimensional Split 1
        split(u, flux, dt, method);
        split(u1, flux1, dt, method);
        updatePhi(phi, u, u1, dt);
        reinitialization(phi);
        //cout << t << endl;

    }while(t<T);

    u = conservedToPrimitive(u);
    u1 = conservedToPrimitive(u1);
    ofstream output1("rho.dat");
    ofstream output2("v.dat");
    ofstream output3("p.dat");
    ofstream output4("phi.dat");

    for (int i = 1; i != nxCells + 1; i++){
        double x = x0 + (i-0.5)*dx;
        output1 << x << " " << u[i][0] << " " << u1[i][0] << endl;
        output2 << x << " " << u[i][1] << " " << u1[i][1] << endl;
        output3 << x << " " << u[i][2] << " " << u1[i][2] << endl;
        output4 << x << " " << phi[i] << endl;
    }

    output1.close();
    output2.close();
    output3.close();
    output4.close();
}


// Initialize phi

void initializePhi(vector<double> &phi)
{
    for(int i = 0; i != phi.size(); i++)
    {
        double x = x0 + (i-0.5)*dx;
        phi[i] = x - xI;
    }
    transmissiveBC(phi);
}

void updatePhi(vector<double> &phi, const mat &u, double &dt)
{
    for(int i = 1; i != nxCells + 1; i++)
    {
        double v = u[i][1]/u[i][0];
        double Dphi;
        Dphi = ((v < 0) ? (phi[i+1] - phi[i])/dx : (phi[i] - phi[i-1])/dx);
        phi[i] = phi[i] - dt*v*Dphi;
    }
    transmissiveBC(phi);
}

void updatePhi(vector<double> &phi, const mat &u, const mat &u1, double &dt)
{
    for(int i = 1; i != nxCells + 1; i++)
    {
        double v;
        if (phi[i] < 0)
            v = u[i][1]/u[i][0];
        else
            v = u1[i][1]/u1[i][0];
        double Dphi;
        Dphi = ((v < 0) ? (phi[i+1] - phi[i])/dx : (phi[i] - phi[i-1])/dx);
        phi[i] = phi[i] - dt*v*Dphi;
    }
    transmissiveBC(phi);
}

void reinitialization(vector<double> &phi)
{
    int index;
    double signLeft;
    for (int i = 0; i != phi.size()-1; i++)
    {
        if (phi[i]*phi[i+1] <= 0)
        {
            index = i;
            signLeft = ((phi[i] < 0) ? -1 : 1);
            break;
        }
    }
    int pt1 = index;
    int pt2 = index + 1;
    
    while (pt1 != 0)
    {
        phi[pt1-1] = phi[pt1] + signLeft*dx;
        pt1--;
    }
    while (pt2 != phi.size())
    {
        phi[pt2+1] = phi[pt2] - signLeft*dx;
        pt2++;
    }
    transmissiveBC(phi);
}