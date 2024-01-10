#include "D:\me\cam\FSM\Practical_3\Ghost_main.H"


// Global Variables

double x0 = 0, x1 = 1, T = 0.18, GammaL = 1.4, GammaR = 1.67;
int nxCells = 400;
double dx = (x1 - x0)/nxCells;
int nVar = 3;
double xI = 0.5;
int Nghost_fluid = 3;


int main()
{
    mat u0(nxCells+2);
    mat u1(nxCells+2);
    vector<double> phi0(nxCells+2);
    initializePhi(phi0);
    for (int i = 0; i != nxCells+2; i++)
    {
        double x = x0 + (i-0.5)*dx;
        u0[i][0] = ((x < 0.25) ? 1.3765 : 1);
        u0[i][1] = ((x < 0.25) ? 0.3948 : 0);
        u0[i][2] = ((x < 0.25) ? 1.57 : 1);
    }
    
    u1 = u0;
    for (int i = 0; i != nxCells+2; i++)
    {
        double x = x0 + (i-0.5)*dx;
        if ((x > 0.4) && (x < 0.6))
        {
            u1[i][0] = 0.138;
            u1[i][1] = 0;
            u1[i][2] = 1;
        }

    }

    solver_Original_Ghost(u0, u1, phi0, GammaL, GammaR, "HLLC");
    // solver_LevelSet(u1, phi0, "HLLC");
    return 0;
}



// MAIN SOLVER

void solver_Original_Ghost(mat &u0, mat &u1, vector<double> &phi0, const double &GammaL, const double &GammaR, string method)
{
    double t = 0;
    double dt = 0;
    mat flux(nxCells+1);
    mat flux1(nxCells+1);
    u0 = primitiveToConserved(u0, GammaL);
    u1 = primitiveToConserved(u1, GammaR);
    vector<double> phi(nxCells+2);
    initializePhi_Helium(phi);

    do{
        computeGhostFluidBoundaries_Helium(u0, u1, phi, GammaL, GammaR);
        computeDomainBoundaries(u0, u1);
        dt = computeTimeStep(u0, u1, GammaL, GammaR);
        t += min(dt, T-t);
        updatePhi(phi, u0, u1, dt);
        fastSweeping(phi);
        //Dimensional Split 1
        split(u0, flux, dt, GammaL, method);
        split(u1, flux1, dt, GammaR, method);
        //cout << t << endl;

    }while(t<T);

    u0 = conservedToPrimitive(u0, GammaL);
    u1 = conservedToPrimitive(u1, GammaR);
    ofstream output1("out.dat");


    for (int i = 1; i != nxCells + 1; i++){
        double x = x0 + (i-0.5)*dx;
        output1 << x << " " << u0[i][0] << " " << u0[i][1] << " " << u0[i][2] << " " << u1[i][0] << " " << u1[i][1] << " " << u1[i][2] << " " << phi[i] << endl;
    }

    output1.close();

}

void reconstructRhoFromConstantEntropy(arr &uI, arr &uG, const double &GammaI, const double &GammaG)
{
    arr uI_pri = indi_conservedToPrimitive(uI, GammaI);
    arr uG_pri = indi_conservedToPrimitive(uG, GammaG);
    uG_pri[0] = uI_pri[0]*pow(uG_pri[2]/uI_pri[2], 1/GammaI);
    uG = indi_primitiveToConserved(uG_pri, GammaI);
}

void isobarixFix(mat &u0, mat &u1, int interface, const double &GammaL, const double &GammaR)
{
    arr u0_pri = indi_conservedToPrimitive(u0[interface], GammaL);
    arr u0_pri_minus = indi_conservedToPrimitive(u0[interface-1], GammaL);
    u0_pri[0] = u0_pri_minus[0]*pow(u0_pri[2]/u0_pri_minus[2], 1/GammaL);

    arr u1_pri = indi_conservedToPrimitive(u1[interface+1], GammaR);
    arr u1_pri_plus = indi_conservedToPrimitive(u1[interface+2], GammaR);
    u1_pri[0] = u1_pri_plus[0]*pow(u1_pri[2]/u1_pri_plus[2], 1/GammaR);


    u0[interface] = indi_primitiveToConserved(u0_pri, GammaL);
    u1[interface+1] = indi_primitiveToConserved(u1_pri, GammaR);
}



void computeGhostFluidBoundaries(mat &u0, mat &u1, vector<double> &phi, const double &GammaL, const double &GammaR)
{
    int interface = positionOfInterface(phi);
    isobarixFix(u0, u1, interface, GammaL, GammaR);
    // Copy conservative rho, rhov, and E;
    for (int i = interface; i != interface - Nghost_fluid; i--)
    {
        u1[i] = u0[i];
        reconstructRhoFromConstantEntropy(u1[interface+1], u1[i], GammaR, GammaL);
    }

    int i = interface - Nghost_fluid;
    while (i != 0)
    {
        u1[i] = u1[interface - Nghost_fluid + 1];
        i--;
    }

    for (int i = interface + 1; i != interface + 1 + Nghost_fluid; i++)
    {
        u0[i] = u1[i];
        reconstructRhoFromConstantEntropy(u0[interface], u0[i], GammaL, GammaR);
    }
    i = interface + 1 + Nghost_fluid;
    while (i != u0.size())
    {
        u0[i] = u0[interface + Nghost_fluid];
        i++;
    }
}


void computeGhostFluidBoundaries2(mat &u0, mat &u1, vector<double> &phi, const double &GammaL, const double &GammaR)
{
    int interface = positionOfInterface(phi);
    // Copy conservative rho, rhov, and E;
    for (int i = interface; i != 0; i--)
    {
        u1[i] = u0[i];
        reconstructRhoFromConstantEntropy(u1[interface+1], u1[i], GammaR, GammaL);
    }


    for (int i = interface + 1; i != u0.size(); i++)
    {
        u0[i] = u1[i];
        reconstructRhoFromConstantEntropy(u0[interface], u0[i], GammaL, GammaR);
    }

}

void computeGhostFluidBoundaries_Helium(mat &u0, mat &u1, vector<double> &phi, const double &GammaL, const double &GammaR)
{
    vector<int> interface = posOfAllInterface(phi);
    for (int i = interface[0] + 1; i != interface[0] + 1 + Nghost_fluid; i++)
    {
        u0[i] = u1[i];
        reconstructRhoFromConstantEntropy(u0[interface[0]], u0[i], GammaL, GammaR);
    }

    for (int i = interface[1]; i != interface[1] - Nghost_fluid; i--)
    {
        u0[i] = u1[i];
        reconstructRhoFromConstantEntropy(u0[interface[1]+1], u0[i], GammaL, GammaR);
    }

    for (int i = interface[0]; i != interface[0] - Nghost_fluid; i--)
    {
        u1[i] = u0[i];
        reconstructRhoFromConstantEntropy(u1[interface[0]+1], u1[i], GammaR, GammaL);
    }

    for (int i = interface[1] + 1; i != interface[1] + 1 + Nghost_fluid; i++)
    {
        u1[i] = u0[i];
        reconstructRhoFromConstantEntropy(u1[interface[1]], u1[i], GammaR, GammaL);
    }


}


void computeDomainBoundaries(mat &u0, mat &u1)
{
    transmissiveBC(u0);
    transmissiveBC(u1);
}


double computeTimeStep(mat &u0, mat &u1, const double &GammaL, const double &GammaR)
{
    double dt = calc_dt(u0, GammaL);
    double dt1 = calc_dt(u1, GammaR);
    return min(dt, dt1);
}

