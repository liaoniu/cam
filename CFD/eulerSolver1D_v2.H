#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;
using mat = vector<array<double, 3>>;   // u[i][var] = var[i]
using arr = array<double, 3>;

// Global Variables

extern double x0, x1, T, Gamma;
extern int nxCells;
extern double dx;
extern int nVar;

// EoS
double computePressureFromEoS(const double& rho, const double& e);
double computeInternalEnergyFromEoS(const double& rho, double& p);
double computeSoundSpeedFromEoS(const double& rho, const double& p);

// Boundary Condition
void transmissiveBC(mat &u);

// Change between Primitive and Conserved
arr indi_primitiveToConserved(const arr &u_pri);
arr indi_conservedToPrimitive(const arr &u_con);
mat conservedToPrimitive(const mat &u);
mat primitiveToConserved(const mat &u);

// Other Solve Routine
double calc_dt(const mat &u_con);
arr f(const arr &u);

// exact Solver Routine
tuple<double, double> fK_and_dfK(double p_star, double rhoK, double vK, double pK, double csK);
double compute_p_star(const arr &uL, const arr &uR, double TOL);
arr exactFindUmid(const arr &uL, const arr &uR);
arr getXFlux(const arr &uL, const arr &uR, double &dt, string method); //u_left and u_right
void update(mat &u, mat &flux, double &dt, string method);

// Data Storing Routine
void storeData(mat &u);




// EoS Routine

double computePressureFromEoS(const double& rho, const double& e)
{
	// Ideal Gas Eos

	double p = rho*e*(Gamma-1.0);

	// Add here other EoS ....

	return p;
}

double computeInternalEnergyFromEoS(const double& rho, double& p)
{
	// Ideal Gas Eos

	double e = p/(rho*(Gamma-1.0));

	// Add here other EoS ....

	return e;
}

double computeSoundSpeedFromEoS(const double& rho, const double& p)
{
	// Ideal Gas Eos

	double a = sqrt(Gamma*p/rho);

	// Add here other EoS ....

	return a;
}


// Boundary Condition

void transmissiveBC(mat &u)
{
    for (int var = 0; var != nVar; var++){
        u[0][var] = u[1][var];
        u[nxCells+1][var] = u[nxCells][var];
    }
}


// Variables inter change

arr indi_primitiveToConserved(const arr &u_pri)
{
    arr u_con = u_pri;
    double rho, v, p;
    tie(rho, v, p) = tie(u_pri[0], u_pri[1], u_pri[2]);
    double rhoV = rho*v;
    double epsilon = computeInternalEnergyFromEoS(rho, p);
    double E = rho*epsilon + rho*(v*v)/2;
    tie(u_con[0], u_con[1], u_con[2]) = tie(rho, rhoV, E);
    return u_con;
}

arr indi_conservedToPrimitive(const arr &u_con)
{
    arr u_pri = u_con;
    double rho, rhoV, E;
    tie(rho, rhoV, E) = tie(u_con[0], u_con[1], u_con[2]);
    double v = rhoV/rho;
    double epsilon = E/rho - 0.5*(v*v);
    double p = computePressureFromEoS(rho, epsilon);
    tie(u_pri[0], u_pri[1], u_pri[2]) = tie(rho, v, p);
    return u_pri;
}

mat conservedToPrimitive(const mat &u)
{
    //u would not be changed
    mat u_pri = u;
    for (int i = 0; i != u.size(); i++){
        u_pri[i] = indi_conservedToPrimitive(u[i]);
    }
    return u_pri;
}

mat primitiveToConserved(const mat &u)
{
    //u would not be changed
    mat u_con = u;
    for (int i = 0; i != u.size(); i++){
        u_con[i] = indi_primitiveToConserved(u[i]);
    }
    return u_con;
}


// Other Routine

double calc_dt(const mat &u_con)
{
    double delta = dx;
    double C = 0.8;
    double a_max = 0;
    int m = u_con.size();
    for (int i = 0; i != m; i++){
            double rho, v, p;
            arr u_pri_indi = indi_conservedToPrimitive(u_con[i]);
            tie(rho, v, p) = tie(u_pri_indi[0], u_pri_indi[1], u_pri_indi[2]);
            double a = sqrt(v*v) + sqrt(Gamma*p/rho);
            a_max = ((a_max < a) ? a : a_max);
    }
    double dt = C*delta/a_max;
    return dt;
}

arr f(const arr &u)
{
    // u is conserved
    double rho, v, p;
    double rhoV, E;
    arr u_pri = indi_conservedToPrimitive(u);
    tie(rho, rhoV, E) = tie(u[0], u[1], u[2]);
    tie(ignore, v, p) = tie(u_pri[0], u_pri[1], u_pri[2]);
    arr f_res;
    f_res[0] = rhoV;
    f_res[1] = rhoV*v + p;
    f_res[2] = (E + p)*v;
    return f_res;
}



arr getXFlux(const arr &uL, const arr &uR, double &dt, string method)
{
    arr flux, u_mid;
    if (method == "FORCE")
    {
        arr fL, fR, LF, RI;
        fL = f(uL);
        fR = f(uR);
        for (int var = 0; var != nVar; var++){
            LF[var] = 0.5*dx/dt*(uL[var] - uR[var]) + 0.5*(fL[var] + fR[var]);
            u_mid[var] = 0.5*(uL[var] + uR[var]) - 0.5*dt/dx*(fR[var] - fL[var]);
        }
        RI = f(u_mid);
        for (int var = 0; var != nVar; var++){
            flux[var] = 0.5*(LF[var] + RI[var]);
        }
    }
    if (method == "HLLC")
    {
        // HLLCFlux(uL, uR, flux);
    }
    if (method == "exact")
    {
        u_mid = exactFindUmid(indi_conservedToPrimitive(uL), indi_conservedToPrimitive(uR));
        u_mid = indi_primitiveToConserved(u_mid);
        flux = f(u_mid);
    }

    return flux;
}


void update(mat &u, mat &flux, double &dt, string method)
{
    for (int i = 0; i != nxCells + 1; i++)
    {
        flux[i] = getXFlux(u[i], u[i+1], dt, method);
    }

    for (int i = 1; i != nxCells + 1; i++)
    {
        for (int var = 0; var != nVar; var++)
        {
            u[i][var] = u[i][var] - dt/dx*(flux[i][var] - flux[i-1][var]);
        }
    }
    transmissiveBC(u);
}

void solver(mat &u, string method)
{
    double t = 0, dt = 0;
    mat flux(nxCells+1);
    u = primitiveToConserved(u);
    do
    {
    // for (int i = 0; i != 10; i++)
    // {
        
        dt = (calc_dt(u));
        t += dt;
        update(u, flux, dt, method);
        transmissiveBC(u);
    // }
    } while (t < T);
    
    u = conservedToPrimitive(u);
    storeData(u);
}






// Data Storing Routing

void storeData(mat &u)
{
    ofstream output("u.dat");
    for (int i = 1; i != nxCells + 1; i++)
    {
        double x = x0 + (i-0.5)*dx;
        output << x << " ";
        for (int var = 0; var != nVar; var++)
        {
            output << u[i][var] << " ";
        }
        output << endl;
    }
    output.close();
}
