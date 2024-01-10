#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;
using mat = vector<array<double, 3>>;   // u[i][var] = var[i]
using arr = array<double, 3>;

bool inCircle(double x, double y, double R);
void transmissiveBC(mat& u);
void transmissiveBC(vector<double> &phi);

arr indi_primitiveToConserved(const arr &u_pri);
arr indi_conservedToPrimitive(const arr &u_con);
mat conservedToPrimitive(const mat &u);
mat primitiveToConserved(const mat &u);
void initializePhi(vector<double> &phi);
double calc_dt(const mat &u_con);
arr f(const arr &u);
double zeta(double delta_minus, double delta_plus, string limiter);
tuple<mat, mat> dataReconstruct(const mat &u, string direction);
tuple<mat, mat> halfTimeUpdate(mat &uL, mat &uR, double dt);

void WaveEstimates(const arr& Ql, const arr& Qr, double& Sl, double& Sr);
void HLLCFlux(const arr& Ql, const arr& Qr, arr& Fhllc);
arr getXFlux(const arr &uL, const arr &uR, double &dt, string method); //u_left and u_right
mat split(mat &u, mat &flux, double dt, string method);


// LevelSet functions
void initializePhi(vector<double> &phi);
void updatePhi(vector<double> &phi, const mat &u, double &dt);
void updatePhi(vector<double> &phi, const mat &u, const mat &u1, double &dt);
void reinitialization(vector<double> &phi);


// Types of Solvers
void solver_LevelSet(mat &u0, vector<double> &phi0, string method);
void solver_LevelSet2Var(mat &u0, mat &u1, vector<double> &phi0, string method);
void solver(mat &u0, string method);


