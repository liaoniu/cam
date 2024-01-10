#ifndef GHOST_H
#define GHOST_H

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
void computeDomainBoundaries(mat &u0, mat &u1);

arr indi_primitiveToConserved(const arr &u_pri, const double &Gamma);
arr indi_conservedToPrimitive(const arr &u_con, const double &Gamma);
mat conservedToPrimitive(const mat &u, const double &Gamma);
mat primitiveToConserved(const mat &u, const double &Gamma);
void initializePhi(vector<double> &phi);
double calc_dt(const mat &u_con, const double &Gamma);
double computeTimeStep(mat &u0, mat &u1, const double &GammaL, const double &GammaR);
arr f(const arr &u, const double &Gamma);
double zeta(double delta_minus, double delta_plus, string limiter);
tuple<mat, mat> dataReconstruct(const mat &u, string direction);
tuple<mat, mat> halfTimeUpdate(mat &uL, mat &uR, double dt, const double &Gamma);

void WaveEstimates(const arr& Ql, const arr& Qr, double& Sl, double& Sr, const double &Gamma);
void HLLCFlux(const arr& Ql, const arr& Qr, arr& Fhllc, const double &Gamma);
arr getXFlux(const arr &uL, const arr &uR, double &dt, const double &Gamma, string method); //u_left and u_right
mat split(mat &u, mat &flux, double dt, const double &Gamma, string method);


// LevelSet functions
void initializePhi(vector<double> &phi);
void initializePhi_Helium(vector<double> &phi);
void updatePhi(vector<double> &phi, const mat &u, double &dt);
void updatePhi(vector<double> &phi, const mat &u, const mat &u1, double &dt);
void reinitialization(vector<double> &phi);
void fastSweeping(vector<double> &phi);
int positionOfInterface(vector<double> &phi);
vector<int> posOfAllInterface(vector<double> &phi);

// Ghost Fluid
void reconstructRhoFromConstantEntropy(arr &uI, arr &uG, const double &GammaI, const double &GammaG);
void isobarixFix(mat &u0, mat &u1, int interface, const double &GammaL, const double &GammaR);
void computeGhostFluidBoundaries(mat &u0, mat &u1, vector<double> &phi, const double &GammaL, const double &GammaR);
void computeGhostFluidBoundaries2(mat &u0, mat &u1, vector<double> &phi, const double &GammaL, const double &GammaR);
void computeGhostFluidBoundaries_Helium(mat &u0, mat &u1, vector<double> &phi, const double &GammaL, const double &GammaR);

// Types of Solvers
void solver_LevelSet(mat &u0, vector<double> &phi0, string method);
void solver_LevelSet2Var(mat &u0, mat &u1, vector<double> &phi0, const double &Gamma, string method);
void solver(mat &u0, const double &Gamma, string method);
void solver_Original_Ghost(mat &u0, mat &u1, vector<double> &phi0, const double &GammaL, const double &GammaR, string method);

#endif