#ifndef GHOST_H
#define GHOST_H

#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;
using mat = vector<array<double, 5>>;   // u[i][var] = var[i]
using arr = array<double, 5>;



bool inCircle(double x, double y, double R);
void transmissiveBC(mat& u);
void transmissiveBC(vector<double> &phi);
void checkAlpha(mat &u);
void checkAlpha(arr &u);
void computeDomainBoundaries(mat &u0, mat &u1);

arr indi_primitiveToConserved(const arr &u_pri);
arr indi_conservedToPrimitive(const arr &u_con);
mat conservedToPrimitive(const mat &u);
mat primitiveToConserved(const mat &u);
void initializePhi(vector<double> &phi);
double computeSoundSpeed_Allaire(const arr &u_con);
tuple<double, double> computeGammaPinfty_Alllaire(double a1);
double calc_dt(const mat &u_con, const double &Gamma);
double calc_dt_Allaire(const mat &u_con);
arr f(const arr &u);


void WaveEstimates(const arr& Ql, const arr& Qr, double& Sl, double& Sr);
void HLLCFlux(const arr& Ql, const arr& Qr, arr& Fhllc, double &S_star);
arr getXFlux(const arr &uL, const arr &uR, double &dt, string method, double &S);
mat split(mat &u, mat &flux, double dt, string method);



// Types of Solvers
void solver_LevelSet(mat &u0, vector<double> &phi0, string method);
void solver_LevelSet2Var(mat &u0, mat &u1, vector<double> &phi0, const double &Gamma, string method);
void solver(mat &u0, string method);
void solver_Original_Ghost(mat &u0, mat &u1, vector<double> &phi0, const double &GammaL, const double &GammaR, string method);

#endif