#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h> 
#include <algorithm> 
#include <vector>

typedef std::vector<double> VectOfDouble;
typedef std::vector<std::vector<double>> VectofVectDouble;

// Routines

void WriteToFile(const VectOfDouble& XCells, const VectofVectDouble& U, const VectofVectDouble& U1, int var);
void ConstructFVDomain(double& dx, VectOfDouble& XCellCentres);
void setInitialConditions(const VectOfDouble& XCellCentres, VectofVectDouble& U, double& Toutput, double& xDiscontinuity);
void ApplyBoundaryConditions(const VectofVectDouble& U, VectofVectDouble& Ubc);
void ComputeDt(const VectofVectDouble& Ubc, const double& dx, const double& time, double& dt);
void UpdatewithFluxes(const VectofVectDouble& UbcOld, VectofVectDouble& Unew, const double& dx, const double& dt);
void HLLFlux(const VectOfDouble& Ql, const VectOfDouble& Qr, VectOfDouble& Fhll);
void HLLCFlux(const VectOfDouble& Ql, const VectOfDouble& Qr, VectOfDouble& Fhllc);
void FORCEFlux(const VectOfDouble& Ql, const VectOfDouble& Qr, VectOfDouble& Fforce, const double& dt, const double& dx);
void WaveEstimates(const VectOfDouble& Ql, const VectOfDouble& Qr, double& Sl, double& Sr);

VectOfDouble EulerFlux(const VectOfDouble& Q);
VectOfDouble conservativeToprimitive(const VectOfDouble& Q);
VectOfDouble primitiveToconservative(const VectOfDouble& W);

double computePressureFromEoS(const double& rho, const double& e);
double computeInternalEnergyFromEoS(const double& rho, double& p);
double computeSoundSpeedFromEoS(const double& rho, const double& p);

void ComputeExactSolution(const VectOfDouble& Wl, const VectOfDouble& Wr, const double& time, const VectOfDouble& X, const double& xd, VectofVectDouble& UExact);
double GuessPressure(const VectOfDouble& Wl, const VectOfDouble& Wr);
double fl_NR(const VectOfDouble W, const double& p);
double fr_NR(const VectOfDouble W, const double& p);
double f_NR(const VectOfDouble Wl, const VectOfDouble Wr, const double& p);
double df_NR(const VectOfDouble Wl, const VectOfDouble Wr, const double& p);
VectOfDouble Sample(const VectOfDouble Wl, const VectOfDouble Wr, const double& Pstar, const double& ustar, const int& Pattern, const double& csi);



