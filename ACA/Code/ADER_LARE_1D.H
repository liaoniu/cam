#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h> 
#include <algorithm> 
#include <vector>
#include <functional>

typedef std::vector<double> VectOfDouble;
typedef std::vector<std::vector<double>> VectOfVectDouble;

enum BcType {TRANS=0,PERIOD=1};

// Routines

void WriteToFile(const VectOfDouble& XCells, const VectOfDouble& U, const VectOfDouble& U1);
void ConstructFVDomain(double& dx, VectOfDouble& XCellCentres);
void setInitialConditionsAndComputeExactSolution(const VectOfDouble& XCellCentres, const double& dx, VectOfDouble& U, VectOfDouble& UExact, double& Toutput);
double ShuTestIcs(double x);
double GaussianIcs(double x);
double GaussLegendreNumericalIntegration(const double& xc, const double& dx, std::function<double(double)> func, int nGaussPoints);
void ApplyBoundaryConditions(const VectOfDouble& U, VectOfDouble& Ubc);
void ComputeDt(const VectOfDouble& Ubc, const double& dx, const double& time, double& dt);
void UpdatewithFluxesAndSource(const VectOfDouble& UbcOld, VectOfDouble& Unew, const double& dx, const double& dt);
double ADER_LARE_Flux(const VectOfDouble& Ubc, const int i_interface, const double& dx, const double& dt);


double ADER_LARE_Source(const VectOfDouble& Ubc, const int i_cell, const double& dx, const double& dt); //Task


double Toro_Titarev_GRPSolver(const VectOfDouble& dqL, const VectOfDouble& dqR, const double& dt);
double CauchyKowaProcedure(const VectOfDouble& dq_x, int der);
void   ExactRiemannSolver(const double& ql, const double& qr, double& qGodunovState);

// Factorial

unsigned int factorial(unsigned int n)
{
    if (n == 0)
        return 1;
    return n * factorial(n - 1);
}

// WENO CLASS

#ifndef _WENO_h_
#define _WENO_h_
class WENO1d 
{		 
public:	
	
	WENO1d(std::string PolyDegree, double dx);
	
	VectOfDouble getStencilforInterface(const VectOfDouble& Ubc, const int i_interface, const int NGhost);
	VectOfDouble getStencilforCenter(const VectOfDouble& Ubc, const int i_cell, const int NGhost);
		
	void reconstructAtInterface(const VectOfDouble& qStencil, VectOfDouble& dqL, VectOfDouble& dqR);
	void reconstructWithinCell(const VectOfDouble& qStencil, VectOfVectDouble& dq);
		
	int m_M;
	
	double m_dx;   

protected:
	
	void computeSmoothIndicators(const VectOfDouble& U, VectOfDouble& SmoothInd) const;

	void computeWeights(const VectOfDouble& U, const VectOfDouble& OptimWeights, const VectOfDouble& SmoothInd, VectOfDouble& Weights) const;
  
private:

};
#endif
