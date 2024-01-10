
/* -------------------------------------------------------------------*/
/*                                                                    */
/*                ADER FINITE VOLUME SCHEMES FOR LA(R)E               */
/*                                                                    */
/*   Name of the program : ADER_LARE_1D.cpp                           */
/*                                                                    */
/*   Purpose : Solve numerically the one dimensional linear advection */
/*             (reaction) equation using a finite volume ADER method. */
/*                                                                    */
/*   Output File : ADER_LARE_1D_Output.dat                            */
/*                                                                    */
/*   Programer: Riccardo Demattè                                      */
/*                                                                    */
/*   Last Revision: 20 November 2023                                  */
/*                                                                    */
/*   Theory is found in Chaps. 19,20 of Prof. Toro's book    		      */
/*                                                                    */
/*   This program is part of the practical sessions for               */
/*   Prof. Toro's lectures.                                           */
/*                                                                    */
/* -------------------------------------------------------------------*/

#include "ADER_LARE_1D.H"

using namespace std ;

// Global variables

int N = 400;

int NGhost = 2; // Note: this must be 3 for a third order method

int TestNumber = 2;

double tOut;
double Cfl = 0.95;
double xLeftDomain =  -1.0;
double xRightDomain = 1.0;

double LAMBDA = 1.0;

double BETA = 0.0; // when 0 the reaction sorce term is switched off


BcType BCTYPE = PERIOD; // can be TRANS or PERIOD

int main(void)
{
	VectOfDouble U(N);
	VectOfDouble UBc(N+2*NGhost);
	VectOfDouble UExact(N);

	// Construct the FV domain

	double dx;

	VectOfDouble XCells(N);

	ConstructFVDomain(dx,XCells);

	// Set Initial conditions & compute exact solution

	setInitialConditionsAndComputeExactSolution(XCells,dx,U,UExact,tOut);

	// Time Loop

	int iter = 0;
	double dt;
	double time = 0;

	while(time < tOut)
	{
		// Apply Boundary Conditions

		ApplyBoundaryConditions(U,UBc);

		// Compute Time Step

		ComputeDt(UBc,dx,time,dt);

		// Compute Numerical Fluxes and Update

		UpdatewithFluxesAndSource(UBc,U,dx,dt);

		iter += 1;
		time += dt;

		std::cout<< "Iteration # " << iter << " Time = " << time << std::endl;
	}

	std::cout<<"Program finished " << std::endl;

	// Write Output file

	WriteToFile(XCells,U,UExact);

	return 0 ;
}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                      PrintUtility Routine                          */
/*                                                                    */
/* -------------------------------------------------------------------*/

void WriteToFile(const VectOfDouble& XCells, const VectOfDouble& U, const VectOfDouble& U1)
{
	std::cout<<"Writing output file" << std::endl;

	std::ofstream outFile ;

	outFile.open( "ADER_LARE_1D_Output.dat" ) ;

    outFile << "# This File is the Output of ADER_LARE_1D.cpp" << std::endl;

	outFile << "# " << std::endl ;

	for(int i=0; i<N; i++)
	{
		outFile << XCells[i] << " " << U[i] << " " << " " << U1[i] << std::endl;
	}

	outFile.close();
}


/* -------------------------------------------------------------------*/
/*                                                                    */
/*                  Domain Construction Routine                       */
/*                                                                    */
/* -------------------------------------------------------------------*/

void ConstructFVDomain(double& dx, VectOfDouble& XCellCentres)
{
	dx = (xRightDomain-xLeftDomain)/N;

	for(int i=0 ; i<N ; i++)
	{
		XCellCentres[i] = xLeftDomain + dx*(i+0.5);
	}
}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                  Initial Conditions Routine                        */
/*                                                                    */
/* -------------------------------------------------------------------*/

void setInitialConditionsAndComputeExactSolution(const VectOfDouble& XCellCentres, const double& dx, VectOfDouble& U, VectOfDouble& UExact, double& Toutput)
{
	double PI = 4.0*atan(1.0);

	function<double(double)> ICfunc;

	switch(TestNumber)
	{
		case(1): // SIN(2PIx)
		{
			Toutput = 10;

			ICfunc = [PI](double X) { return sin(PI*X); };

			break;
		}
		case(2): //Shu test
		{
			Toutput = 8.0;

			ICfunc = [PI](double X) { return ShuTestIcs(X); };

			break;
		}
		case(3): //Gaussian
		{
			Toutput = 50.0;

			ICfunc = [PI](double X) { return GaussianIcs(X); };

			break;
		}
		case(4): //SIN4(PIX)
		{
			Toutput = 2.0;

			ICfunc = [PI](double X) { return pow(sin(PI*X),4); };

			break;
		}
		// Add here other tests for LARE...
		default:
		{
			std::cerr << "Test not recognized";
			std::exit(EXIT_FAILURE);
		}
	}

	// Integrate & set IC for the selected test

	for(int i=0; i<N; i++)
	{
		double xcentre = XCellCentres[i];

		U[i] = GaussLegendreNumericalIntegration(xcentre,dx,ICfunc,2);
	}

	// Compute Exact Solution for the selected test

	function<double(double)> EXfunc = [Toutput,ICfunc](double X) { return ICfunc(X)*exp(BETA*Toutput); };

	for(int i=0; i<N; i++)
	{
		double xcentre;

		if(BCTYPE == PERIOD)
		{
			int periods = (int) std::abs(LAMBDA)*Toutput/(xRightDomain-xLeftDomain);

			double Tperiod = (xRightDomain-xLeftDomain)/std::abs(LAMBDA);

			xcentre = XCellCentres[i] - LAMBDA*(Toutput-periods*Tperiod);

		}
		else
		{
			xcentre = XCellCentres[i]-LAMBDA*Toutput;
		}

		UExact[i] = GaussLegendreNumericalIntegration(xcentre,dx,EXfunc,2);
	}
}

double ShuTestIcs(double x)
{
	double IC;

	double teta = 0.005;
	double z = -0.7;
	double gam = log(2.0)/(36.0*teta*teta);
	double a = 0.5;
	double alpha = 10.0;

	function<double(double,double,double)> G = [](double X, double GAM, double Z)   { return  exp(-GAM*pow(X-Z,2)); };
	function<double(double,double,double)> F = [](double X, double ALPHA, double A) { return  sqrt(std::max(1.0-ALPHA*ALPHA*pow(X-A,2),0.0)); };

	if(-0.8 <= x && x <= -0.6)
	{
		IC = 1.0/6.0*(G(x,gam,z-teta)+G(x,gam,z+teta)+4.0*G(x,gam,z));
	}
	else if(-0.4 <= x && x <= -0.2)
	{
		IC = 1.0;
	}
	else if(0.0 <= x && x <= 0.2)
	{
		IC = 1.0-std::abs(10.0*(x-0.1));
	}
	else if(0.4 <= x && x <= 0.6)
	{
		IC = 1.0/6.0*(F(x,alpha,a-teta)+F(x,alpha,a+teta)+4.0*F(x,alpha,a));
	}
	else
	{
		IC = 0.0;
	}

	return IC;
}

double GaussianIcs(double x)
{
	double IC;

	double alpha = 1.0;
	double beta  = 2.0;

	IC = alpha*exp(-beta*x*x);

	return IC;
}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                 Numerical Integration Routines                     */
/*                                                                    */
/* -------------------------------------------------------------------*/

double GaussLegendreNumericalIntegration(const double& xc, const double& dx, function<double(double)> func, int nGaussPoints)
{
	VectOfDouble Weights(nGaussPoints);
	VectOfDouble Location(nGaussPoints);

	switch(nGaussPoints)
	{
		case(1):
		{
			Location = {0.0};
			Weights  = {1.0};

			break;
		}
		case(2):
		{
			Location = {-dx/2.0/sqrt(3.0),dx/2.0/sqrt(3.0)};
			Weights  = {0.5,0.5};
			break;
		}

		// TWO POINTS ARE ENOUGH FOR A THIRD ORDER SCHEME ... FOR HIGHER ORDER SCHEMES (4TH-5TH ETC) USE 3 GAUSS-LEGENDRE POINTS (see Wikipedia)

		// Add here more accurate GL quadratures

		default:
		{
			std::cerr << "Number of GaussLegendre quadrature points not implemented";
			std::exit(EXIT_FAILURE);
			break;
		}
	}

	double integral = 0;

	for(int igauss=0; igauss<nGaussPoints; igauss++)
	{
		integral += Weights[igauss]*func(xc+Location[igauss]);
	}

	return integral;

}


/* -------------------------------------------------------------------*/
/*                                                                    */
/*                 Boundary Conditions Routine                        */
/*                                                                    */
/* -------------------------------------------------------------------*/

void ApplyBoundaryConditions(const VectOfDouble& U, VectOfDouble& Ubc)
{

	// Copy the valid state

	for(int i=0; i<N; i++)
	{
		Ubc[i+NGhost] = U[i];
	}

	// Transmissive or Periodic boundary conditions

	switch(BCTYPE)
	{
		case(TRANS):
		{
			for(int ig=0; ig<NGhost; ig++)
			{
				Ubc[ig] = U[0];
				Ubc[N+NGhost+ig] = U[N-1];
			}
			break;
		}
		case(PERIOD):
		{
			for(int ig=0; ig<NGhost; ig++)
			{
				Ubc[NGhost-1-ig] = U[N-1-ig];
				Ubc[N+NGhost+ig] = U[ig];
			}
			break;
		}
		default:
		{
			std::cerr << "Boundary conditions not recognized";
			std::exit(EXIT_FAILURE);
			break;
		}
	}
}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                  Compute Time Step Routine                         */
/*                                                                    */
/* -------------------------------------------------------------------*/

void ComputeDt(const VectOfDouble& Ubc, const double& dx, const double& time, double& dt)
{
	// Based only on advection velocity
	double Smax= std::abs(LAMBDA);

	dt = std::min(Cfl*dx/Smax,tOut-time);
}


/* -------------------------------------------------------------------*/
/*                                                                    */
/*                       Update Routine                               */
/*                                                                    */
/* -------------------------------------------------------------------*/


void UpdatewithFluxesAndSource(const VectOfDouble& UbcOld, VectOfDouble& Unew, const double& dx, const double& dt)
{

	VectOfDouble NumericalFluxes(N+1);
	VectOfDouble NumericalSources(N,0.0);

	// Compute the numerical Fluxes

	for(int i=0; i<N+1; i++)
	{
		NumericalFluxes[i] = ADER_LARE_Flux(UbcOld,i,dx,dt);
	}

	// Compute numerical Sources

	if(BETA != 0.0)
	{
		for(int i=0; i<N; i++)
		{
			NumericalSources[i] = ADER_LARE_Source(UbcOld,i,dx,dt);
		}
	}

	// Update

	double Fl,Fr,S;

	for(int i=0; i<N; i++)
	{
		Fl = NumericalFluxes[i];
		Fr = NumericalFluxes[i+1];
		S  = NumericalSources[i];

		Unew[i] = UbcOld[i+NGhost] - dt/dx*(Fr-Fl) + dt*S;
	}
}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                  ADER Numerical Flux Routines                      */
/*                                                                    */
/* -------------------------------------------------------------------*/

double ADER_LARE_Flux(const VectOfDouble& Ubc, const int i_interface, const double& dx, const double& dt)
{
	// WENO reconstruction

	VectOfDouble dqWENO_L, dqWENO_R;

	WENO1d Weno("LINEAR",dx); // Note: use parabolas for third order method

	VectOfDouble qStencil = Weno.getStencilforInterface(Ubc,i_interface,NGhost);

	Weno.reconstructAtInterface(qStencil,dqWENO_L,dqWENO_R);

	//! Note: dqWENO_L/R contains reconstructed values and spatial derivatives of q at the left/right of the interface

	// Solve the Generalized Riemann Problem (GRP) : the well-known Toro-Titarev GRP solver is here implemented

	return Toro_Titarev_GRPSolver(dqWENO_L,dqWENO_R,dt);

}

double Toro_Titarev_GRPSolver(const VectOfDouble& dqL, const VectOfDouble& dqR, const double& dt)
{
	// Solve a non-linear RP for the leading term (for LARE equation also this RP is obviousely linear)

	double qLead;

	ExactRiemannSolver(dqL[0],dqR[0],qLead);

	// First order term

	double d1q_x;

	ExactRiemannSolver(dqL[1],dqR[1],d1q_x);

	double d1q_t = CauchyKowaProcedure({qLead,d1q_x},1);

	// Note: Add here second order term for constructing a 3rd order scheme and so on	...

	// Taylor series expansion at the interface

	function<double(double)> qTaylor = [qLead,d1q_t](double tau){ return qLead+d1q_t*tau;}; // Note: Add here second order term for constructing a 3rd order scheme and so on	...

	// Flux integration

	function<double(double)> Fluxfuction = [qTaylor](double tau){ return LAMBDA*qTaylor(tau); };

	double TTFlux = GaussLegendreNumericalIntegration(dt/2.0,dt,Fluxfuction,1); // 1 point (mid rule) is sufficient for 2nd order ADER. Use 2 points for 3rd order

	return TTFlux;

}

double CauchyKowaProcedure(const VectOfDouble& dq_x, int der)
{
	return pow(-1,der)*pow(LAMBDA,der)*dq_x[der];

	// Note: modify the CK procedure accordingly when the reaction term is present (BETA != 0)

}

void ExactRiemannSolver(const double& ql, const double& qr, double& qGodunovState)
{
	qGodunovState = (LAMBDA >=0) ? ql : qr;
}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                 ADER Numerical Source Routines                     */
/*                                                                    */
/* -------------------------------------------------------------------*/


double ADER_LARE_Source(const VectOfDouble& Ubc, const int i_cell, const double& dx, const double& dt)
{
	double Source = 0 ;

	// WENO reconstruction

	VectOfVectDouble dq;

	WENO1d Weno("LINEAR",dx); // Note: use parabolas for third order method

	VectOfDouble qStencil = Weno.getStencilforCenter(Ubc,i_cell,NGhost);

	Weno.reconstructWithinCell(qStencil,dq);

	// TO BE COMPLETED...

	// OUTLINE OF THE MAIN PASSAGES

	// construct taylor expansion at the gauss quadrature point within the cell

	// CK procedure needed to convert from spatial derivatives to time derivatives

	// the computation of the spatial derivatives doesn't require the solution of any Riemann problem since the data is smooth within the cell (no jumps). The value returned by WENO can be thus used directly.

	// Numerically integrate the source term in both space and time:
	//                                                              for a second order method this is trival (mid-point rule). use 1 point located in the middle of the cell and at tn+dt/2
	//                                                              for a third order method this requires the construction of 2 Taylor expansions within the cell (at xGauss1 and xGauss2)
	//                                                              each of which is then evaluated at tGauss1 and tGauss2.


	return Source;

}


/* -------------------------------------------------------------------*/
/*                                                                    */
/*                  WENO Reconstruction Routines                      */
/*                                                                    */
/* -------------------------------------------------------------------*/


// Note: the WENO reconstruction is implemented for linear and parabolic polynomials which can be used to construct 2nd and 3rd order ADER methods, respectively.
//       Increase the order of the reconstruction if you want to build higher order schemes! Not a trivial task :)

WENO1d::WENO1d(std::string PolyDegree, double dx)
{
	m_dx = dx;

	if(PolyDegree == "LINEAR")
	{
		m_M = 1;
	}
	else if(PolyDegree == "PARABOLIC")
	{
		m_M = 2;
	}
	else
	{
		std::cerr << "WENO polynomial not implemented";
		std::exit(EXIT_FAILURE);
	}
}

VectOfDouble WENO1d::getStencilforInterface(const VectOfDouble& Ubc, const int i_interface, const int NGhost)
{
	int numberofGhostRequired = m_M+1;

	if(NGhost < numberofGhostRequired)
	{
		std::cerr << "The WENO reconstruction requires at least " << numberofGhostRequired << " ghost cells. Change the number of ghost cells used.";
		std::exit(EXIT_FAILURE);
	}

	VectOfDouble Stencil(2*m_M+2);

	for(int i=0; i<Stencil.size(); i++)
	{
		Stencil[i] = Ubc[i_interface+NGhost-numberofGhostRequired+i];
	}

	return Stencil;
}

VectOfDouble WENO1d::getStencilforCenter(const VectOfDouble& Ubc, const int i_cell, const int NGhost)
{
	int numberofGhostRequired = m_M;

	if(NGhost < numberofGhostRequired)
	{
		std::cerr << "The WENO reconstruction requires at least " << numberofGhostRequired << " ghost cells. Change the number of ghost cells used.";
		std::exit(EXIT_FAILURE);
	}

	VectOfDouble Stencil(2*m_M+1);

	for(int i=0; i<Stencil.size(); i++)
	{
		Stencil[i] = Ubc[i_cell+NGhost-numberofGhostRequired+i];
	}

	return Stencil;
}

void WENO1d::reconstructAtInterface(const VectOfDouble& qStencil, VectOfDouble& dqL, VectOfDouble& dqR)
{
	VectOfDouble U1l = {qStencil.begin(),qStencil.end()-1};
	VectOfDouble U1r = {qStencil.begin()+1,qStencil.end()};

	VectOfDouble Weights(m_M+1);
	VectOfDouble SmoothInd(m_M+1);
	VectOfDouble OptimWeights(m_M+1);

	dqL.resize(m_M+1);
	dqR.resize(m_M+1);

    switch(m_M)
    {
		case 1:
		{
			double Up0,Up1;

			// First compute the recons value at the left of the interf.

			OptimWeights[0] = 1.0/3.0;
			OptimWeights[1] = 2.0/3.0;

			computeSmoothIndicators(U1l,SmoothInd);

			computeWeights(U1l,OptimWeights,SmoothInd,Weights);

			//! Reconstructed value

			Up0 = (-1.0*U1l[0] + 3.0*U1l[1])/2.0;
			Up1 = ( 1.0*U1l[1] + 1.0*U1l[2])/2.0;

			dqL[0] = Weights[0]*Up0+Weights[1]*Up1;

			//! First derivative

			Up0 = (-U1l[0] + U1l[1])/m_dx;
			Up1 = (-U1l[1] + U1l[2])/m_dx;

			dqL[1] = Weights[0]*Up0+Weights[1]*Up1;

			// Then compute the recons value at the right of the interf.

			OptimWeights[0] = 2.0/3.0;
			OptimWeights[1] = 1.0/3.0;

			computeSmoothIndicators(U1r,SmoothInd);

			computeWeights(U1r,OptimWeights,SmoothInd,Weights);

			//! Reconstructed value

			Up0 = (1.0*U1r[0] + 1.0*U1r[1])/2.0;
			Up1 = (3.0*U1r[1] - 1.0*U1r[2])/2.0;

			dqR[0] = Weights[0]*Up0+Weights[1]*Up1;

			//! First derivative

			Up0 = (-U1r[0] + U1r[1])/m_dx;
			Up1 = (-U1r[1] + U1r[2])/m_dx;

			dqR[1] = Weights[0]*Up0+Weights[1]*Up1;

			break;
		}
		case 2:
		{
			double Up0,Up1,Up2 ;

			// First compute the recons value at the left of the interf.

			OptimWeights[0] = 1.0/10.0;
			OptimWeights[1] = 3.0/5.0 ;
			OptimWeights[2] = 3.0/10.0;

			computeSmoothIndicators(U1l,SmoothInd);

			computeWeights(U1l,OptimWeights,SmoothInd,Weights);

			//! Reconstructed value

			Up0 = ( 2.0*U1l[0] - 7.0*U1l[1] + 11.0*U1l[2])/6.0 ;
			Up1 = (-1.0*U1l[1] + 5.0*U1l[2] +  2.0*U1l[3])/6.0 ;
			Up2 = ( 2.0*U1l[2] + 5.0*U1l[3] -  1.0*U1l[4])/6.0 ;

			dqL[0] = Weights[0]*Up0+Weights[1]*Up1+Weights[2]*Up2;

			//! First derivative

			Up0 = ( 1.0*U1l[0] - 3.0*U1l[1] + 2.0*U1l[2])/m_dx;
			Up1 = (-0.0*U1l[1] - 1.0*U1l[2] + 1.0*U1l[3])/m_dx;
			Up2 = (-1.0*U1l[2] + 1.0*U1l[3] + 0.0*U1l[4])/m_dx;

			dqL[1] = Weights[0]*Up0+Weights[1]*Up1+Weights[2]*Up2;

			//! Second derivative

			Up0 = (U1l[0] - 2.0*U1l[1] + U1l[2])/pow(m_dx,2);
			Up1 = (U1l[1] - 2.0*U1l[2] + U1l[3])/pow(m_dx,2);
			Up2 = (U1l[2] - 2.0*U1l[3] + U1l[4])/pow(m_dx,2);

			dqL[2] = Weights[0]*Up0+Weights[1]*Up1+Weights[2]*Up2;

			// Then compute the recons value at the right of the interf.

			OptimWeights[0] = 3.0/10.0;
			OptimWeights[1] = 3.0/5.0 ;
			OptimWeights[2] = 1.0/10.0;

			computeSmoothIndicators(U1r,SmoothInd);

			computeWeights(U1r,OptimWeights,SmoothInd,Weights);

			//! Reconstructed value

			Up0 = (-1.0*U1r[0] + 5.0*U1r[1] + 2.0*U1r[2])/6.0 ;
			Up1 = ( 2.0*U1r[1] + 5.0*U1r[2] - 1.0*U1r[3])/6.0 ;
			Up2 = (11.0*U1r[2] - 7.0*U1r[3] + 2.0*U1r[4])/6.0 ;

			dqR[0] = Weights[0]*Up0+Weights[1]*Up1+Weights[2]*Up2;

			//! First derivative

			Up0 = (-0.0*U1r[0] - 1.0*U1r[1] + 1.0*U1r[2])/m_dx;
			Up1 = (-1.0*U1r[1] + 1.0*U1r[2] - 0.0*U1r[3])/m_dx;
			Up2 = (-2.0*U1r[2] + 3.0*U1r[3] - 1.0*U1r[4])/m_dx;

			dqR[1] = Weights[0]*Up0+Weights[1]*Up1+Weights[2]*Up2;

			//! Second derivative

			Up0 = (U1r[0] - 2.0*U1r[1] + U1r[2])/pow(m_dx,2);
			Up1 = (U1r[1] - 2.0*U1r[2] + U1r[3])/pow(m_dx,2);
			Up2 = (U1r[2] - 2.0*U1r[3] + U1r[4])/pow(m_dx,2);

			dqR[2] = Weights[0]*Up0+Weights[1]*Up1+Weights[2]*Up2;

			break;
		}
		// add here higher order reconstructions
		default:
		{
			std::cerr << "WENO polynomial not implemented";
			std::exit(EXIT_FAILURE);
			break;
		}
	}
}

void WENO1d::reconstructWithinCell(const VectOfDouble& qStencil, VectOfVectDouble& dq)
{
	//! For Linear polynomial the weno reconstruction is evaluated at the cell centre
	//! For Parabolic polynomial the weno reconstruction is evaluated at 2 Gauss Legendre points
	//! Also for cubic polynomials 2 Gauss Legendre points are sufficient

	VectOfDouble U = qStencil;

	VectOfDouble Delta(m_M+1);
	VectOfDouble Weights(m_M+1);
	VectOfDouble SmoothInd(m_M+1);
	VectOfDouble OptimWeights(m_M+1);

	int npoints;

	if(m_M == 1)
	{
		npoints = 1;
	}
	else if( m_M == 2 )
	{
		npoints = 2;
	}

	dq.resize(m_M+1,VectOfDouble(npoints));

	VectOfDouble dqR(m_M+1);
	VectOfDouble dqL(m_M+1);

	double sqrt3 = sqrt(3.0);

    switch(m_M)
    {
		case 1:
		{
			double Up0,Up1;

			// Smoothness Indicators

			computeSmoothIndicators(U,SmoothInd);

			OptimWeights[0] = 0.5;
			OptimWeights[1] = 0.5;

			computeWeights(U,OptimWeights,SmoothInd,Weights);

			Up0 = (1.0*U[0] + 1.0*U[1])/2.0;
			Up1 = (3.0*U[1] - 1.0*U[2])/2.0;

			dq[0][0] = Weights[0]*Up0+Weights[1]*Up1;

			// First derivative

			Up0 = (-U[0] + U[1])/m_dx;
			Up1 = (-U[1] + U[2])/m_dx;

			dq[1][0] = Weights[0]*Up0+Weights[1]*Up1;

			dq[0][0] += dq[1][0]*m_dx/2.0;

			break;

		}
		case 2:
		{
			double Up0,Up1,Up2 ;

			// Smoothness Indicators

			computeSmoothIndicators(U,SmoothInd);

			// Deltas

			Delta[0] = ( U[0] - 4.0*U[1] + 3.0*U[2])*sqrt3/12.0;
			Delta[1] = (-U[1] + U[3])*sqrt3/12.0;
			Delta[2] = (-3.0*U[2] + 4.0*U[3] - U[4])*sqrt3/12.0;

			// First compute the recons value at the lower Gauss point G = -h/2/SQRT(3)

			OptimWeights[0] = sqrt3*(70.0*sqrt3+1.0)/1080.0;
			OptimWeights[1] = 11.0/18.0;
			OptimWeights[2] = sqrt3*(70.0*sqrt3-1.0)/1080.0;

			computeWeights(U,OptimWeights,SmoothInd,Weights);

			// Reconstructed value

			Up0 = U[2] - Delta[0];
			Up1 = U[2] - Delta[1];
			Up2 = U[2] - Delta[2];

			dq[0][0] = Weights[0]*Up0+Weights[1]*Up1+Weights[2]*Up2;

			// First derivative

			Up0 =  (9.0*U[2]-12.0*U[1]+3.0*U[0]-sqrt3*U[2]+2.0*sqrt3*U[1]-sqrt3*U[0])/(6.0*m_dx);
			Up1 = -(3.0*U[1]-3.0*U[3]-2.0*sqrt3*U[2]+sqrt3*U[1]+sqrt3*U[3])/(6.0*m_dx);
			Up2 = -(9.0*U[2]-12.0*U[3]+3.0*U[4]+sqrt3*U[2]-2.0*sqrt3*U[3]+sqrt3*U[4])/(6.0*m_dx);

			dq[1][0] = Weights[0]*Up0+Weights[1]*Up1+Weights[2]*Up2;

			// Second derivative

			Up0 = (U[0] - 2.0*U[1] + U[2])/pow(m_dx,2);
			Up1 = (U[1] - 2.0*U[2] + U[3])/pow(m_dx,2);
			Up2 = (U[2] - 2.0*U[3] + U[4])/pow(m_dx,2);

			dq[2][0] = Weights[0]*Up0+Weights[1]*Up1+Weights[2]*Up2;

			// Then compute the recons value at the upper Gauss point G = +h/2/SQRT(3)

			OptimWeights[0] = sqrt3*(70.0*sqrt3-1.0)/1080.0;
			OptimWeights[1] = 11.0/18.0;
			OptimWeights[2] = sqrt3*(70.0*sqrt3+1.0)/1080.0;

			computeWeights(U,OptimWeights,SmoothInd,Weights);

			// Reconstructed value

			Up0 = U[2] + Delta[0];
			Up1 = U[2] + Delta[1];
			Up2 = U[2] + Delta[2];

			dq[0][1] = Weights[0]*Up0+Weights[1]*Up1+Weights[2]*Up2;

			// First derivative

			Up0 = -(-9.0*U[2]+12.0*U[1]-3.0*U[0]-sqrt3*U[2]+2.0*sqrt3*U[1]-sqrt3*U[0])/(6.0*m_dx);
			Up1 =  (-3.0*U[1]+3.0*U[3]-2.0*sqrt3*U[2]+sqrt3*U[1]+sqrt3*U[3])/(6.0*m_dx);
			Up2 =  (-9.0*U[2]+12.0*U[3]-3.0*U[4]+sqrt3*U[2]-2.0*sqrt3*U[3]+sqrt3*U[4])/(6.0*m_dx);

			dq[1][1] = Weights[0]*Up0+Weights[1]*Up1+Weights[2]*Up2;

			// Second derivative

			Up0 = (U[0] - 2.0*U[1] + U[2])/pow(m_dx,2);
			Up1 = (U[1] - 2.0*U[2] + U[3])/pow(m_dx,2);
			Up2 = (U[2] - 2.0*U[3] + U[4])/pow(m_dx,2);

			dq[2][1] = Weights[0]*Up0+Weights[1]*Up1+Weights[2]*Up2;

			break;
		}
		// add here higher order reconstructions
		default:
		{
			std::cerr << "WENO polynomial not implemented";
			std::exit(EXIT_FAILURE);
			break;
		}
	}
}

void WENO1d::computeSmoothIndicators(const VectOfDouble& U, VectOfDouble& SmoothInd) const
{
	//! U has size 2*M+1
	//! SmoothInd has size M+1

	switch(m_M)
    {
		case 1:
		{
			SmoothInd[0] = pow(-U[0]+U[1],2);
			SmoothInd[1] = pow(-U[1]+U[2],2);

			break;
		}

		case 2:
		{
			SmoothInd[0] = 13.0/12.0*pow(U[0] - 2.0*U[1] + U[2],2) + 1.0/4.0*pow(U[0] - 4.0*U[1] + 3.0*U[2],2);
			SmoothInd[1] = 13.0/12.0*pow(U[1] - 2.0*U[2] + U[3],2) + 1.0/4.0*pow(U[1] - U[3],2);
			SmoothInd[2] = 13.0/12.0*pow(U[2] - 2.0*U[3] + U[4],2) + 1.0/4.0*pow(3.0*U[2] - 4.0*U[3] + U[4],2);

			break;
		}
		// add here higher order reconstructions
		default:
		{
			std::cerr << "WENO polynomial not implemented";
			std::exit(EXIT_FAILURE);
		}
	}
}

void WENO1d::computeWeights(const VectOfDouble& U, const VectOfDouble& OptimWeights, const VectOfDouble& SmoothInd, VectOfDouble& Weights) const
{
	//! U has size 2*M+1
	//! SmoothInd has size M+1
	//! OptimWeights has size M+1

	double Eps = 1.0E-6;
	double WeightsSum ;

	switch(m_M)
    {
		case 1:
		{
			Weights[0] = OptimWeights[0]/pow(Eps+SmoothInd[0],2);
			Weights[1] = OptimWeights[1]/pow(Eps+SmoothInd[1],2);

			WeightsSum = Weights[0]+Weights[1];

			Weights[0] /= WeightsSum;
			Weights[1] /= WeightsSum;

			break;
		}

		case 2:
		{
			// Weno non-linear weights

			Weights[0] = OptimWeights[0]/pow(Eps+SmoothInd[0],2);
			Weights[1] = OptimWeights[1]/pow(Eps+SmoothInd[1],2);
			Weights[2] = OptimWeights[2]/pow(Eps+SmoothInd[2],2);

			WeightsSum = Weights[0]+Weights[1]+Weights[2];

			Weights[0] /= WeightsSum;
			Weights[1] /= WeightsSum;
			Weights[2] /= WeightsSum;

			break;
		}
		// add here higher order reconstructions
		default:
		{
			std::cerr << "WENO polynomial not implemented";
			std::exit(EXIT_FAILURE);
		}
	}
}




/* ------------------------------------------------------------------ */
/*                                                                    */
/*                       END OF THE FILE                              */
/*                                                                    */
/* ------------------------------------------------------------------ */
