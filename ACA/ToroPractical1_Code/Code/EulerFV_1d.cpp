
/* -------------------------------------------------------------------*/
/*                                                                    */
/*               FINITE VOLUME SCHEMES FOR THE EULER EQs              */
/*                                                                    */
/*   Name of the program : EulerFV_1d.C                               */
/*                                                                    */
/*   Purpose : Solve numerically the one dimensional Euler equations  */
/*             using a finite volume method. Centred and upwind       */
/*             numerical fluxes are implemented and explored.         */
/*                                                                    */
/*   Output File : EulerFV_1dOutput.dat                               */
/*                                                                    */
/*   Programer: Riccardo Demattè                                      */
/*                                                                    */
/*   Last Revision: 9 November 2023                                   */
/*                                                                    */
/*   Theory is found in Chaps. 3,4,5 and 6 of Prof. Toro's book       */
/*                                                                    */
/*   This program is part of the practical sessions for               */
/*   Prof. Toro's lectures.                                           */
/*                                                                    */
/* -------------------------------------------------------------------*/

#include "EulerFV_1d.H"

using namespace std ;

// Global variables

int TestNumber = 1;
int N = 200;
int NGhost = 1;
int nVar = 3;

double tOut,RP_Xdisc;
double Cfl = 0.8;
double xLeftDomain = 0;
double xRightDomain = 1.0;
double Gamma = 1.4;

VectOfDouble RP_LeftState;
VectOfDouble RP_RightState;

int main(void)
{
	VectofVectDouble U(N,VectOfDouble(nVar));
	VectofVectDouble UBc(N+2*NGhost,VectOfDouble(nVar));
	VectofVectDouble UExact(N,VectOfDouble(nVar));

	// Construct the FV domain

	double dx;

	VectOfDouble XCells(N);

	ConstructFVDomain(dx,XCells);

	// Set Initial conditions

	setInitialConditions(XCells,U,tOut,RP_Xdisc);

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

		UpdatewithFluxes(UBc,U,dx,dt);

		iter += 1;
		time += dt;

		std::cout<< "Iteration # " << iter << " Time = " << time << std::endl;
	}

	// Compute Exact Solution

	ComputeExactSolution(RP_LeftState,RP_RightState,time,XCells,RP_Xdisc,UExact);

	// Write Output file

	WriteToFile(XCells,U,UExact,0);

	return 0 ;
}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                      PrintUtility Routine                          */
/*                                                                    */
/* -------------------------------------------------------------------*/

void WriteToFile(const VectOfDouble& XCells, const VectofVectDouble& U, const VectofVectDouble& U1, int var)
{
	std::ofstream outFile ;

	outFile.open( "EulerFV_1dOutput.dat" ) ;

    outFile << "# This File is the Output of Schemes_FV.C" << std::endl;

	outFile << "# " << std::endl ;

	for(int i=0; i<N; i++)
	{
		outFile << XCells[i] << " " << U[i][var] << " " << " " << U1[i][var] << std::endl;
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

void setInitialConditions(const VectOfDouble& XCellCentres, VectofVectDouble& U, double& Toutput, double& xDiscontinuity)
{

	// Toro's tests 1,2,3

	VectOfDouble Wl(nVar),Wr(nVar),Ql(nVar),Qr(nVar);

	switch(TestNumber)
	{
		case(1):
		{
			Toutput = 0.2;
			xDiscontinuity = 0.3;

			Wl = {1.0,0.75,1.0};
			Wr = {0.125,0.0,0.1};

			break;
		}

		case(2):
		{
			Toutput = 0.15;
			xDiscontinuity = 0.5;

			Wl = {1.0,-2.0,0.4};
			Wr = {1.0,+2.0,0.4};

			break;
		}

		case(3):
		{
			Toutput = 0.012;
			xDiscontinuity = 0.5;

			Wl = {1.0,0.0,1000.0};
			Wr = {1.0,0.0,0.01};

			break;
		}

		// Add here other Riemann tests .. 4,5,6 from Toro's book

	}

	Ql = primitiveToconservative(Wl);
	Qr = primitiveToconservative(Wr);

	for(int i=0; i<N; i++)
	{
		if(XCellCentres[i] <= xDiscontinuity)
		{
			U[i] = Ql;
		}
		else
		{
			U[i] = Qr;
		}
	}

	// Save Left and Right states of the RP problem

	RP_LeftState  = Wl;
	RP_RightState = Wr;

}


/* -------------------------------------------------------------------*/
/*                                                                    */
/*                 Boundary Conditions Routine                        */
/*                                                                    */
/* -------------------------------------------------------------------*/

void ApplyBoundaryConditions(const VectofVectDouble& U, VectofVectDouble& Ubc)
{

	// Copy the valid state

	for(int i=0; i<N; i++)
	{
		Ubc[i+NGhost] = U[i];
	}

	// Transmissive boundary conditions

	for(int ig=0; ig<NGhost; ig++)
	{
		Ubc[ig] = U[0];
		Ubc[N+NGhost+ig] = U[N-1];
	}

	// Implement here other boundary conditions ... (reflective for example)

}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                  Compute Time Step Routine                         */
/*                                                                    */
/* -------------------------------------------------------------------*/

void ComputeDt(const VectofVectDouble& Ubc, const double& dx, const double& time, double& dt)
{
	double Sl,Sr,Smax=0;
	VectOfDouble Ql(nVar),Qr(nVar);

	for(int i=0; i<N+1; i++)
	{
		Ql = Ubc[i+NGhost-1];
		Qr = Ubc[i+NGhost];

		WaveEstimates(Ql,Qr,Sl,Sr);

		Smax = std::max(Smax,std::max(std::abs(Sl),std::abs(Sr)));
	}

	dt = std::min(Cfl*dx/Smax,tOut-time);
}


/* -------------------------------------------------------------------*/
/*                                                                    */
/*                       Update Routine                               */
/*                                                                    */
/* -------------------------------------------------------------------*/


void UpdatewithFluxes(const VectofVectDouble& UbcOld, VectofVectDouble& Unew, const double& dx, const double& dt)
{

	VectofVectDouble NumericalFluxes(N+1,VectOfDouble(nVar));

	// Compute the numerical FLux

	VectOfDouble Ql(nVar),Qr(nVar);

	for(int i=0; i<N+1; i++)
	{
		Ql = UbcOld[i+NGhost-1];
		Qr = UbcOld[i+NGhost];

		HLLCFlux(Ql,Qr,NumericalFluxes[i]);  // Implement and Replace here other numerical fluxes... for example HLLC,FORCE or TV
	}

	// Update

	VectOfDouble Fl(nVar),Fr(nVar);

	for(int i=0; i<N; i++)
	{
		Fl = NumericalFluxes[i];
		Fr = NumericalFluxes[i+1];

		for(int n=0; n<nVar; n++)
		{
			Unew[i][n] = UbcOld[i+NGhost][n] + dt/dx*(Fl[n]-Fr[n]);
		}
	}
}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                    Numerical Flux Routines                         */
/*                                                                    */
/* -------------------------------------------------------------------*/

void HLLFlux(const VectOfDouble& Ql, const VectOfDouble& Qr, VectOfDouble& Fhll)
{
	double Sl,Sr;

	WaveEstimates(Ql,Qr,Sl,Sr);

	VectOfDouble Fl(nVar),Fr(nVar);

	Fl = EulerFlux(Ql);
	Fr = EulerFlux(Qr);

	if(Sl >= 0)
	{
		Fhll = Fl;
		return;
	}

	if(Sr <= 0)
	{	Fhll = Fr;
		return;
	}

	for(int n=0; n<nVar; n++)
	{
		Fhll[n] = (Sr*Fl[n]-Sl*Fr[n]+Sl*Sr*(Qr[n]-Ql[n]))/(Sr-Sl);
	}
}


void HLLCFlux(const VectOfDouble& Ql, const VectOfDouble& Qr, VectOfDouble& Fhllc)
{
	double Sl, Sr, S_star;

	WaveEstimates(Ql, Qr, Sl, Sr);

	VectOfDouble Fl(nVar), Fr(nVar), Fl_star(nVar), Fr_star(nVar);

	Fl = EulerFlux(Ql);
	Fr = EulerFlux(Qr);
	
	if(Sl >= 0)
	{
		Fhllc = Fl;
		return;
	}

	if(Sr <= 0)
	{
		Fhllc = Fr;
		return;
	}

	
	double pl, pr, rhol, rhor, ul, ur;
	VectOfDouble Wl(nVar), Wr(nVar);

	Wl = conservativeToprimitive(Ql);
	Wr = conservativeToprimitive(Qr);
	tie(rhol, ul, pl) = tie(Wl[0], Wl[1], Wl[2]);
	tie(rhor, ur, pr) = tie(Wr[0], Wr[1], Wr[2]);
	S_star = (pr - pl + rhol*ul*(Sl - ul) - rhor*ur*(Sr - ur))/(rhol*(Sl - ul) - rhor*(Sr - ur));
	

	VectOfDouble D_star = {0, 1, S_star}, Ql_star(nVar), Qr_star(nVar);			// Change if move to higher dimensions
	double pstarl, pstarr;
	pstarl = pl + rhol*(Sl - ul)*(S_star - ul);
	pstarr = pr + rhor*(Sr - ur)*(S_star - ur);

	for(int i = 0; i != nVar; i++)
	{
		Ql_star[i] = (Sl*Ql[i] - Fl[i] + pstarl*D_star[i])/(Sl - S_star);
		Qr_star[i] = (Sr*Qr[i] - Fr[i] + pstarr*D_star[i])/(Sr - S_star);
	}

	for(int i = 0; i != nVar; i++)
	{
		Fl_star[i] = Fl[i] + Sl*(Ql_star[i] - Ql[i]);
		Fr_star[i] = Fr[i] + Sr*(Qr_star[i] - Qr[i]);
	}


	if(S_star >= 0)
	{
		Fhllc = Fl_star;
		return;
	}

	else
	{
		Fhllc = Fr_star;
		return;
	}
	
}

void FORCEFlux(const VectOfDouble& Ql, const VectOfDouble& Qr, VectOfDouble& Fforce, const double& dt, const double& dx)
{
	VectOfDouble Qhalf(nVar), Fl(nVar), Fr(nVar);

	Fl = EulerFlux(Ql);
	Fr = EulerFlux(Qr);
	
	for(int i = 0; i != nVar; i++)
	{
		Qhalf[i] = 0.5*(Ql[i] + Qr[i]) + dt/2/dx*(Fl[i] - Fr[i]);
	}

	VectOfDouble Fhalf(nVar);

	Fhalf = EulerFlux(Qhalf);
	

	for(int i = 0; i != nVar; i++)
	{
		Fforce[i] = 0.5*(Fhalf[i] + 0.5*(Fl[i] + Fr[i])) + dx/dt/4*(Ql[i] - Qr[i]);
	}

}

//void TVFlux(const VectOfDouble& Ql, const VectOfDouble& Qr, VectOfDouble& Fhllc){}

void WaveEstimates(const VectOfDouble& Ql, const VectOfDouble& Qr, double& Sl, double& Sr)
{
	VectOfDouble Wl(nVar),Wr(nVar);

	Wl = conservativeToprimitive(Ql);
	Wr = conservativeToprimitive(Qr);

	double rhol = Wl[0];
	double rhor = Wr[0];

	double ul = Wl[1];
	double ur = Wr[1];

	double pl = Wl[2];
	double pr = Wr[2];

	double al = computeSoundSpeedFromEoS(rhol,pl);
	double ar = computeSoundSpeedFromEoS(rhor,pr);

	// Pressure–Based Wave Speed Estimates (ideal gases)

	double ql,qr;

	// Two-rarefaction Riemann solver TRRS for computing Pstar

	double z = (Gamma-1)/(2.0*Gamma);

	double pLR = pow(pl/pr,z);

	double ustar = (pLR*ul/al+ur/ar+2.0*(pLR-1.0)/(Gamma-1.0))/(pLR/al+1.0/ar);

	double pstar = 0.5*(pl*pow(1.0+(Gamma-1.0)/(2.0*al)*(ul-ustar),1.0/z)+pr*pow(1.0+(Gamma-1.0)/(2.0*ar)*(ustar-ur),1.0/z));

	if(pstar <= pl)
	{
		ql = 1.0;
	}
	else
	{
		ql = sqrt(1.0+(Gamma+1.0)/(2.0*Gamma)*(pstar/pl-1.0));
	}

	if(pstar <= pr)
	{
		qr = 1.0;
	}
	else
	{
		qr = sqrt(1.0+(Gamma+1.0)/(2.0*Gamma)*(pstar/pr-1.0));
	}

	Sl = ul-al*ql;
	Sr = ur+ar*qr;

	// Add here different Wave speed estimates ....


}

VectOfDouble EulerFlux(const VectOfDouble& Q)
{
	VectOfDouble F(nVar);

	VectOfDouble W = conservativeToprimitive(Q);

	double rho = W[0];
	double u   = W[1];
	double p   = W[2];
	double E   = Q[2];

	F[0] = rho*u;
	F[1] = rho*u*u+p;
	F[2] = u*(E+p);

	return F;

}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*          Primitive and Conservative Variables Routines             */
/*                                                                    */
/* -------------------------------------------------------------------*/

VectOfDouble conservativeToprimitive(const VectOfDouble& Q)
{
	VectOfDouble W(nVar);

	double rho = Q[0];
	double u   = Q[1]/Q[0];
	double E   = Q[2];
  double kin = 0.5*rho*u*u;
	double e   = (E-kin)/rho;

	double p = computePressureFromEoS(rho,e);

	W[0] = rho;
	W[1] = u;
	W[2] = p;

	return W;

}

VectOfDouble primitiveToconservative(const VectOfDouble& W)
{
	VectOfDouble Q(nVar);

	double rho = W[0];
	double u   = W[1];
	double p   = W[2];
    double kin = 0.5*rho*u*u;
	double e   = computeInternalEnergyFromEoS(rho,p);

	double E = rho*e+kin;

	Q[0] = rho;
	Q[1] = rho*u;
	Q[2] = E;

	return Q;

}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                    Equation of State Routines                      */
/*                                                                    */
/* -------------------------------------------------------------------*/

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

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                   RP: Exact Solution Routines                      */
/*                                                                    */
/* -------------------------------------------------------------------*/

// This implementation is specific for ideal gases

void ComputeExactSolution(const VectOfDouble& Wl, const VectOfDouble& Wr, const double& time, const VectOfDouble& X, const double& xd, VectofVectDouble& UExact)
{
	double G1 = (Gamma-1.0)/(2.0*Gamma);
	double G4 = 2.0/(Gamma-1.0);

	double rhol = Wl[0];
	double rhor = Wr[0];

	double ul = Wl[1];
	double ur = Wr[1];

	double pl = Wl[2];
	double pr = Wr[2];


	double al = sqrt(Gamma*pl/rhol);
	double ar = sqrt(Gamma*pr/rhor);

    // Check that Vacuum is not generated

    double DeltaU = ur-ul;

    bool Vacuum = G4*(al+ar) <= DeltaU;

    if(Vacuum)
    {
		std::cerr << "The initial data is such that vacuum is generated. Program stopped." ;
		std::_Exit(EXIT_FAILURE) ;
	}

    // Newton Raphson method

    int iter = 0;
    int iterMax = 100;

    double Err = 1.0;
    double Tol = 1.0E-9;

	double Pold = GuessPressure(Wl,Wr);

	double Piter;

	while(Err > Tol && iter < iterMax)
	{
		Piter = Pold-(f_NR(Wl,Wr,Pold)/df_NR(Wl,Wr,Pold));

		Err = 2.0*std::abs((Piter-Pold)/(Piter+Pold));

		Pold = std::max(0.0,Piter) ;

		iter += 1 ;
	}

	double Pstar = Pold;
	double ustar = 0.5*(ul+ur)+0.5*(fr_NR(Wr,Pstar)-fl_NR(Wl,Pstar));

	// Wave pattern

	int Pattern;

	if(Pstar < pl && Pstar < pr)
	{
		// The pattern is rarefaction-contact-rarefaction
		Pattern = 1;
	}
	else if(Pstar < pl && Pstar >= pr)
	{
		// The pattern is rarefaction-contact-shock
		Pattern = 2;
	}
	else if(Pstar >= pl && Pstar < pr)
	{
		// The pattern is shock-contact-rarefaction
		Pattern = 3 ;
	}
	else if(Pstar >= pl && Pstar >= pr)
	{
		// The pattern is shock-contact-shock
		Pattern = 4 ;
	}

	//std::cout<< "Pattern = " << Pattern << " pstar = " << Pstar << " ustar = " << ustar << std::endl;

	// Sample the solution

	for(int i=0; i<X.size(); i++)
	{
		double csi = (X[i]-xd)/time;

		VectOfDouble WExact(nVar);

		WExact = Sample(Wl,Wr,Pstar,ustar,Pattern,csi);

		UExact[i] = primitiveToconservative(WExact);
	}
}

double GuessPressure(const VectOfDouble& Wl, const VectOfDouble& Wr)
{
	double PGuess;

	double G1,G3,G4,G5,G6,G7;

	G1 = (Gamma-1)/(2*Gamma) ;
	G3 = 2*Gamma/(Gamma-1) ;
	G4 = 2/(Gamma-1) ;
	G5 = 2/(Gamma+1) ;
	G6 = (Gamma-1)/(Gamma+1) ;
	G7 = (Gamma-1)/2 ;

	double rhol = Wl[0];
	double rhor = Wr[0];

	double ul = Wl[1];
	double ur = Wr[1];

	double pl = Wl[2];
	double pr = Wr[2];

    double al = sqrt(Gamma*pl/rhol);
	double ar = sqrt(Gamma*pr/rhor);

    // Guess pressure based on PVRS RP Solver

    int Quser = 2;
    double cup,ppv,pmin,pmax,qmax;

    cup = 0.25*(rhol+rhor)*(al+ar);
    ppv = std::max(0.0,0.5*(pl+pr)+0.5*(ul-ur)*cup);

    pmin = std::min(pl,pr);
    pmax = std::max(pl,pr);

    qmax = pmax/pmin;

    bool Cond1 = (qmax <= Quser) && ((pmin <= ppv) && (pmax >= ppv));

    if(Cond1)
    {
		// Select PVRS RP Solver
		PGuess = ppv;
	}
	else
	{
		if(ppv < pmin)
		{
			// Select two-rarefaction RP Solver

			double pq,um,ptl,ptr;

			pq = pow(pl/pr,G1);
			um = (pq*ul/al+ur/ar+G4*(pq-1.0))/(pq/al+1.0/ar);

			ptl = 1.0+G7*(ul-um)/al;
			ptr = 1.0+G7*(um-ur)/ar;

			PGuess = 0.5*(pl*pow(ptl,G3)+pr*pow(ptr,G3));
		}
		else
		{
			// Select two-shock RP solver with PVRS as estimate

			double gel,ger;

			gel = sqrt((G5/rhol)/(G6*pl+ppv));
			ger = sqrt((G5/rhor)/(G6*pr+ppv));

			PGuess = (gel*pl+ger*pr-(ur-ul))/(gel+ger);
		}
	}

	return PGuess;
}

double fl_NR(const VectOfDouble Wl, const double& p)
{
	double fl;

	double G1,G4,G5,G6;

	G1 = (Gamma-1)/(2*Gamma) ;
	G4 = 2/(Gamma-1) ;
	G5 = 2/(Gamma+1) ;
	G6 = (Gamma-1)/(Gamma+1) ;

    double rhol = Wl[0];
	double ul   = Wl[1];
	double pl   = Wl[2];

	double al = sqrt(Gamma*pl/rhol);

	double Al = G5/rhol;
    double Bl = G6*pl;

    // Newton Raphson function

    fl = (p > pl) ? (p-pl)*sqrt(Al/(p+Bl)) : G4*al*(pow(p/pl,G1)-1.0);

	return fl;
}

double fr_NR(const VectOfDouble Wr, const double& p)
{
	double fr;

	double G1,G4,G5,G6;

	G1 = (Gamma-1)/(2*Gamma) ;
	G4 = 2/(Gamma-1) ;
	G5 = 2/(Gamma+1) ;
	G6 = (Gamma-1)/(Gamma+1) ;

    double rhor = Wr[0];
	double ur   = Wr[1];
	double pr   = Wr[2];

	double ar = sqrt(Gamma*pr/rhor);

	double Ar = G5/rhor;
    double Br = G6*pr;

    // Newton Raphson function

    fr = (p > pr) ? (p-pr)*sqrt(Ar/(p+Br)) : G4*ar*(pow(p/pr,G1)-1.0);

	return fr;
}

double f_NR(const VectOfDouble Wl, const VectOfDouble Wr, const double& p)
{
	double f;

	double ul = Wl[1];
	double ur = Wr[1];

	double fl = fl_NR(Wl,p);
	double fr = fr_NR(Wr,p);

	// Newton Raphson function

	f = fl+fr+(ur-ul);

	return f;

}

double df_NR(const VectOfDouble Wl, const VectOfDouble Wr, const double& p)
{
	double df;

	double G2,G5,G6;

	G2 = (Gamma+1)/(2*Gamma) ;
	G5 = 2/(Gamma+1) ;
	G6 = (Gamma-1)/(Gamma+1) ;

	double rhol = Wl[0];
	double rhor = Wr[0];

	double ul = Wl[1];
	double ur = Wr[1];

	double pl = Wl[2];
	double pr = Wr[2];

    double al = sqrt(Gamma*pl/rhol);
	double ar = sqrt(Gamma*pr/rhor);

	double Al = G5/rhol;
	double Ar = G5/rhor;

    double Bl = G6*pl;
    double Br = G6*pr;

    // Newton Raphson d_function

    double dfL = (p > pl) ? sqrt(Al/(Bl+p))*(1-0.5*(p-pl)/(Bl+p))
                          : 1.0/(rhol*al)*pow(p/pl,-G2);
    double dfR = (p > pr) ? sqrt(Ar/(Br+p))*(1-0.5*(p-pr)/(Br+p))
                          : 1.0/(rhor*ar)*pow(p/pr,-G2);

	df = dfL+dfR;

	return df;

}

VectOfDouble Sample(const VectOfDouble Wl, const VectOfDouble Wr, const double& pstar, const double& ustar, const int& Pattern, const double& csi)
{
	VectOfDouble Wexact(nVar);

	double G1,G2,G3,G4,G5,G6,G7 ;

	G1 = (Gamma-1)/(2*Gamma) ;
	G2 = (Gamma+1)/(2*Gamma) ;
	G3 = 2*Gamma/(Gamma-1) ;
	G4 = 2/(Gamma-1) ;
	G5 = 2/(Gamma+1) ;
	G6 = (Gamma-1)/(Gamma+1) ;
	G7 = (Gamma-1)/2 ;

    double rhol = Wl[0];
	double rhor = Wr[0];

	double ul = Wl[1];
	double ur = Wr[1];

	double pl = Wl[2];
	double pr = Wr[2];

    double al = sqrt(Gamma*pl/rhol);
	double ar = sqrt(Gamma*pr/rhor);

    // Sampling the solution

    double rholstar,alstar,Shl,Stl,Sl;
    double rhorstar,arstar,Shr,Str,Sr;

    if(Pattern == 1 || Pattern == 2)
    {
		rholstar = rhol*pow(pstar/pl,1.0/Gamma);
		alstar   = al*pow(pstar/pl,G1);
		Shl      = ul-al;
		Stl      = ustar-alstar ;
	}
	else
	{
		rholstar = rhol*(pstar/pl+G6)/(1.0+G6*pstar/pl);
        Sl       = ul-al*sqrt(G1+G2*pstar/pl);
	}

	if(Pattern == 1 || Pattern == 3)
	{
		rhorstar = rhor*pow(pstar/pr,1.0/Gamma);
        arstar   = ar*pow(pstar/pr,G1);
        Shr      = ur+ar;
		Str      = ustar+arstar;
	}
	else
	{
		rhorstar = rhor*(pstar/pr+G6)/(1.0+G6*pstar/pr) ;
        Sr       = ur+ar*sqrt(G1+G2*pstar/pr) ;
	}

    switch(Pattern)
    {
		case 1:
		{
			if(csi < Shl)
			{
				Wexact[0] = Wl[0];
				Wexact[1] = Wl[1];
				Wexact[2] = Wl[2];
		    }
		    else if(csi >= Shl && csi < Stl)
		    {
				Wexact[0] = rhol*pow(G5+G6/al*(ul-csi),G4);
				Wexact[1] = G5*(al+G7*ul+csi);
				Wexact[2] = pl*pow(G5+G6/al*(ul-csi),G3);
		    }
		    else if(csi >= Stl && csi < ustar)
		    {
				Wexact[0] = rholstar;
				Wexact[1] = ustar;
				Wexact[2] = pstar;
		    }
		    else if(csi >= ustar && csi < Str)
		    {
				Wexact[0] = rhorstar;
				Wexact[1] = ustar;
				Wexact[2] = pstar;
		    }
		    else if(csi >= Str && csi < Shr)
		    {
				Wexact[0] = rhor*pow(G5-G6/ar*(ur-csi),G4);
				Wexact[1] = G5*(-ar+G7*ur+csi);
				Wexact[2] = pr*pow(G5-G6/ar*(ur-csi),G3);
		    }
		    else if(csi >= Shr)
		    {
				Wexact[0] = Wr[0];
				Wexact[1] = Wr[1];
				Wexact[2] = Wr[2];
		    }

			break;
		}

		case 2:
		{
			if (csi < Shl)
			{
				Wexact[0] = Wl[0];
				Wexact[1] = Wl[1];
				Wexact[2] = Wl[2];
		    }
		    else if(csi >= Shl && csi < Stl)
		    {
				Wexact[0] = rhol*pow(G5+G6/al*(ul-csi),G4);
				Wexact[1] = G5*(al+G7*ul+csi);
				Wexact[2] = pl*pow(G5+G6/al*(ul-csi),G3);
		    }
		    else if(csi >= Stl && csi < ustar)
		    {
				Wexact[0] = rholstar;
				Wexact[1] = ustar;
				Wexact[2] = pstar;
		    }
		    else if(csi >= ustar && csi < Sr)
		    {
				Wexact[0] = rhorstar;
				Wexact[1] = ustar;
				Wexact[2] = pstar;
		    }
		    else if(csi >= Sr)
		    {
				Wexact[0] = Wr[0];
				Wexact[1] = Wr[1];
				Wexact[2] = Wr[2];
		    }

			break;
		}

		case 3:
		{
			if(csi < Sl)
			{
				Wexact[0] = Wl[0];
				Wexact[1] = Wl[1];
				Wexact[2] = Wl[2];
		    }
		    else if(csi >= Sl && csi < ustar)
		    {
				Wexact[0] = rholstar;
				Wexact[1] = ustar;
				Wexact[2] = pstar;
		    }
		    else if(csi >= ustar && csi < Str)
		    {
				Wexact[0] = rhorstar;
				Wexact[1] = ustar;
				Wexact[2] = pstar;
		    }
		    else if(csi >= Str && csi < Shr)
		    {
				Wexact[0] = rhor*pow(G5-G6/ar*(ur-csi),G4);
				Wexact[1] = G5*(-ar+G7*ur+csi);
				Wexact[2] = pr*pow(G5-G6/ar*(ur-csi),G3);
		    }
		    else if(csi >= Shr)
		    {
				Wexact[0] = Wr[0];
				Wexact[1] = Wr[1];
				Wexact[2] = Wr[2];
		    }

			break;
		}

		case 4:
		{
			if(csi < Sl)
			{
				Wexact[0] = Wl[0];
				Wexact[1] = Wl[1];
				Wexact[2] = Wl[2];
		    }

		    else if(csi >= Sl && csi < ustar)
		    {
				Wexact[0] = rholstar;
				Wexact[1] = ustar;
				Wexact[2] = pstar;
		    }
		    else if(csi >= ustar && csi < Sr )
		    {
				Wexact[0] = rhorstar;
				Wexact[1] = ustar;
				Wexact[2] = pstar;
		    }
		    else if(csi >= Sr)
		    {
				Wexact[0] = Wr[0];
				Wexact[1] = Wr[1];
				Wexact[2] = Wr[2];
		    }

			break;
		}
	}

	return Wexact;

}




/* ------------------------------------------------------------------ */
/*                                                                    */
/*                       END OF THE FILE                              */
/*                                                                    */
/* ------------------------------------------------------------------ */
