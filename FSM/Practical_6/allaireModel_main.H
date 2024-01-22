#include "allaireModel.H"



// Global Variables

extern double x0, x1, T, Gamma1, Gamma2, p_infty1, p_infty2;
extern int nxCells;
extern double dx, C;
extern int nVar;



// MAIN SOLVER

void solver(mat &u0, string method = "HLLC"){
    double t = 0;
    double dt = 0;
    mat flux;
    mat u = u0;
    flux.resize(nxCells+1);
    u = primitiveToConserved(u);


    do
    // {
    // for (int i = 0; i != 3; i++)
    {
        dt = min(calc_dt_Allaire(u), T-t);
        t += dt;
        //Transmissive BC
        transmissiveBC(u);
        //Dimensional Split 1
        split(u, flux, dt, method);
        //cout << t << endl;


    // }
    }while(t<T);

        // dt = calc_dt_Allaire(u);
        // t += dt;
        // //Transmissive BC
        // transmissiveBC(u);
        // //Dimensional Split 1
        // split(u, flux, dt, method);
        // //cout << t << endl;


    u = conservedToPrimitive(u);
    ofstream output1("a1.dat");
    ofstream output2("rho1.dat");
    ofstream output3("rho2.dat");
    ofstream output4("v.dat");
    ofstream output5("p.dat");
    ofstream output6("rho.dat");


    for (int i = 1; i != nxCells + 1; i++){
        double x = x0 + (i-0.5)*dx;
        output1 << x << " " << u[i][0] << endl;
        output2 << x << " " << u[i][1] << endl;
        output3 << x << " " << u[i][2] << endl;
        output4 << x << " " << u[i][3] << endl;
        output5 << x << " " << u[i][4] << endl;
        output6 << x << " " << u[i][0]*u[i][1] + (1-u[i][0])*u[i][2] << endl;
        }

    output1.close();
    output2.close();
    output3.close();
    output4.close();
    output5.close();
    output6.close();

}






// EoS Routine

double computePressureFromEoS(const double& rho, const double& e, const double &Gamma, const double &p_infty)
{
	// Ideal Gas Eos

	// Add here other EoS ....

    double p = rho*e*(Gamma-1.0) - Gamma*p_infty;

	return p;
}

double computeInternalEnergyFromEoS(const double& rho, double& p, const double &Gamma, const double &p_infty)
{
	// Ideal Gas Eos

	// Add here other EoS ....

    double e = (p+Gamma*p_infty)/((Gamma-1)*rho);

	return e;
}

double computeSoundSpeedFromEoS(const double& rho, const double& p, const double &Gamma)
{
	// Ideal Gas Eos

	double a = sqrt(Gamma*p/rho);

	// Add here other EoS ....



	return a;
}

double computeSoundSpeed_Allaire(const arr &u_con)
{
    double cs;
    double a1, rho1, rho2, v, p;
    double Gamma, p_infty;
    double rho, E;
    tie(Gamma, p_infty) = computeGammaPinfty_Alllaire(u_con[0]);
    arr u_pri = indi_conservedToPrimitive(u_con);
    E = u_con[4];
    tie(a1, rho1, rho2, v, p) = tie(u_pri[0], u_pri[1], u_pri[2], u_pri[3], u_pri[4]);
    rho = a1*rho1 + (1-a1)*rho2;
    double zeta = a1/(Gamma1 - 1) + (1-a1)/(Gamma2 - 1);
    double e = E/rho - 0.5*v*v;
    double h = e + p/rho;
    return sqrt(h/zeta);
    // return sqrt(Gamma*(p + p_infty)/rho);
}

tuple<double, double> computeGammaPinfty_Alllaire(double a1)
{
    double zeta = a1/(Gamma1 - 1) + (1-a1)/(Gamma2 - 1);
    double Gamma = 1 + 1/zeta;
    double p_infty = (Gamma-1)/Gamma*(a1*p_infty1*Gamma1/(Gamma1-1) + (1-a1)*p_infty2*Gamma2/(Gamma2-1));
    return tie(Gamma, p_infty);
}



// Variables inter change

arr indi_primitiveToConserved(const arr &u_pri)
{
    double Gamma, p_infty;
    tie(Gamma, p_infty) = computeGammaPinfty_Alllaire(u_pri[0]);
    arr u_con = u_pri;
    double a1, rho1, rho2, v, p;
    tie(a1, rho1, rho2, v, p) = tie(u_pri[0], u_pri[1], u_pri[2], u_pri[3], u_pri[4]);
    double a2 = 1 - a1;
    double rho = a1*rho1 + (1-a1)*rho2;
    double rhov = rho*v;
    double a1rho1 = a1*rho1;
    double a2rho2 = (1-a1)*rho2;
    double epsilon = computeInternalEnergyFromEoS(rho, p, Gamma, p_infty);
    double E = rho*epsilon + rho*(v*v)/2;
    tie(u_con[0], u_con[1], u_con[2], u_con[3], u_con[4]) = tie(a1, a1rho1, a2rho2, rhov, E);
    return u_con;
}

arr indi_conservedToPrimitive(const arr &u_con)
{
    double Gamma, p_infty;
    tie(Gamma, p_infty) = computeGammaPinfty_Alllaire(u_con[0]);
    arr u_pri = u_con;
    double a1, a1rho1, a2rho2, rhov, E;
    tie(a1, a1rho1, a2rho2, rhov, E) = tie(u_con[0], u_con[1], u_con[2], u_con[3], u_con[4]);
    double rho1 = a1rho1/a1;
    double rho2 = a2rho2/(1-a1);
    double rho = a1rho1 + a2rho2;
    double v = rhov/rho;
    double epsilon = E/rho - 0.5*(v*v);
    double p = computePressureFromEoS(rho, epsilon, Gamma, p_infty);
    tie(u_pri[0], u_pri[1], u_pri[2], u_pri[3], u_pri[4]) = tie(a1, rho1, rho2, v, p);
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

double calc_dt(const mat &u_con, const double &Gamma)
{
    double delta = dx;
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

double calc_dt_Allaire(const mat &u_con)
{
    double a_max = 0;
    for (int i = 0; i != u_con.size(); i++)
    {
        double cs;
        double a1, rho1, rho2, v, p;
        arr u_pri = indi_conservedToPrimitive(u_con[i]);
        tie(a1, rho1, rho2, v, p) = tie(u_pri[0], u_pri[1], u_pri[2], u_pri[3], u_pri[4]);
        cs = computeSoundSpeed_Allaire(u_con[i]);
        double a = sqrt(v*v) + cs;
        a_max = max(a, a_max);
    }
    double dt = C*dx/a_max;
    return dt;
}

arr f(const arr &u)
{
    double a1, rho1, rho2, v, p;
    double a1rho1, a2rho2, rhov, E;
    arr u_pri = indi_conservedToPrimitive(u);
    tie(a1, a1rho1, a2rho2, rhov, E) = tie(u[0], u[1], u[2], u[3], u[4]);
    tie(a1, rho1, rho2, v, p) = tie(u_pri[0], u_pri[1], u_pri[2], u_pri[3], u_pri[4]);
    arr f_res;
    f_res[0] = v*a1;
    f_res[1] = a1rho1*v;
    f_res[2] = a2rho2*v;
    f_res[3] = rhov*v + p;
    f_res[4] = (E+p)*v;
    return f_res;
}



// Get Flux

arr getXFlux(const arr &uL, const arr &uR, double &dt, string method, double &S)   //u_left and u_right
{ 
    arr flux, u_mid;
    if (method == "FORCE")
    {
        // arr fL, fR, LF, RI;
        // fL = f(uL, Gamma);
        // fR = f(uR, Gamma);
        // for (int var = 0; var != nVar; var++){
        //     LF[var] = 0.5*dx/dt*(uL[var] - uR[var]) + 0.5*(fL[var] + fR[var]);
        //     u_mid[var] = 0.5*(uL[var] + uR[var]) - 0.5*dt/dx*(fR[var] - fL[var]);
        // }
        // RI = f(u_mid, Gamma);
        // for (int var = 0; var != nVar; var++){
        //     flux[var] = 0.5*(LF[var] + RI[var]);
        // }
    }
    else if (method == "HLLC")
    {
        HLLCFlux(uL, uR, flux, S);
    }
    return flux;
}

void HLLCFlux(const arr& Ql, const arr& Qr, arr& Fhllc, double &S_star)
{
	double Sl, Sr;

	WaveEstimates(Ql, Qr, Sl, Sr);

	arr Fl, Fr, Fl_star, Fr_star;

	Fl = f(Ql);
	Fr = f(Qr);
	
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

	
	double a1l, rho1l, rho2l, ul, pl;
    double a1r, rho1r, rho2r, ur, pr;
    double El, Er;
	arr Wl, Wr;

    El = Ql[4];
    Er = Qr[4];
	Wl = indi_conservedToPrimitive(Ql);
	Wr = indi_conservedToPrimitive(Qr);
	tie(a1l, rho1l, rho2l, ul, pl) = tie(Wl[0], Wl[1], Wl[2], Wl[3], Wl[4]);
	tie(a1r, rho1r, rho2r, ur, pr) = tie(Wr[0], Wr[1], Wr[2], Wr[3], Wr[4]);
    double rhol = a1l*rho1l + (1-a1l)*rho2l;
    double rhor = a1r*rho1r + (1-a1r)*rho2r;
	S_star = (pr - pl + rhol*ul*(Sl - ul) - rhor*ur*(Sr - ur))/(rhol*(Sl - ul) - rhor*(Sr - ur));
	

	arr Ql_star, Qr_star;			// Change if move to higher dimensions


    // Ql_star[0] = a1l;
    // Qr_star[0] = a1r;
    Ql_star[0] = a1l*(Sl - ul)/(Sl - S_star);
    Qr_star[0] = a1r*(Sr - ur)/(Sr - S_star);
    Ql_star[1] = a1l*rho1l*(Sl - ul)/(Sl - S_star);
    Qr_star[1] = a1r*rho1r*(Sr - ur)/(Sr - S_star);
    Ql_star[2] = (1-a1l)*rho2l*(Sl - ul)/(Sl - S_star);
    Qr_star[2] = (1-a1r)*rho2r*(Sr - ur)/(Sr - S_star);
    Ql_star[3] = rhol*(Sl - ul)/(Sl - S_star)*S_star;
    Qr_star[3] = rhor*(Sr - ur)/(Sr - S_star)*S_star;
    Ql_star[4] = rhol*(Sl - ul)/(Sl - S_star)*(El/rhol + (S_star - ul)*(S_star + pl/rhol/(Sl - ul)));
    Qr_star[4] = rhor*(Sr - ur)/(Sr - S_star)*(Er/rhor + (S_star - ur)*(S_star + pr/rhor/(Sr - ur)));

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


// Wave Estimate

void WaveEstimates(const arr& Ql, const arr& Qr, double& Sl, double& Sr)
{
	arr Wl, Wr;

	Wl = indi_conservedToPrimitive(Ql);
	Wr = indi_conservedToPrimitive(Qr);

    double a1l, rho1l, rho2l, ul, pl;
    double a1r, rho1r, rho2r, ur, pr;
    tie(a1l, rho1l, rho2l, ul, pl) = tie(Wl[0], Wl[1], Wl[2], Wl[3], Wl[4]);
	tie(a1r, rho1r, rho2r, ur, pr) = tie(Wr[0], Wr[1], Wr[2], Wr[3], Wr[4]);

	double rhol = a1l*rho1l + (1-a1l)*rho2l;
    double rhor = a1r*rho1r + (1-a1r)*rho2r;


	double al = computeSoundSpeed_Allaire(Ql);
	double ar = computeSoundSpeed_Allaire(Qr);

	// Pressure–Based Wave Speed Estimates (ideal gases)

	double ql,qr;

	// Two-rarefaction Riemann solver TRRS for computing Pstar
    double GammaL, GammaR;
    tie(GammaL, ignore) = computeGammaPinfty_Alllaire(Ql[0]);
    tie(GammaR, ignore) = computeGammaPinfty_Alllaire(Qr[0]);
    double Gamma = (GammaL + GammaR)/2;

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

	Sl = ul-al;
	Sr = ur+ar;

	// Add here different Wave speed estimates ....
}




mat split(mat &u, mat &flux, double dt, string method)
{
    mat uL, uR;
    if (method == "SLIC"){
        // tie(uL, uR) = dataReconstruct(u);
        // halfTimeUpdate(uL, uR, dt, Gamma);
        // //getXFLux
        // for (int i = 0; i != nxCells+1; i++){
        //         flux[i] = getXFlux(uR[i], uL[i+1], dt, Gamma, "FORCE");
        // }
        // for (int i = 1; i != nxCells+1; i++){
        //     for (int k = 0; k != 3; k++){
        //         u[i][k] = u[i][k] - dt/dx*(flux[i][k] - flux[i-1][k]);
        //     }
        // }
    }
    else {
        // checkAlpha(u);
        //getXFLux
        vector<double> S(nxCells+1);
        // getXFlux(u[73], u[74], dt, method, S[73]);
        for (int i = 0; i != nxCells+1; i++){
            flux[i] = getXFlux(u[i], u[i+1], dt, method, S[i]);
        }
        for (int i = 1; i != nxCells+1; i++){
            double a1n = u[i][0];
            for (int k = 0; k != nVar; k++){
                u[i][k] = u[i][k] - dt/dx*(flux[i][k] - flux[i-1][k]);
            }
            // checkAlpha(u[i]);
            for (int j = 0; j != 5; j++)
            {
                u[i][0] = u[i][0] + dt/dx/5*a1n*(S[i] - S[i-1]);
                a1n = u[i][0];
            }

        }
        // checkAlpha(u);
    }

    //Need ensure boundary condition here
    transmissiveBC(u);
    return u;
}


// Boundary Condition

void transmissiveBC(mat &u)
{
    for (int var = 0; var != nVar; var++){
        u[0][var] = u[1][var];
        u[nxCells+1][var] = u[nxCells][var];
    }
}

void transmissiveBC(vector<double> &phi)
{
    int end = phi.size()-1;
    phi[0] = phi[1];
    phi[end] = phi[end-1];
}

void checkAlpha(mat &u)
{
    for (int i = 0; i != u.size(); i++)
    {
        checkAlpha(u[i]);
    }
}

void checkAlpha(arr &u)
{
    u[0] = ((u[0] > 1) ? 1-1e-2 : u[0]);
    u[0] = ((u[0] < 0) ? 1e-2 : u[0]);
}



