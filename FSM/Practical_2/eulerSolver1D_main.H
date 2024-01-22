#ifndef EULERSOLVER1D_H
#define EULERSOLVER1D_H

#include "eulerSolver1D.H"

// Global Variables

extern double x0, x1, T, Gamma;
extern int nxCells;
extern double dx;
extern int nVar;

// MAIN SOLVER

void solver(mat &u0, string method = "FORCE"){
    double t = 0;
    double dt = 0;
    mat flux;
    mat u = u0;
    flux.resize(nxCells+1);
    u = primitiveToConserved(u);

    do{
        dt = calc_dt(u);
        t += dt;
        //Transmissive BC
        transmissiveBC(u);
        //Dimensional Split 1
        mat u_temp = u;
        split(u, flux, dt, method);
        //cout << t << endl;
    }while(t<T);

    u = conservedToPrimitive(u);
    ofstream output1("rho.dat");
    ofstream output2("v.dat");
    ofstream output3("p.dat");


    for (int i = 1; i != nxCells + 1; i++){
        double x = x0 + (i-0.5)*dx;
        output1 << x << " " << u[i][0] << endl;
        output2 << x << " " << u[i][1] << endl;
        output3 << x << " " << u[i][2] << endl;
        }

    output1.close();
    output2.close();
    output3.close();

}






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


// SLIC Routine

double zeta(double delta_minus, double delta_plus, string limiter = "Minbee")
{
    if (limiter == "Minbee"){
        if (delta_plus == 0){
            return 0;
        }
        else{
            double r = delta_minus/delta_plus;
            return ((r <= 0) ? 0 : ((r <= 1) ? r : ((1 < 2/(1+r)) ? 1 : 2/(1+r))));
        }
    }
    return 0;
}

tuple<mat, mat> dataReconstruct(const mat &u)
{
    mat uL = u, uR = u;
    int nxCells = u.size()-2;
    double w = 0;
    for (int i = 1; i != nxCells+1; i++){
        for (int var = 0; var != nVar; var++){
            double delta_minus = u[i][var] - u[i-1][var];
            double delta_plus = u[i+1][var] - u[i][var];
            double delta = 0.5*(1+w)*delta_minus + 0.5*(1-w)*delta_plus;
            uL[i][var] = u[i][var] - 0.5*zeta(delta_minus, delta_plus)*delta;
            uR[i][var] = u[i][var] + 0.5*zeta(delta_minus, delta_plus)*delta;
        }
    }

    return tie(uL, uR);
}

tuple<mat, mat> halfTimeUpdate(mat &uL, mat &uR, double dt)
{
    int nxCells = uL.size()-2;
    for (int i = 0; i != nxCells+2; i++){
        arr fL = f(uL[i]);
        arr fR = f(uR[i]);
        for (int var = 0; var != nVar; var++){
            uL[i][var] = uL[i][var] - 0.5*dt/dx*(fR[var] - fL[var]);
            uR[i][var] = uR[i][var] - 0.5*dt/dx*(fR[var] - fL[var]);
        }
        
    }
    return tie(uL, uR);
}


// Get Flux

arr getXFlux(const arr &uL, const arr &uR, double &dt, string method)   //u_left and u_right
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
    else if (method == "HLLC")
    {
        HLLCFlux(uL, uR, flux);
    }
    return flux;
}

void HLLCFlux(const arr& Ql, const arr& Qr, arr& Fhllc)
{
	double Sl, Sr, S_star;

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

	
	double pl, pr, rhol, rhor, ul, ur;
	arr Wl, Wr;

	Wl = indi_conservedToPrimitive(Ql);
	Wr = indi_conservedToPrimitive(Qr);
	tie(rhol, ul, pl) = tie(Wl[0], Wl[1], Wl[2]);
	tie(rhor, ur, pr) = tie(Wr[0], Wr[1], Wr[2]);
	S_star = (pr - pl + rhol*ul*(Sl - ul) - rhor*ur*(Sr - ur))/(rhol*(Sl - ul) - rhor*(Sr - ur));
	

	arr D_star = {0, 1, S_star}, Ql_star, Qr_star;			// Change if move to higher dimensions
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


// Wave Estimate

void WaveEstimates(const arr& Ql, const arr& Qr, double& Sl, double& Sr)
{
	arr Wl, Wr;

	Wl = indi_conservedToPrimitive(Ql);
	Wr = indi_conservedToPrimitive(Qr);

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




mat split(mat &u, mat &flux, double dt, string method)
{
    mat uL, uR;
    if (method == "SLIC"){
        tie(uL, uR) = dataReconstruct(u);
        halfTimeUpdate(uL, uR, dt);
        //getXFLux
        for (int i = 0; i != nxCells+1; i++){
                flux[i] = getXFlux(uR[i], uL[i+1], dt, "FORCE");
        }
        for (int i = 1; i != nxCells+1; i++){
            for (int k = 0; k != 3; k++){
                u[i][k] = u[i][k] - dt/dx*(flux[i][k] - flux[i-1][k]);
            }
        }
    }
    else {
        //getXFLux
        for (int i = 0; i != nxCells+1; i++){
                flux[i] = getXFlux(u[i], u[i+1], dt, method);
        }
        for (int i = 1; i != nxCells+1; i++){
            for (int k = 0; k != 3; k++){
                u[i][k] = u[i][k] - dt/dx*(flux[i][k] - flux[i-1][k]);
            }

        }
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



#endif