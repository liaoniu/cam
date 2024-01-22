#include "Ghost.H"



// Global Variables

extern double x0, x1, T, GammaL, GammaR, xI;
extern int nxCells;
extern double dx;
extern int nVar;
extern int Nghost_fluid;


// MAIN SOLVER

void solver(mat &u0, const double &Gamma, string method = "FORCE"){
    double t = 0;
    double dt = 0;
    mat flux;
    mat u = u0;
    flux.resize(nxCells+1);
    u = primitiveToConserved(u, Gamma);

    do{
        dt = calc_dt(u, Gamma);
        t += dt;
        //Transmissive BC
        transmissiveBC(u);
        //Dimensional Split 1
        mat u_temp = u;
        split(u, flux, dt, Gamma, method);
        //cout << t << endl;
    }while(t<T);

    u = conservedToPrimitive(u, Gamma);
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

double computePressureFromEoS(const double& rho, const double& e, const double &Gamma)
{
	// Ideal Gas Eos

	double p = rho*e*(Gamma-1.0);

	// Add here other EoS ....

	return p;
}

double computeInternalEnergyFromEoS(const double& rho, double& p, const double &Gamma)
{
	// Ideal Gas Eos

	double e = p/(rho*(Gamma-1.0));

	// Add here other EoS ....

	return e;
}

double computeSoundSpeedFromEoS(const double& rho, const double& p, const double &Gamma)
{
	// Ideal Gas Eos

	double a = sqrt(Gamma*p/rho);

	// Add here other EoS ....

	return a;
}


// Variables inter change

arr indi_primitiveToConserved(const arr &u_pri, const double &Gamma)
{
    arr u_con = u_pri;
    double rho, v, p;
    tie(rho, v, p) = tie(u_pri[0], u_pri[1], u_pri[2]);
    double rhoV = rho*v;
    double epsilon = computeInternalEnergyFromEoS(rho, p, Gamma);
    double E = rho*epsilon + rho*(v*v)/2;
    tie(u_con[0], u_con[1], u_con[2]) = tie(rho, rhoV, E);
    return u_con;
}

arr indi_conservedToPrimitive(const arr &u_con, const double &Gamma)
{
    arr u_pri = u_con;
    double rho, rhoV, E;
    tie(rho, rhoV, E) = tie(u_con[0], u_con[1], u_con[2]);
    double v = rhoV/rho;
    double epsilon = E/rho - 0.5*(v*v);
    double p = computePressureFromEoS(rho, epsilon, Gamma);
    tie(u_pri[0], u_pri[1], u_pri[2]) = tie(rho, v, p);
    return u_pri;
}

mat conservedToPrimitive(const mat &u, const double &Gamma)
{
    //u would not be changed
    mat u_pri = u;
    for (int i = 0; i != u.size(); i++){
        u_pri[i] = indi_conservedToPrimitive(u[i], Gamma);
    }
    return u_pri;
}

mat primitiveToConserved(const mat &u, const double &Gamma)
{
    //u would not be changed
    mat u_con = u;
    for (int i = 0; i != u.size(); i++){
        u_con[i] = indi_primitiveToConserved(u[i], Gamma);
    }
    return u_con;
}

double calc_dt(const mat &u_con, const double &Gamma)
{
    double delta = dx;
    double C = 0.8;
    double a_max = 0;
    int m = u_con.size();
    for (int i = 0; i != m; i++){
            double rho, v, p;
            arr u_pri_indi = indi_conservedToPrimitive(u_con[i], Gamma);
            tie(rho, v, p) = tie(u_pri_indi[0], u_pri_indi[1], u_pri_indi[2]);
            double a = sqrt(v*v) + sqrt(Gamma*p/rho);
            a_max = ((a_max < a) ? a : a_max);
    }
    double dt = C*delta/a_max;
    return dt;
}

arr f(const arr &u, const double &Gamma)
{
    double rho, v, p;
    double rhoV, E;
    arr u_pri = indi_conservedToPrimitive(u, Gamma);
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

tuple<mat, mat> halfTimeUpdate(mat &uL, mat &uR, double dt, const double &Gamma)
{
    int nxCells = uL.size()-2;
    for (int i = 0; i != nxCells+2; i++){
        arr fL = f(uL[i], Gamma);
        arr fR = f(uR[i], Gamma);
        for (int var = 0; var != nVar; var++){
            uL[i][var] = uL[i][var] - 0.5*dt/dx*(fR[var] - fL[var]);
            uR[i][var] = uR[i][var] - 0.5*dt/dx*(fR[var] - fL[var]);
        }
        
    }
    return tie(uL, uR);
}


// Get Flux

arr getXFlux(const arr &uL, const arr &uR, double &dt, const double &Gamma, string method)   //u_left and u_right
{ 
    
    arr flux, u_mid;
    if (method == "FORCE")
    {
        arr fL, fR, LF, RI;
        fL = f(uL, Gamma);
        fR = f(uR, Gamma);
        for (int var = 0; var != nVar; var++){
            LF[var] = 0.5*dx/dt*(uL[var] - uR[var]) + 0.5*(fL[var] + fR[var]);
            u_mid[var] = 0.5*(uL[var] + uR[var]) - 0.5*dt/dx*(fR[var] - fL[var]);
        }
        RI = f(u_mid, Gamma);
        for (int var = 0; var != nVar; var++){
            flux[var] = 0.5*(LF[var] + RI[var]);
        }
    }
    else if (method == "HLLC")
    {
        HLLCFlux(uL, uR, flux, Gamma);
    }
    return flux;
}

void HLLCFlux(const arr& Ql, const arr& Qr, arr& Fhllc, const double &Gamma)
{
	double Sl, Sr, S_star;

	WaveEstimates(Ql, Qr, Sl, Sr, Gamma);

	arr Fl, Fr, Fl_star, Fr_star;

	Fl = f(Ql, Gamma);
	Fr = f(Qr, Gamma);
	
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

	Wl = indi_conservedToPrimitive(Ql, Gamma);
	Wr = indi_conservedToPrimitive(Qr, Gamma);
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

void WaveEstimates(const arr& Ql, const arr& Qr, double& Sl, double& Sr, const double &Gamma)
{
	arr Wl, Wr;

	Wl = indi_conservedToPrimitive(Ql, Gamma);
	Wr = indi_conservedToPrimitive(Qr, Gamma);

	double rhol = Wl[0];
	double rhor = Wr[0];

	double ul = Wl[1];
	double ur = Wr[1];

	double pl = Wl[2];
	double pr = Wr[2];

	double al = computeSoundSpeedFromEoS(rhol,pl,Gamma);
	double ar = computeSoundSpeedFromEoS(rhor,pr,Gamma);

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




mat split(mat &u, mat &flux, double dt, const double &Gamma, string method)
{
    mat uL, uR;
    if (method == "SLIC"){
        tie(uL, uR) = dataReconstruct(u);
        halfTimeUpdate(uL, uR, dt, Gamma);
        //getXFLux
        for (int i = 0; i != nxCells+1; i++){
                flux[i] = getXFlux(uR[i], uL[i+1], dt, Gamma, "FORCE");
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
                flux[i] = getXFlux(u[i], u[i+1], dt, Gamma, method);
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





// Level Set Functions

// WENO

double divededD_pos(vector<double> &phi, int i)
{
    double D = (phi[i+1] - phi[i])/dx;
    return D;
}

double divededD_neg(vector<double> &phi, int i)
{
    double D = (phi[i] - phi[i-1])/dx;
    return D;
}

tuple<double, double> weno(vector<double> &phi, int i)
{
    double v1, v2, v3, v4, v5;
    double S1, S2, S3;
    double alpha1, alpha2, alpha3;
    double w1, w2, w3;
    double a;
    v1 = divededD_neg(phi, i-2);
    v2 = divededD_neg(phi, i-1);
    v3 = divededD_neg(phi, i);
    v4 = divededD_neg(phi, i+1);
    v5 = divededD_neg(phi, i+2);
    double phi_neg_1 = v1/3 - 7*v2/6 + 11*v3/6;
    double phi_neg_2 = -v2/6 - 5*v3/6 + v4/3;
    double phi_neg_3 = v3/3 + 5*v4/6 - v5/6;
    S1 = 13/12*pow(v1 - 2*v2 + v3, 2) + 1/4*pow(v1 - 4*v2 + 3*v3, 2);
    S2 = 13/12*pow(v2 - 2*v3 + v4, 2) + 1/4*pow(v2 - v4, 2);
    S3 = 13/12*pow(v3 - 2*v4 + v5, 2) + 1/4*pow(3*v3 - 4*v4 + v5, 2);
    alpha1 = 0.1/pow(S1 + 1e-6, 2);
    alpha2 = 0.6/pow(S2 + 1e-6, 2);
    alpha3 = 0.3/pow(S3 + 1e-6, 2);
    a = alpha1 + alpha2 + alpha3;
    w1 = alpha1/a;
    w2 = alpha2/a;
    w3 = alpha3/a;
    double phi_neg = w1*phi_neg_1 + w2*phi_neg_2 + w3*phi_neg_3;

    v1 = divededD_pos(phi, i+2);
    v2 = divededD_pos(phi, i+1);
    v3 = divededD_pos(phi, i);
    v4 = divededD_pos(phi, i-1);
    v5 = divededD_pos(phi, i-2);
    double phi_plus_1 = v1/3 - 7*v2/6 + 11*v3/6;
    double phi_plus_2 = -v2/6 - 5*v3/6 + v4/3;
    double phi_plus_3 = v3/3 + 5*v4/6 - v5/6;
    S1 = 13/12*pow(v1 - 2*v2 + v3, 2) + 1/4*pow(v1 - 4*v2 + 3*v3, 2);
    S2 = 13/12*pow(v2 - 2*v3 + v4, 2) + 1/4*pow(v2 - v4, 2);
    S3 = 13/12*pow(v3 - 2*v4 + v5, 2) + 1/4*pow(3*v3 - 4*v4 + v5, 2);
    alpha1 = 0.1/pow(S1 + 1e-6, 2);
    alpha2 = 0.6/pow(S2 + 1e-6, 2);
    alpha3 = 0.3/pow(S3 + 1e-6, 2);
    a = alpha1 + alpha2 + alpha3;
    w1 = alpha1/a;
    w2 = alpha2/a;
    w3 = alpha3/a;
    double phi_plus = w1*phi_plus_1 + w2*phi_plus_2 + w3*phi_plus_3;
    

    return tie(phi_plus, phi_neg);
}



// Initialize phi

void initializePhi(vector<double> &phi)
{
    for(int i = 0; i != phi.size(); i++)
    {
        double x = x0 + (i-0.5)*dx;
        phi[i] = x - xI;
    }
    transmissiveBC(phi);
}

void initializePhi_Helium(vector<double> &phi)
{
    for(int i = 0; i != phi.size(); i++)
    {
        double x = x0 + (i-0.5)*dx;
        phi[i] = 0.1 - fabs(x - 0.5);
    }
    // transmissiveBC(phi);
}

void updatePhi(vector<double> &phi, const mat &u, double &dt)
{
    vector<double> phi_new = phi;
    for(int i = 1; i != nxCells + 1; i++)
    {
        double v = u[i][1]/u[i][0];
        double Dphi;
        Dphi = ((v < 0) ? (phi[i+1] - phi[i])/dx : (phi[i] - phi[i-1])/dx);
        phi_new[i] = phi[i] - dt*v*Dphi;
    }
    phi = phi_new;
    transmissiveBC(phi);
}

void updatePhi(vector<double> &phi, const mat &u, const mat &u1, double &dt)
{
    vector<double> phi_new = phi;
    for(int i = 1; i != nxCells + 1; i++)
    {
        double v;
        if (phi[i] < 0)
            v = u[i][1]/u[i][0];
        else
            v = u1[i][1]/u1[i][0];
        double Dphi;
        double phi_plus, phi_neg;
        Dphi = ((v < 0) ? (phi[i+1] - phi[i])/dx : (phi[i] - phi[i-1])/dx);

        // tie(phi_plus, phi_neg) = weno(phi, i);
        // Dphi = ((v < 0) ? phi_plus : phi_neg);
        phi_new[i] = phi[i] - dt*v*Dphi;
    }

    phi = phi_new;

    // int index = positionOfInterface(phi);
    // int i = index;
    // double v = ((phi[index] < 0) ? u[index][1]/u[index][0] : u1[index][1]/u1[index][0]);
    // double Dphi, phi_plus, phi_minus;
    // Dphi = ((v < 0) ? (phi[i+1] - phi[i])/dx : (phi[i] - phi[i-1])/dx);

    // transmissiveBC(phi);
}

void reinitialization(vector<double> &phi)
{
    int index;
    double signLeft;
    for (int i = 0; i != phi.size()-1; i++)
    {
        if (phi[i]*phi[i+1] <= 0)
        {
            index = i;
            signLeft = ((phi[i] < 0) ? -1 : 1);
            break;
        }
    }
    int pt1 = index;
    int pt2 = index + 1;
    
    while (pt1 != 0)
    {
        phi[pt1-1] = phi[pt1] + signLeft*dx;
        pt1--;
    }
    while (pt2 != phi.size())
    {
        phi[pt2+1] = phi[pt2] - signLeft*dx;
        pt2++;
    }
    transmissiveBC(phi);
}

void fastSweeping(vector<double> &phi)
{
    bool implement = false;
    int i = 0;
    // Direction 1
    while (i != phi.size()-1)
    {
        if (implement)
        {
            double phi_guess;
            if (phi[i] >= 0)
            {
                double phi_x = min(phi[i-1], phi[i+1]);
                phi_guess = phi_x + dx;
            }
            else
            {
                double phi_x = max(phi[i-1], phi[i+1]);
                phi_guess = phi_x - dx;
            }
            if (fabs(phi_guess) < fabs(phi[i]))
            {
                phi[i] = phi_guess;
            }
        }

        if (phi[i]*phi[i+1] < 0)
        {
            implement = true;
            i = i+2;
        }
        else
            i++;
        
        if (i == phi.size()-1)
        {
            phi[i] = ((phi[i-1] < 0) ? phi[i-1] - dx : phi[i-1] + dx);
        }
    }

    // Direction 2
    implement = false;
    i = phi.size() - 1;
    while (i != 0)
    {
        if (implement)
        {
            double phi_guess;
            if (phi[i] >= 0)
            {
                double phi_x = min(phi[i-1], phi[i+1]);
                phi_guess = phi_x + dx;
            }
            else
            {
                double phi_x = max(phi[i-1], phi[i+1]);
                phi_guess = phi_x - dx;
            }
            if (fabs(phi_guess) < fabs(phi[i]))
            {
                phi[i] = phi_guess;
            }
        }

        if (phi[i]*phi[i-1] < 0)
        {
            implement = true;
            i = i-2;
        }
        else
            i--;
        
        if (i == 0)
        {
            phi[i] = ((phi[i+1] < 0) ? phi[i+1] - dx : phi[i+1] + dx);
        }
    }
    // transmissiveBC(phi);
}


// Phi Operations

int positionOfInterface(vector<double> &phi)
{
    int index;
    for (int i = 0; i != phi.size()-1; i++)
    {
        if (phi[i]*phi[i+1] <= 0)
        {
            index = i;
            return index;
        }
    }
    return phi.size()-1;
}

vector<int> posOfAllInterface(vector<double> &phi)
{
    vector<int> pos;
    for (int i = 0; i != phi.size()-1; i++)
    {
        if (phi[i]*phi[i+1] < 0)
        {
            pos.push_back(i);
        }
        else if (phi[i]*phi[i+1] == 0)
        {
            pos.push_back(i);
            i++;
        }
    }
    return pos;
}



void solver_LevelSet2Var(mat &u0, mat &u1, vector<double> &phi0, const double &Gamma, string method)
{
    double t = 0;
    double dt = 0;
    mat flux(nxCells+1);
    mat flux1(nxCells+1);
    mat u = u0;
    u = primitiveToConserved(u, Gamma);
    u1 = primitiveToConserved(u1, Gamma);
    vector<double> phi(nxCells+2);
    initializePhi(phi);

    do{
        dt = min(calc_dt(u, Gamma), calc_dt(u1, Gamma));
        t += min(dt, T-t);
        //Transmissive BC
        transmissiveBC(u);
        transmissiveBC(u1);
        //Dimensional Split 1
        split(u, flux, dt, Gamma, method);
        split(u1, flux1, dt, Gamma, method);
        updatePhi(phi, u, u1, dt);
        reinitialization(phi);
        //cout << t << endl;

    }while(t<T);

    u = conservedToPrimitive(u, Gamma);
    u1 = conservedToPrimitive(u1, Gamma);
    ofstream output1("rho.dat");
    ofstream output2("v.dat");
    ofstream output3("p.dat");
    ofstream output4("phi.dat");

    for (int i = 1; i != nxCells + 1; i++){
        double x = x0 + (i-0.5)*dx;
        output1 << x << " " << u[i][0] << " " << u1[i][0] << endl;
        output2 << x << " " << u[i][1] << " " << u1[i][1] << endl;
        output3 << x << " " << u[i][2] << " " << u1[i][2] << endl;
        output4 << x << " " << phi[i] << endl;
    }

    output1.close();
    output2.close();
    output3.close();
    output4.close();
}



