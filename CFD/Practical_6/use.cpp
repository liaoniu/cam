
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;
using mat = vector<vector<double>>;
using myTuple = tuple<vector<double>, vector<double>, vector<double>>;

tuple<double, double, double> indi_primitiveToConserved(double rho, double v, double p, double gamma){
    double rhoV = rho*v;
    double epsilon = p/(gamma - 1)/rho;
    double E = rho*epsilon + rho*v*v/2;
    return tie(rho, rhoV, E);
}

tuple<double, double, double> indi_conservedToPrimitive(double rho, double rhov, double E, double gamma){
    double v = rhov/rho;
    double epsilon = E/rho - 0.5*v*v;
    double p = (gamma-1)*rho*epsilon;
    return tie(rho, v, p);
}

myTuple primitiveToConserved(vector<double> rho, vector<double> v, vector<double> p, double gamma){
    vector<double> rho_vec, rhov_vec, E_vec;
    for (int i = 0; i != rho.size(); i++){
        double r, rhov, E;
        tie(r,rhov,E) = indi_primitiveToConserved(rho[i], v[i], p[i], gamma);
        rho_vec.push_back(r);
        rhov_vec.push_back(rhov);
        E_vec.push_back(E);
    }
    return tie(rho_vec, rhov_vec, E_vec);
}

myTuple conservedToPrimitive(vector<double> rho, vector<double> rhov, vector<double> E, double gamma){
    vector<double> rho_vec, v_vec, p_vec;
    for (int i = 0; i != rho.size(); i++){
        double r, v, p;
        tie(r, v, p) = indi_conservedToPrimitive(rho[i], rhov[i], E[i], gamma);
        rho_vec.push_back(r);
        v_vec.push_back(v);
        p_vec.push_back(p);
    }
    return tie(rho_vec, v_vec, p_vec);
}

double calc_dt(vector<double> rho, vector<double> rhov, vector<double> E, double gamma, double dx){
    vector<double> cs, v, p;
    tie(rho, v, p) = conservedToPrimitive(rho, rhov, E, gamma);
    for (int i = 0; i != rho.size(); i++){
        cs.push_back(sqrt(gamma*p[i]/rho[i]));
    }
    double a_max = 0;
    for (int i = 0; i != v.size(); i++){
        a_max = ((a_max < fabs(v[i]) + cs[i]) ? fabs(v[i]) + cs[i] : a_max);
    }
    return 0.8*dx/a_max;
}

tuple<double, double, double> f(double rho, double rhov, double E, double gamma){
    double f1, f2, f3, v, p;
    v = rhov/rho;
    p = (gamma - 1)*(E - 0.5*rhov*v);
    f1 = rhov;
    f2 = rhov*v + p;
    f3 = (E + p)*v;
    return tie(f1, f2, f3);
}

tuple<double, double> fK_and_dfK(double p_star, tuple<double, double, double> uK, double gamma){
    double rhoK, vK, pK;
    tie(rhoK, vK, pK) = uK;
    double AK = 2/((gamma + 1)*rhoK);
    double BK = (gamma - 1)/(gamma + 1)*pK;
    double fK, dfK;
    if (p_star > pK){
        fK = (p_star - pK)*sqrt(AK/(p_star + BK));
        dfK = sqrt(AK/(BK + p_star))*(1 - (p_star - pK)/2/(BK + p_star));
    }
    else{
        double cs = sqrt(gamma*pK/rhoK);
        fK = 2*cs/(gamma - 1)*(pow(p_star/pK, (gamma - 1)/2/gamma) - 1);
        dfK = 1/rhoK/cs*pow(p_star/pK, -(gamma+1)/2/gamma);
    }
    return tie(fK, dfK);
}

tuple<double, double> f_and_df(double p_star, tuple<double, double, double> uL, tuple<double, double, double> uR, double gamma){
    double f, df;
    double vL, vR, dv;
    tie(ignore, vL, ignore) = uL;
    tie(ignore, vR, ignore) = uR;
    dv = vR - vL;
    double fL, fR, dfL, dfR;
    tie(fL, dfL) = fK_and_dfK(p_star, uL, gamma);
    tie(fR, dfR) = fK_and_dfK(p_star, uR, gamma);
    f = fL + fR + dv;
    df = dfL + dfR;
    return tie(f, df);
}


tuple<double, double, double, double, double, bool, bool> compute_stars_shock(tuple<double, double, double> uL, tuple<double, double, double> uR, double gamma, double p_star){
    double v_star, rhoL_star, rhoR_star, SL = 0, SR = 0;
    double rhoL, vL, pL, rhoR, vR, pR;
    bool left_shock = false, right_shock = false;
    tie(rhoL, vL, pL) = uL;
    tie(rhoR, vR, pR) = uR;
    double csL, csR;
    csL = sqrt(gamma*pL/rhoL);
    csR = sqrt(gamma*pR/rhoR);
    if (p_star > pL){
        SL = vL - csL*sqrt((gamma+1)/2/gamma*p_star/pL + (gamma-1)/2/gamma);
        rhoL_star = rhoL*((p_star/pL + (gamma-1)/(gamma+1))/((gamma - 1)/(gamma+1)*(p_star/pL)+1));
        v_star = (1 - rhoL/rhoL_star)*SL + vL*rhoL/rhoL_star;
        left_shock = true;
    }
    else{
        rhoL_star = rhoL*pow(p_star/pL, 1/gamma);
        v_star = vL - 2*csL/(gamma - 1)*(pow(p_star/pL, (gamma - 1)/2/gamma) - 1);
    }
    if (p_star > pR){
        SR = vR + csR*sqrt((gamma+1)/2/gamma*p_star/pR + (gamma-1)/2/gamma);
        rhoR_star = rhoR*((p_star/pR + (gamma-1)/(gamma+1))/((gamma - 1)/(gamma+1)*(p_star/pR)+1));
        v_star = (1 - rhoR/rhoR_star)*SR + vR*rhoR/rhoR_star;
        right_shock = true;
    }
    else{
        rhoR_star = rhoR*pow(p_star/pR, 1/gamma);
        v_star = vR + 2*csR/(gamma - 1)*(pow(p_star/pR, (gamma - 1)/2/gamma) - 1);
    }

    return tie(v_star, rhoL_star, rhoR_star, SL, SR, left_shock, right_shock);
}

tuple<double, double, double> primitiveInRarefaction(double gamma, tuple<double, double, double> uK, bool left){
    double rhoK, vK, pK;
    double rhoK_raref, vK_raref, pK_raref;
    tie(rhoK, vK, pK) = uK;
    double csK = sqrt(gamma*pK/rhoK);
    if (left){
        rhoK_raref = rhoK*pow(2/(gamma+1) + (gamma-1)/(gamma+1)/csK*(vK), 2/(gamma - 1));
        vK_raref = 2/(gamma + 1)*(csK + (gamma - 1)/2*vK);
        pK_raref = pK*pow(2/(gamma+1) + (gamma-1)/(gamma+1)/csK*(vK), 2*gamma/(gamma - 1));
    }
    else{
        rhoK_raref = rhoK*pow(2/(gamma+1) - (gamma-1)/(gamma+1)/csK*(vK), 2/(gamma - 1));
        vK_raref = 2/(gamma + 1)*(-csK + (gamma - 1)/2*vK);
        pK_raref = pK*pow(2/(gamma+1) - (gamma-1)/(gamma+1)/csK*(vK), 2*gamma/(gamma - 1));
    }
    return tie(rhoK_raref, vK_raref, pK_raref);
}

tuple<double, double, double> calc_central_u(tuple<double, double, double> uL, tuple<double, double, double> uR, double gamma, double p_star){
    double rhoL, vL, pL, rhoR, vR, pR, v_star, rhoL_star, rhoR_star, SL, SR;
    bool left_shock, right_shock;
    tie(rhoL, vL, pL) = uL;
    tie(rhoR, vR, pR) = uR;
    tie(v_star, rhoL_star, rhoR_star, SL, SR, left_shock, right_shock) = compute_stars_shock(uL, uR, gamma, p_star);
    double csL_star = sqrt(gamma*p_star/rhoL_star);
    double csR_star = sqrt(gamma*p_star/rhoR_star);
    double csL = sqrt(gamma*pL/rhoL);
    double csR = sqrt(gamma*pR/rhoR);
    
    if (left_shock){
        if (SL > 0){
            return tie(rhoL, vL, pL);
        }
    }
    else{
        if (vL - csL > 0){
            return tie(rhoL, vL, pL);
        }
    }
    if (right_shock){
        if (SR < 0){
            return tie(rhoR, vR, pR);
        }
    }
    else{
        if (vR + csR < 0){
            return tie(rhoR, vR, pR);
        }
    }
    if (!left_shock && (v_star - csL_star > 0)){
        return primitiveInRarefaction(gamma, uL, true);
    }
    if (!right_shock && (v_star + csR_star < 0)){
        return primitiveInRarefaction(gamma, uR, false);
    }
    else{
        if (v_star < 0){
            return tie(rhoR_star, v_star, p_star);
        }
        else{
            return tie(rhoL_star, v_star, p_star);
        }
    }
    
}



double compute_p_star(tuple<double, double, double> uL, tuple<double, double, double> uR, double gamma, double epsilon = 1e-3){
    double rhoL, vL, pL, rhoR, vR, pR;
    tie(rhoL, vL, pL) = uL;
    tie(rhoR, vR, pR) = uR;
    double csL = sqrt(gamma*pL/rhoL);
    double csR = sqrt(gamma*pR/rhoR);
    double p_init, p_PV;
    p_PV = (pR + pL)/2 - 1/8*(vR - vL)*(rhoL + rhoR)*(csL + csR);
    p_init = pow((csL+ csR - (gamma - 1)*(vR - vL)/2)/(csL/pow(pL, (gamma - 1)/2/gamma) + csR/pow(pR, (gamma - 1)/2/gamma)), 2*gamma/(gamma - 1));
    double p_new = p_init;
    double p_old;
    double f, df;
    double err;
    do{
        p_old = p_new;
        tie(f, df) = f_and_df(p_new, uL, uR, gamma);
        p_new = p_old - f/df;
        err = fabs(p_new - p_old)/2/fabs(p_new + p_old);
    }while(err > epsilon);
    return p_new;
}


tuple<double, double, double> getFlux(double rhoL, double rhoR, double rhovL, double rhovR, double EL, double ER, double gamma, double dx, double dt, string method){
    vector<double> flux, LF, RI, mid_u;
    vector<double> fL, fR;
    flux.resize(3);
    LF.resize(3);
    RI.resize(3);
    fL.resize(3);
    fR.resize(3);
    mid_u.resize(3);
    if (method == "Lax_Fredrichs"){
        tie(fL[0], fL[1], fL[2]) = f(rhoL, rhovL, EL, gamma);
        tie(fR[0], fR[1], fR[2]) = f(rhoR, rhovR, ER, gamma);
        LF[0] = 0.5*(dx/dt)*(rhoL - rhoR) + 0.5*(fL[0] + fR[0]);
        LF[1] = 0.5*(dx/dt)*(rhovL - rhovR) + 0.5*(fL[1] + fR[1]);
        LF[2] = 0.5*(dx/dt)*(EL - ER) + 0.5*(fL[2] + fR[2]);
        mid_u[0] = 0.5*(rhoR + rhoL) - 0.5*dt/dx*(fR[0] - fL[0]);
        mid_u[1] = 0.5*(rhovR + rhovL) - 0.5*dt/dx*(fR[1] - fL[1]);
        mid_u[2] = 0.5*(ER + EL) - 0.5*dt/dx*(fR[2] - fL[2]);
        tie(RI[0], RI[1], RI[2]) = f(mid_u[0], mid_u[1], mid_u[2], gamma);
        for(int i = 0; i != 3; i++){
            flux[i] = 0.5*(LF[i] + RI[i]);
        }
    }

    if (method == "Godunov"){
        double vL, pL, vR, pR, rho_mid, v_mid, p_mid;
        tie(rhoL, vL, pL) = indi_conservedToPrimitive(rhoL, rhovL, EL, gamma);
        tie(rhoR, vR, pR) = indi_conservedToPrimitive(rhoR, rhovR, ER, gamma);
        double p_star = compute_p_star(tie(rhoL, vL, pL), tie(rhoR, vR, pR), gamma);
        tie(rho_mid, v_mid, p_mid) = calc_central_u(tie(rhoL, vL, pL), tie(rhoR, vR, pR), gamma, p_star);
        tie(mid_u[0], mid_u[1], mid_u[2]) = indi_primitiveToConserved(rho_mid, v_mid, p_mid, gamma);
        tie(flux[0], flux[1], flux[2]) = f(mid_u[0], mid_u[1], mid_u[2], gamma);
    }
    return tie(flux[0], flux[1], flux[2]);
}



void solver(mat u0, double T, double x0, double x1, int nCells = 100, double gamma = 1.4, string method = "Lax_Fredrichs"){
    double dx = (x1-x0)/nCells;
    double t = 0;
    double dt = 0;
    vector<double> rho, rhov, v, p, E, rho_new, v_new, p_new, rhov_new, E_new;
    rho.resize(nCells+2);
    rhov.resize(nCells+2);
    v.resize(nCells+2);
    p.resize(nCells+2);
    E.resize(nCells+2);
    tie(rho, v, p) = tie(u0[0], u0[1], u0[2]);
    rho_new.resize(nCells+2);
    rhov_new.resize(nCells+2);
    v_new.resize(nCells+2);
    p_new.resize(nCells+2);
    E_new.resize(nCells+2);
    vector<double> flux_1, flux_2, flux_3;
    flux_1.resize(nCells+1);
    flux_2.resize(nCells+1);
    flux_3.resize(nCells+1);
    tie(rho, rhov, E) = primitiveToConserved(rho, v, p, gamma);

    
        dt = calc_dt(rho, rhov, E, gamma, dx);
        t += dt;
        rho[0] = rho[1];
        rho[nCells+1] = rho[nCells];
        rhov[0] = rhov[1];
        rhov[nCells+1] = rhov[nCells];
        E[0] = E[1];
        E[nCells+1] = E[nCells];
        for (int i = 0; i != flux_1.size(); i++){
            tie(flux_1[i], flux_2[i], flux_3[i]) = getFlux(rho[i], rho[i+1], rhov[i], rhov[i+1], E[i], E[i+1], gamma, dx, dt, method);
            
        }
        for (int i = 1; i != nCells+1; i++){
            rho_new[i] = rho[i] - dt/dx*(flux_1[i] - flux_1[i-1]);
            rhov_new[i] = rhov[i] - dt/dx*(flux_2[i] - flux_2[i-1]);
            E_new[i] = E[i] - dt/dx*(flux_3[i] - flux_3[i-1]);
        }
        rho = rho_new;
        rhov = rhov_new;
        E = E_new;
        cout << t << endl;
        

    tie(rho, v, p) = conservedToPrimitive(rho, rhov, E, gamma);
    ofstream output1("rho.dat");
    ofstream output2("v.dat");
    ofstream output3("p.dat");
    for (int i = 1; i != nCells + 1; i++){
        double x = x0 + (i-0.5)*dx;
        output1 << x << " " << rho[i] << endl;
        output2 << x << " " << v[i] << endl;
        output3 << x << " " << p[i] << endl;
    }
    output1.close();
    output2.close();
    output3.close();


}








int main()
{
    double x0 = 0, x1 = 1;
    int nCells = 2;
    double dx = (x1 - x0)/nCells;
    double gamma = 1.4;
    vector<double> rho, v, p;
    rho.resize(nCells+2);
    v.resize(nCells+2);
    p.resize(nCells+2);
    for (int i = 0; i != rho.size(); i++){
        double x = x0 + (i-0.5)*dx;
        rho[i] = ((x < 0.5) ? 1 : 0.5);
        v[i] = ((x < 0.5) ? 0.5 : 0.5);
        p[i] = ((x < 0.5) ? 1 : 1);
    }
    double T = 0.25;
    mat u0;
    u0.push_back(rho);
    u0.push_back(v);
    u0.push_back(p);
    solver(u0, T, x0, x1, nCells, gamma, "Godunov");
    return 0;
}