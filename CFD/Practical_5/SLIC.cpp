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

double zeta(double delta_minus, double delta_plus, string limiter = "Minbee"){
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

tuple<vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>> dataReconstruct(vector<double> rho, vector<double> rhov, vector<double> E){
    vector<double> rhoL, rhoR, rhovL, rhovR, EL, ER;
    double w = 0;
    int nCells = rho.size() - 2;
    rhoL.resize(nCells+2);
    rhoR.resize(nCells+2);
    rhovL.resize(nCells+2);
    rhovR.resize(nCells+2);
    EL.resize(nCells+2);
    ER.resize(nCells+2);
    for (int i = 1; i != nCells+1; i++){
        double delta_rho, delta_rhov, delta_E;
        double delta_minus, delta_plus;
        //Find rho_bar
        delta_minus = rho[i] - rho[i-1];
        delta_plus = rho[i+1] - rho[i];
        delta_rho = 0.5*(1+w)*delta_minus + 0.5*(1-w)*delta_plus;
        rhoL[i] = rho[i] - 0.5*zeta(delta_minus, delta_plus)*delta_rho;
        rhoR[i] = rho[i] + 0.5*zeta(delta_minus, delta_plus)*delta_rho;
        //Find rhov_bar
        delta_minus = rhov[i] - rhov[i-1];
        delta_plus = rhov[i+1] - rhov[i];
        delta_rhov = 0.5*(1+w)*delta_minus + 0.5*(1-w)*delta_plus;
        rhovL[i] = rhov[i] - 0.5*zeta(delta_minus, delta_plus)*delta_rhov;
        rhovR[i] = rhov[i] + 0.5*zeta(delta_minus, delta_plus)*delta_rhov;
        //Find E_bar
        delta_minus = E[i] - E[i-1];
        delta_plus = E[i+1] - E[i];
        delta_E = 0.5*(1+w)*delta_minus + 0.5*(1-w)*delta_plus;
        EL[i] = E[i] - 0.5*zeta(delta_minus, delta_plus)*delta_E;
        ER[i] = E[i] + 0.5*zeta(delta_minus, delta_plus)*delta_E;
    }
    rhoL[0] = rho[0];
    rhoR[0] = rho[0];
    rhovL[0] = rhov[0];
    rhovR[0] = rhov[0];
    EL[0] = E[0];
    ER[0] = E[0];
    rhoL[nCells+1] = rho[nCells+1];
    rhoR[nCells+1] = rho[nCells+1];
    rhovL[nCells+1] = rhov[nCells+1];
    rhovR[nCells+1] = rhov[nCells+1];
    EL[nCells+1] = E[nCells+1];
    ER[nCells+1] = E[nCells+1];
    return tie(rhoL, rhoR, rhovL, rhovR, EL, ER);
}

tuple<vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>> halfTimeUpdate(vector<double> rhoL, vector<double> rhovL, vector<double> EL, vector<double> rhoR, vector<double> rhovR, vector<double> ER, double gamma, double dx, double dt){
    vector<double> rhoL_HT, rhoR_HT, rhovL_HT, rhovR_HT, EL_HT, ER_HT;
    int nCells = rhoL.size() - 2;
    rhoL_HT.resize(nCells+2);
    rhoR_HT.resize(nCells+2);
    rhovL_HT.resize(nCells+2);
    rhovR_HT.resize(nCells+2);
    EL_HT.resize(nCells+2);
    ER_HT.resize(nCells+2);
    for (int i = 0; i != rhoL_HT.size(); i++){
        double update_rho_R, update_rhov_R, update_E_R, update_rho_L, update_rhov_L, update_E_L, update_rho, update_rhov, update_E;
        tie(update_rho_R, update_rhov_R, update_E_R) = f(rhoR[i], rhovR[i], ER[i], gamma);
        tie(update_rho_L, update_rhov_L, update_E_L) = f(rhoL[i], rhovL[i], EL[i], gamma);
        update_rho = 0.5*dt/dx*(update_rho_R - update_rho_L);
        update_rhov = 0.5*dt/dx*(update_rhov_R - update_rhov_L);
        update_E = 0.5*dt/dx*(update_E_R - update_E_L);
        rhoL_HT[i] = rhoL[i] - update_rho;
        rhoR_HT[i] = rhoR[i] - update_rho;
        rhovL_HT[i] = rhovL[i] - update_rhov;
        rhovR_HT[i] = rhovR[i] - update_rhov;
        EL_HT[i] = EL[i] - update_E;
        ER_HT[i] = ER[i] - update_E;
    }
    return tie(rhoL_HT, rhoR_HT, rhovL_HT, rhovR_HT, EL_HT, ER_HT);
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
    if (method == "FORCE"){
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




    return tie(flux[0], flux[1], flux[2]);
}



void solver(mat u0, double T, double x0, double x1, int nCells = 100, double gamma = 1.4, string method = "SLIC"){
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

    do{
        dt = calc_dt(rho, rhov, E, gamma, dx);
        t += dt;
        rho[0] = rho[1];
        rho[nCells+1] = rho[nCells];
        rhov[0] = rhov[1];
        rhov[nCells+1] = rhov[nCells];
        E[0] = E[1];
        E[nCells+1] = E[nCells];
        if (method != "SLIC"){
            for (int i = 0; i != flux_1.size(); i++){
                tie(flux_1[i], flux_2[i], flux_3[i]) = getFlux(rho[i], rho[i+1], rhov[i], rhov[i+1], E[i], E[i+1], gamma, dx, dt, method);
            }
        }
        else{
            vector<double> rhoL, rhoR, rhovL, rhovR, EL, ER;
            tie(rhoL, rhoR, rhovL, rhovR, EL, ER) = dataReconstruct(rho, rhov, E);
            tie(rhoL, rhoR, rhovL, rhovR, EL, ER) = halfTimeUpdate(rhoL, rhovL, EL, rhoR, rhovR, ER, gamma, dx, dt);
            for (int i = 0; i != flux_1.size(); i++){
                tie(flux_1[i], flux_2[i], flux_3[i]) = getFlux(rhoR[i], rhoL[i+1], rhovR[i], rhovL[i+1], ER[i], EL[i+1], gamma, dx, dt, "FORCE");
            }
        }



        for (int i = 1; i != nCells+1; i++){
            rho_new[i] = rho[i] - dt/dx*(flux_1[i] - flux_1[i-1]);
            rhov_new[i] = rhov[i] - dt/dx*(flux_2[i] - flux_2[i-1]);
            E_new[i] = E[i] - dt/dx*(flux_3[i] - flux_3[i-1]);
        }
        rho = rho_new;
        rhov = rhov_new;
        E = E_new;

    }while(t<T);

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
    int nCells = 1000;
    double dx = (x1 - x0)/nCells;
    double gamma = 1.4;
    vector<double> rho, v, p;
    rho.resize(nCells+2);
    v.resize(nCells+2);
    p.resize(nCells+2);
    for (int i = 0; i != rho.size(); i++){
        double x = x0 + (i-0.5)*dx;
        rho[i] = ((x < 0.5) ? 1 : 0.125);
        v[i] = ((x < 0.5) ? 0 : 0);
        p[i] = ((x < 0.5) ? 1 : 0.1);
    }
    double T = 0.25;
    mat u0;
    u0.push_back(rho);
    u0.push_back(v);
    u0.push_back(p);
    solver(u0, T, x0, x1, nCells, gamma, "SLIC");
    return 0;
}


