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

tuple<double, double, double> getFlux(double rhoL, double rhoR, double rhovL, double rhovR, double EL, double ER, double gamma, double dx, double dt){
    vector<double> flux, LF, RI, mid_u;
    vector<double> fL, fR;
    flux.resize(3);
    LF.resize(3);
    RI.resize(3);
    fL.resize(3);
    fR.resize(3);
    mid_u.resize(3);
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
    return tie(flux[0], flux[1], flux[2]);
}



void solver(mat u0, double T, double x0, double x1, int nCells = 100, double gamma = 1.4){
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
        for (int i = 0; i != flux_1.size(); i++){
            tie(flux_1[i], flux_2[i], flux_3[i]) = getFlux(rho[i], rho[i+1], rhov[i], rhov[i+1], E[i], E[i+1], gamma, dx, dt);
            
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
    solver(u0, T, x0, x1, nCells, gamma);
    return 0;
}


