#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;

tuple<double, double, double> rhovpToOriginal(double rho, double v, double p, double gamma){
    double rhoV = rho*v;
    double epsilon = p/(gamma - 1)/rho;
    double E = rho*epsilon + rho*v*v/2;
    return tie(rho, rhoV, E);
}

tuple<double, double, double> rhovpToV(double rho, double v, double p, double cs){
    double v0 = rho - p/cs/cs;
    double v1 = v + p/rho/cs;
    double v2 = v - p/rho/cs;
    return tie(v0, v1, v2);
}

tuple<vector<double>, vector<double>, vector<double>> rhoToV_vec(vector<double> rho, vector<double> v, vector<double> p, double cs){
    vector<double> v0_vec, v1_vec, v2_vec;
    for (int i = 0; i != rho.size(); i++){
        double v0, v1, v2;
        tie(v0, v1, v2) = rhovpToV(rho[i], v[i], p[i], cs);
    }
    return tie(v0_vec, v1_vec, v2_vec);
}


double computeTimeStep(vector<double> u, double dx){
    double dt;
    double a_max = 0;
    double C = 0.8;
    for (int i = 0; i != u.size(); i++){
        a_max = ((a_max < fabs(u[i])) ? fabs(u[i]) : a_max);
    }
    dt = C*dx/a_max;
    return dt;
}


double getFlux(double uL, double uR, double dt, double dx, string method, double (*f)(double, double, double), double a, double cs){
    double flux;
    if (method == "Lax_Friedrichs"){
        flux = 0.5*dx/dt*(uL - uR) + 0.5*(f(a,uR,cs) + f(a,uL,cs));
    }
    else if (method == "FORCE"){
        double LF = 0.5*dx/dt*(uL - uR) + 0.5*(f(a,uR,cs) + f(a,uL,cs));
        double mid_u = 0.5*(uL+uR) - 0.5*dt/dx*(f(a,uR,cs) - f(a,uL,cs));
        double RI = mid_u*mid_u/2;
        flux = 0.5*(LF+RI);
    }
    else if (method == "Godunov"){
        double s = 0.5*(uL+uR);
        double mid_u;
        if (uL > uR){
            if (s>0)
                mid_u = uL;
            else
                mid_u = uR;
        }
        else {
            if (uL > 0)
                mid_u = uL;
            else if (uR < 0)
                mid_u = uR;
            else
                mid_u = 0;
        }
        flux = mid_u*mid_u/2;
    }

    return flux;
}


double f1(double a, double u, ...){
    return a*u;
}

double f2(double a, double u, double cs){
    return a*u + cs;
}

double f3(double a, double u, double cs){
    return a*u - cs;
}


vector<double> solver(vector<double> u0, double T, double x0, double x1, string method, int nCells = 100, double (*f)(double, double, double)){
    double dx = (x1-x0)/nCells;
    double t = 0;
    double dt = 0;
    vector<double> u;
    u.resize(nCells+2);
    u = u0;
    vector<double> flux;
    flux.resize(nCells+1);
    vector<double> u_new;
    u_new.resize(nCells+2);
    int index = 0;
    do{

        dt = computeTimeStep(u, dx);
        t += dt;

        u[0] = u[1];
        u[nCells+1] = u[nCells];
        for (int i = 0; i != nCells + 1; i++){
            flux[i] = getFlux(u[i], u[i+1], dt, dx, method, f1, a, cs);
        }
        for (int i = 1; i != nCells + 1; i++){
            u_new[i] = u[i] - dt/dx*(flux[i] - flux[i-1]);
        }
        u = u_new;
        

    }while(t<T);


    return u;
}


int main()
{
    double x0 = 0, x1 = 1;
    int nCells = 10;
    double dx = (x1 - x0)/nCells;
    double gamma = 1.4;
    vector<double> rho, v, p;
    rho.resize(nCells+2);
    v.resize(nCells+2);
    p.resize(nCells+2);
    for (int i = 0; i != rho.size(); i++){
        double x = x0 - (i-0.5)*dx;
        rho[i] = ((x < 0.5) ? 1 : 0.125);
        v[i] = ((x < 0.5) ? 0 : 0);
        p[i] = ((x < 0.5) ? 1 : 0.1);
    }
    double T = 0.25;
    



}

