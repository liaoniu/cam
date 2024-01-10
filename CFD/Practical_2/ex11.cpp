#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
using namespace std;

void solver(double nPoints, double a, double T, string method, vector<double> u0){
    
    double x0 = 0, x1 = 2, t0 = 0;
    double dx = (x1-x0)/nPoints;
    double C = 0.8;
    double dt = C*dx/a;
    double t = t0;
    vector<double> u;
    u.resize(nPoints + 2);
    
    vector<double> u_new;
    u_new.resize(nPoints+2);
    if (method == "Lax_Friedrichs"){
        u = u0;
        do{
            t += dt;
            u[0] = u[nPoints];
            u[1] = u[nPoints+1];
            for (int i=1; i != nPoints+1; i++){
                u_new[i] = 0.5*(1+C)*u[i-1] + 0.5*(1-C)*u[i+1];
            }
            u = u_new;
        }while(t < T);
    }
    else if (method == "Lax_Wendroff"){
        u = u0;
        do{
            t += dt;
            u[0] = u[nPoints];
            u[1] = u[nPoints+1];
            for (int i=1; i != nPoints+1; i++){
                u_new[i] = 0.5*C*(1+C)*u[i-1] + (1-C*C)*u[i] - 0.5*C*(1-C)*u[i+1];
            }
            u = u_new;
        }while(t < T);
    }

    else if (method == "Warming_Beam"){
        u = u0;
        u.resize(nPoints+3);
        u_new.resize(nPoints+3);
        do{
            t += dt;
            u[0] = u[nPoints];
            u[1] = u[nPoints+1];
            u[2] = u[nPoints+2];
            for (int i=2; i != nPoints+1; i++){
                u_new[i] = -0.5*C*(1-C)*u[i-2] + C*(2-C)*u[i-1] + 0.5*(C-1)*(C-2)*u[i];
            }
            u = u_new;
        }while(t < T);

    }





    ofstream output("results.dat");
        for (int i = 1; i != nPoints+1; i++){
            double x = x0 + (i-1)*dx;
            output << x << " " << u[i] << endl;
        }
        output.close();

}

void originalPlot(int nPoints, double x0, double dx){
    vector<double> s;
    s.resize(nPoints);
    for (int i = 0; i != s.size(); i++){
        double x = x0 + i*dx;
        s[i] = ((0.3 <= x && x <= 0.7) ? 1 : 0);
    }
    ofstream original("original.dat");
    for (int i = 0; i != s.size(); i++){
        double x = x0 + i*dx;
        original << x << " " << s[i] << endl;
    }
    original.close();
}




int main()
{
    int nPoints = 100;
    double x0 = 0, x1 = 2;
    double dx = (x1 - x0)/nPoints;
    vector<double> u0;
    u0.resize(nPoints+2);
    for (int i = 0; i != u0.size(); i++){
        double x = x0 + (i-1)*dx;
        u0[i] = ((0.3 <= x && x <= 0.7) ? 1 : 0);
    }

    /*
    u0.resize(nPoints+3);
    for (int i = 0; i != u0.size(); i++){
        double x = x0 + (i-2)*dx;
        u0[i] = exp(-8*x*x);
    }
    */
    solver(nPoints, 1, 1, "Lax_Friedrichs", u0);

    originalPlot(nPoints, x0, dx);

    return 0;
}