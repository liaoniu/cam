#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
using namespace std;

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


double getFlux(double uL, double uR, double dt, double dx, string method){
    double flux;
    if (method == "Lax_Friedrichs"){
        flux = 0.5*dx/dt*(uL - uR) + 0.5*(uR*uR/2 + uL*uL/2);
    }
    else if (method == "FORCE"){
        double LF = 0.5*dx/dt*(uL - uR) + 0.5*(uR*uR/2 + uL*uL/2);
        double mid_u = 0.5*(uL+uR) - 0.5*dt/dx*(uR*uR/2 - uL*uL/2);
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


void solver(vector<double> u0, double T, double x0, double x1, string method, int nCells = 100){
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
        if (t > index*T/5 || index == 0){
            
            string str = "middle_";
            stringstream num;
            num << index;
            str += num.str();
            str += ".dat";
            index += 1;
            ofstream output(str);
            for (int i = 1; i != nCells+1; i++){
                double x = x0 + (i-0.5)*dx;
                output << x << " " << u[i] << endl;
            }
            output.close();
        }
        dt = computeTimeStep(u, dx);
        t += dt;

        u[0] = u[1];
        u[nCells+1] = u[nCells];
        for (int i = 0; i != nCells + 1; i++){
            flux[i] = getFlux(u[i], u[i+1], dt, dx, method);
        }
        for (int i = 1; i != nCells + 1; i++){
            u_new[i] = u[i] - dt/dx*(flux[i] - flux[i-1]);
        }
        u = u_new;
        

    }while(t<T);

    ofstream output("final.dat");
    for (int i = 1; i != nCells + 1; i++){
        double x = x0 + (i-0.5)*dx;
        output << x << " " << u[i] << endl;
    }
    output.close();

}



int main()
{
    double x0 = 0, x1 = 1;
    int nCells = 1000;
    double dx = (x1 - x0)/nCells;
    vector<double> u0;
    u0.resize(nCells+2);
    for (int i = 0; i != u0.size(); i++){
        double x = x0 + (i-0.5)*dx;
        u0[i] = ((x < 0.5) ? 2 : 1);
    }



    solver(u0, 0.2, x0, x1, "Godunov", nCells);
    return 0;
}



