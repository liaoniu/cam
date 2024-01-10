#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
using namespace std;


int main()
{
    double a = 1.0;      // Advection coefficient
    int nPoints = 100;
    double dx = 0.01;
    double dt = dx;   // Time step size
    double T = 1.0;      // Total simulation time
    int nx = nPoints + 2;  // Number of spatial points
    double x1 = 1.0, x0 = 0.0, t0 = 0.0;

    vector<double> u(nPoints+2,0.0);
    
    //double dx = (x1 - x0)/nPoints;

    for (int i = 0; i != u.size(); i++){
        double x = x0 + (i-1)*dx;
        u[i] = sin(x);
    }
    ofstream initial("initial.dat");
    for (int i = 1; i != nPoints+1; i++){
        double x = x0 + (i-1)*dx;
        initial << x << " " << u[i] << endl;
    }
    initial.close();


    vector<double> u_new(nPoints+2, 0.0);
    double t = t0;
    do{
        t += dt;
        u[0] = u[nPoints];
        u[nPoints+1] = u[1];
        for (int i = 1; i != nPoints+1; i++){
            double du_dx = (u[i] - u[(i-1)]) / dx;
            u_new[i] = u[i] - a * du_dx * dt;
        }
        u = u_new;
    }while(t < T);
    ofstream output("advectionsResults.dat");
    for (int i = 1; i != nPoints+1; i++){
        double x = x0 + (i-1)*dx;
        output << x << " " << u[i] << endl;
    
    }
    output.close();
    return 0;
}