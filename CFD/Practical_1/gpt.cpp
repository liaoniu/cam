#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;

int main() {
    // Define parameters
    double a = 1.0;      // Advection coefficient
    int nPoints = 100;
    double dx = 0.01;
    double dt = dx;   // Time step size
    double T = 1.0;      // Total simulation time
    int nx = nPoints + 2;  // Number of spatial points
    
    // Create a vector for the solution
    std::vector<double> u(nx, 0.0);

    // Set initial condition
    for (int i = 0; i < nPoints+1; i++) {
        double x = (i-1)*dx;
        u[i] = sin(x);
    }
    double t = 0;
    // Time-stepping loop



    do{
        t += dt;
        std::vector<double> u_new(nx, 0.0);
        for (int i = 0; i < nx - 1; i++) {
            // Compute the spatial derivative using forward difference
            double du_dx = (u[i + 1] - u[i]) / dx;

            // Update the solution using the advection equation
            u_new[i] = u[i] - a * du_dx * dt;
        }
        u_new[nx - 1] = u[nx - 1];  // Maintain the boundary condition
        u = u_new;
    }while(t<T);




    ofstream output("gpt.dat");
    for (int i = 1; i != nPoints+1; i++){
        double x = 0 + (i-1)*dx;
        output << x << " " << u[i] << endl;
    }
    return 0;
}
