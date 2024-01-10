#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;


int main() {
    // Define parameters
    double a = 1.0;      // Advection coefficient
    double dx = 0.01;    // Spatial step size
    double dt = 0.01;   // Time step size
    double T = 1.0;      // Total simulation time
    int nx = static_cast<int>(1.0 / dx) + 2;  // Number of spatial points
    int nt = static_cast<int>(T / dt) + 1;    // Number of time steps

    // Create a vector for the solution
    std::vector<double> u(nx, 0.0);

    // Set initial condition
    for (int i = 0; i < nx; i++) {
        u[i] = sin(i * dx);
    }

    // Time-stepping loop
    for (int n = 0; n < nt - 1; n++) {
        std::vector<double> u_new(nx, 0.0);
        for (int i = 0; i < nx - 1; i++) {
            // Compute the spatial derivative using forward difference
            double du_dx = (u[i + 1] - u[i]) / dx;

            // Update the solution using the advection equation
            u_new[i] = u[i] - a * du_dx * dt;
        }
        u_new[nx - 1] = u[nx - 1];  // Maintain the boundary condition

        u = u_new;
    }

    ofstream output("gpt2.dat");
    for (int i = 1; i != nx; i++){
        double x = 0 + (i-1)*dx;
        output << x << " " << u[i] << endl;
    }

    return 0;
}
