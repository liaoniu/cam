#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;
// calculating timestep
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
// Define the initial condition u0(x)
double u0(double x) {
    if (x < 0.5)
        return 2.0;
    else
        return 1.0;
}

int main() {
    int nPoints = 100; // Number of spatial points
    double x0 = 0; // Initial spatial coordinate
    double x1 = 1.0; // Final spatial coordinate
    double tStart = 0.0; // Initial time
    double tStop = 0.2; // Final time
    double dx = (x1 - x0) / nPoints; // Spatial grid spacing
    double dt = 0.8 * dx; // Time step size
    double c = dt / dx; // Courant number
    

    vector<double> u(nPoints + 2); // Array to store the data at each spatial point
    vector<double> uPlus1(nPoints + 2); // Temporary array for the updated data
    vector<double> flux(nPoints + 1); // Array to store flux values

    // Initialize the data based on the initial condition u0(x)
    for (int i = 0; i < u.size(); i++) {
        double x = x0 + (i - 0.5) * dx;
        u[i] = u0(x);
    }
    dt = computeTimeStep(u, dx);
    tStart += dt;
    do{
           // Apply periodic boundary conditions
        u[0] = u[1];
        u[nPoints+1] = u[nPoints];
        for (int i = 0; i != nPoints + 1; i++){
            flux[i] = 0.5*dx/dt*(u[i] - u[i+1]) + 0.5*(u[i+1]*u[i+1]/2 + u[i]*u[i]/2);
        }
        for (int i = 1; i != nPoints + 1; i++){
            uPlus1[i] = u[i] - dt/dx*(flux[i] - flux[i-1]);
        }
        // Update the data for the next time step
        u = uPlus1;
        dt = computeTimeStep(u, dx);
        tStart += dt;
         }while(tStart<tStop);
    


    // Output the results to a file named "ShockWaveResults.dat"
    ofstream outfile("ShockWaveResults.dat");
    if (!outfile) return -1; // Check if the file is opened successfully

    for (int i = 1; i <= nPoints; i++) {
        double x = x0 + (i - 0.5) * dx;
        outfile << x << " " << u[i] << endl;
    }

    outfile.close(); // Close the output file

    return 0; // Exit the program
}