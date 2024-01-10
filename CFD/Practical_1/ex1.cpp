#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
using namespace std;

class Solver{
public:
    double x0 = 0, x1 = 1, t0 = 0, T = 1, a = 1;
    int nPoints = 100;
    void forward(double dt, double a=1.0, double t0=0.0, double T=1.0){
        
        vector<double> u;
        u.resize(nPoints+2);
        double dx = (x1 - x0)/nPoints;
        dt = dx;
        for (int i = 0; i != u.size(); i++){
            double x = x0 + (i-1)*dx;
            u[i] = sin(x);
        }
        int nx = nPoints + 2;  // Number of spatial points
        ofstream initial("initial.dat");
        for (int i = 1; i != nPoints+1; i++){
            double x = x0 + (i-1)*dx;
            initial << x << " " << u[i] << endl;
        }


        vector<double> u_new(nPoints+2, 0.0);
        double t = t0;
        do{
            t += dt;
            u[0] = u[nPoints];
            u[nPoints+1] = u[1];
            for (int i = 1; i != nPoints+1; i++){
                double du_dx = (u[i+1] - u[((i) + nx) % nx]) / dx;
                u_new[i] = u[i] - a * du_dx * dt;
            }
            u = u_new;
        }while(t < T);
        ofstream output("caonimab.dat");
        for (int i = 1; i != nPoints+1; i++){
            double x = x0 + (i-1)*dx;
            output << x << " " << u[i] << endl;
        }
        output.close();
    }

    void backward(double dt, double a=1.0, double t0=0.0, double T=1.0){
        
        vector<double> u;
        u.resize(nPoints+2);
        double dx = (x1 - x0)/nPoints;
        dt = dx;
        for (int i = 0; i != u.size(); i++){
            double x = x0 + (i-1)*dx;
            u[i] = sin(x);
        }

        ofstream initial("initial.dat");
        for (int i = 1; i != nPoints+1; i++){
            double x = x0 + (i-1)*dx;
            initial << x << " " << u[i] << endl;
        }


        vector<double> u_new(nPoints+2, 0);
        double t = t0;
        do{
            t += dt;
            u[0] = u[nPoints];
            u[nPoints+1] = u[1];
            for (int i = 1; i != nPoints+1; i++){
                u_new[i] = u[i] - a*(dt/dx)*(u[i] - u[i-1]);
            }
            u = u_new;
        }while(t < T);
        ofstream output("advectionsResults.dat");
        for (int i = 1; i != nPoints+1; i++){
            double x = x0 + (i-1)*dx;
            output << x << " " << u[i] << endl;
        }
        output.close();
    }

    void central(double dt, double a=1.0, double t0=0.0, double T=1.0){
        vector<double> u;
        u.resize(nPoints+2);
        double dx = (x1 - x0)/nPoints;
        dt = dx;
        for (int i = 0; i != u.size(); i++){
            double x = x0 + (i-1)*dx;
            u[i] = sin(x);
        }

        ofstream initial("initial.dat");
        for (int i = 1; i != nPoints+1; i++){
            double x = x0 + (i-1)*dx;
            initial << x << " " << u[i] << endl;
        }

        vector<double> u_new(nPoints+2, 0.0);

        double t = t0;
        do{
            t += dt;
            u[0] = u[nPoints];
            u[nPoints+1] = u[1];
            for (int i = 1; i != nPoints+1; i++){
                u_new[i] = u[i] - a*(dt/dx/2)*(u[i+1] - u[i-1]);
            }
            u = u_new;

        }while(t<T);
        ofstream output("advectionsResults.dat");
        for (int i = 1; i != nPoints+1; i++){
            double x = x0 + (i-1)*dx;
            output << x << " " << u[i] << endl;
        }
        output.close();
    }

};




int main()
{
    Solver s;
    s.backward(0.001);
    return 0;
}

