#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;
using mat = vector<vector<double>>;

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

tuple<vector<double>, vector<double>, vector<double>> primitiveToConserved(vector<double> rho, vector<double> v, vector<double> p, double gamma){
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


tuple<vector<double>, vector<double>, vector<double>> conservedToPrimitive(vector<double> rho, vector<double> rhov, vector<double> E, double gamma){
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



tuple<double, double, double> rhovpToV(double rho, double v, double p, double cs){
    double v0 = rho - p/cs/cs;
    double v1 = v + p/rho/cs;
    double v2 = v - p/rho/cs;
    return tie(v0, v1, v2);
}

tuple<double, double, double> VToRhovp(double v0, double v1, double v2, double cs){
    double rho, v, p;
    v = (v1+v2)/2;
    rho = v0/(1-(v1-v)/cs);
    p = (v1 - v)*rho*cs;
    return tie(rho, v, p);
}


tuple<vector<double>, vector<double>, vector<double>> primitiveToCharacteristic(vector<double> rho, vector<double> v, vector<double> p, double cs){
    vector<double> v0_vec, v1_vec, v2_vec;
    for (int i = 0; i != rho.size(); i++){
        double v0, v1, v2;
        tie(v0, v1, v2) = rhovpToV(rho[i], v[i], p[i], cs);
        v0_vec.push_back(v0);
        v1_vec.push_back(v1);
        v2_vec.push_back(v2);
    }
    return tie(v0_vec, v1_vec, v2_vec);
}

tuple<vector<double>, vector<double>, vector<double>> characteristicToPrimitive(vector<double> v0_vec, vector<double> v1_vec, vector<double> v2_vec, double cs){
    vector<double> rho_vec, v_vec, p_vec;
    for (int i = 0; i != v0_vec.size(); i++){
        double rho, v, p;
        tie(rho, v, p) = VToRhovp(v0_vec[i], v1_vec[i], v2_vec[i],cs);
        rho_vec.push_back(rho);
        v_vec.push_back(v);
        p_vec.push_back(p);
    }
    return tie(rho_vec, v_vec, p_vec);
}

vector<double> find_cs(vector<double> p, vector<double> rho, double gamma){
    vector<double> cs;
    for (int i = 0; i != p.size(); i++){
        cs.push_back(sqrt(gamma*p[i]/rho[i]));
    }
    return cs;
}

double computeTimeStep(vector<double> v, double dx, vector<double> cs){
    double dt;
    double a_max = 0;
    double C = 0.8;
    for (int i = 0; i != v.size(); i++){
        a_max = ((a_max < fabs(v[i]) + cs[i]) ? fabs(v[i]) + cs[i] : a_max);
    }
    
    dt = C*dx/a_max;
    return dt;
}

vector<double> f(vector<double> u, double gamma){
    vector<double> res;
    double rho, v, p, E;
    E = u[2];
    tie(rho, v, p) = indi_conservedToPrimitive(u[0], u[1], u[2], gamma);
    res.push_back(rho*v);
    res.push_back(rho*v*v + p);
    res.push_back((E+p)*v);
    return res;
}



vector<double> getFlux(vector<double> uL_vec, vector<double> uR_vec, double dt, double dx, string method, double gamma){
    vector<double> flux;

    if (method == "FORCE"){
        vector<double> LF, RI, mid_u;
        for (int i = 0; i != uL_vec.size(); i++){
            double uL = uL_vec[i], uR = uR_vec[i];     
            LF.push_back(0.5*dx/dt*(uL - uR) + 0.5*(f(uR_vec, gamma)[i] + f(uL_vec, gamma)[i]));
            mid_u.push_back(0.5*(uL+uR) - 0.5*dt/dx*(f(uR_vec,gamma)[i] - f(uL_vec,gamma)[i]));
        }
        RI = f(mid_u, gamma);
        for (int i = 0; i != LF.size(); i++){
            flux.push_back(0.5*(LF[i] + RI[i]));
        }
    }
        

    return flux;
}


void solver(mat u0, double T, double x0, double x1, int nCells = 100, double gamma = 1.4){
    double dx = (x1-x0)/nCells;
    double t = 0;
    double dt = 0;

    mat u;
    u = u0;
    mat u_new;
    mat flux;
    u_new.resize(3);
    for(int i = 0; i != u_new.size(); i++){
        u_new[i].resize(nCells+2);
    }
    
    flux.resize(3);
    for(int i = 0; i != flux.size(); i++){
        flux[i].resize(nCells+1);
    }
    vector<double> rho, v, p;
    do{

        tie(rho, v, p) = tie(u[0], u[1], u[2]);
        vector<double> cs = find_cs(p, rho, gamma);
        dt = computeTimeStep(v, dx, cs);
        t += dt;
        tie(u[0], u[1], u[2]) = primitiveToConserved(u[0], u[1], u[2], gamma);

        for (int i = 0; i != u.size(); i++){
            u[i][0] = u[i][1];
            u[i][nCells+1] = u[i][nCells];
        }
        
        for (int j = 0; j != nCells + 1; j++){
            vector<double> ui, uiP;
            for (int i = 0; i != u.size(); i++){
                ui.push_back(u[i][j]);
                uiP.push_back(u[i][j+1]);
            }
            
            vector<double> flux_vec = getFlux(ui, uiP, dt, dx, "FORCE", gamma);
            for (int i = 0; i != flux.size(); i++){
                flux[i][j] = flux_vec[i];
            }
        }

        for (int j = 1; j != nCells + 1; j++){
            for (int i = 0; i != u.size(); i++){
            u_new[i][j] = u[i][j] - dt/dx*(flux[i][j] - flux[i][j-1]);
            }
        }
        u = u_new;
        tie(u[0], u[1], u[2]) = conservedToPrimitive(u[0], u[1], u[2], gamma);
        
        

    }while(t<T);

    ofstream output1("rho.dat");
    ofstream output2("v.dat");
    ofstream output3("p.dat");
    for (int i = 1; i != nCells + 1; i++){
        double x = x0 + (i-0.5)*dx;
        output1 << x << " " << u[0][i] << endl;
        output2 << x << " " << u[1][i] << endl;
        output3 << x << " " << u[2][i] << endl;
    }
    output1.close();
    output2.close();
    output3.close();



}


int main()
{
    double x0 = 0, x1 = 1;
    int nCells = 100;
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