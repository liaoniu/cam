#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;
using mat = vector<array<double, 3>>;   // u[i][var] = var[i]


array<double, 3> indi_primitiveToConserved(array<double, 3> &u_pri, double &gamma);
array<double, 3> indi_conservedToPrimitive(array<double, 3> &u_con, double &gamma);
mat conservedToPrimitive(mat &u, double &gamma);
mat primitiveToConserved(mat &u, double &gamma);
double calc_dt(mat &u_con, double &gamma, double &dx, double &dy);
array<double, 3> f(array<double, 3> &u, double &gamma);
double zeta(double delta_minus, double delta_plus, string limiter);
tuple<mat, mat> dataReconstruct(mat &u, string direction);
array<double, 3> getXFlux(array<double, 3> &uL, array<double, 3> &uR); //u_left and u_right
mat split(mat &u, mat &flux, double dt, double &dx, int &nxCells, double &gamma, string direction, string method);
bool inCircle(double x, double y, double R);
void transmissiveBC(mat& u, int nxCells);


bool inCircle(double x, double y, double R){
    return ((x-1)*(x-1) + (y-1)*(y-1) < R) ? true : false;
}

void transmissiveBC(mat &u, int nxCells){
    for (int var = 0; var != 3; var++){
        u[0][var] = u[1][var];
        u[nxCells+1][var] = u[nxCells][var];
    }
}

array<double, 3> indi_primitiveToConserved(array<double, 3> &u_pri, double &gamma){
    array<double, 3> u_con = u_pri;
    double rho, v, p;
    tie(rho, v, p) = tie(u_pri[0], u_pri[1], u_pri[2]);
    double rhoV = rho*v;
    double epsilon = p/(gamma - 1)/rho;
    double E = rho*epsilon + rho*(v*v)/2;
    tie(u_con[0], u_con[1], u_con[2]) = tie(rho, rhoV, E);
    return u_con;
}

array<double, 3> indi_conservedToPrimitive(array<double, 3> &u_con, double &gamma){
    array<double, 3> u_pri = u_con;
    double rho, rhoV, E;
    tie(rho, rhoV, E) = tie(u_con[0], u_con[1], u_con[2]);
    double v = rhoV/rho;
    double epsilon = E/rho - 0.5*(v*v);
    double p = (gamma-1)*rho*epsilon;
    tie(u_pri[0], u_pri[1], u_pri[2]) = tie(rho, v, p);
    return u_pri;
}

mat conservedToPrimitive(mat &u, double &gamma){
    //u would be changed
    for (int i = 0; i != u.size(); i++){
        u[i] = indi_conservedToPrimitive(u[i], gamma);
    }
    return u;
}

mat primitiveToConserved(mat &u, double &gamma){
    //u would be changed

    for (int i = 0; i != u.size(); i++){
        u[i] = indi_primitiveToConserved(u[i], gamma);
    }
    
    return u;
}

double calc_dt(mat &u_con, double &gamma, double &dx){
    double delta = dx;
    double C = 0.8;
    double a_max = 0;
    int m = u_con.size();
    for (int i = 0; i != m; i++){
            double rho, v, p;
            array<double, 3> u_pri_indi = indi_conservedToPrimitive(u_con[i], gamma);
            tie(rho, v, p) = tie(u_pri_indi[0], u_pri_indi[1], u_pri_indi[2]);
            double a = sqrt(v*v) + sqrt(gamma*p/rho);
            a_max = ((a_max < a) ? a : a_max);
    }
    double dt = C*delta/a_max;
    return dt;
}

array<double, 3> f(array<double, 3> &u, double &gamma){
    double rho, v, p;
    double rhoV, E;
    tie(rho, rhoV, E) = tie(u[0], u[1], u[2]);
    v = rhoV/rho;
    p = (gamma - 1)*(E - 0.5*(rhoV*v));
    array<double, 3> f_res;
    f_res[0] = rhoV;
    f_res[1] = rhoV*v + p;
    f_res[2] = (E + p)*v;
    return f_res;
}


double zeta(double delta_minus, double delta_plus, string limiter = "Minbee"){
    if (limiter == "Minbee"){
        if (delta_plus == 0){
            return 0;
        }
        else{
            double r = delta_minus/delta_plus;
            return ((r <= 0) ? 0 : ((r <= 1) ? r : ((1 < 2/(1+r)) ? 1 : 2/(1+r))));
        }
    }
    return 0;
}

tuple<mat, mat> dataReconstruct(mat &u){
    mat uL = u, uR = u;
    int nxCells = u.size()-2;
    double w = 0;
    for (int i = 1; i != nxCells+1; i++){
        for (int var = 0; var != 3; var++){
            double delta_minus = u[i][var] - u[i-1][var];
            double delta_plus = u[i+1][var] - u[i][var];
            double delta = 0.5*(1+w)*delta_minus + 0.5*(1-w)*delta_plus;
            uL[i][var] = u[i][var] - 0.5*zeta(delta_minus, delta_plus)*delta;
            uR[i][var] = u[i][var] + 0.5*zeta(delta_minus, delta_plus)*delta;
        }
    }

    return tie(uL, uR);
}

tuple<mat, mat> halfTimeUpdate(mat &uL, mat &uR, double gamma, double dx, double dt){
    int nxCells = uL.size()-2;
    for (int i = 0; i != nxCells+2; i++){
        array<double, 3> fL = f(uL[i], gamma);
        array<double, 3> fR = f(uR[i], gamma);
        for (int var = 0; var != 3; var++){
            uL[i][var] = uL[i][var] - 0.5*dt/dx*(fR[var] - fL[var]);
            uR[i][var] = uR[i][var] - 0.5*dt/dx*(fR[var] - fL[var]);
        }
        
    }
    return tie(uL, uR);
}




array<double, 3> getXFlux(array<double, 3> &uL, array<double, 3> &uR, double &gamma, double &dx, double &dt, string method){ //u_left and u_right
    array<double, 3> flux, u_mid;
    if (method == "FORCE"){
        array<double, 3> fL, fR, LF, RI;
        fL = f(uL, gamma);
        fR = f(uR, gamma);
        for (int var = 0; var != 3; var++){
            LF[var] = 0.5*dx/dt*(uL[var] - uR[var]) + 0.5*(fL[var] + fR[var]);
            u_mid[var] = 0.5*(uL[var] + uR[var]) - 0.5*dt/dx*(fR[var] - fL[var]);
        }
        RI = f(u_mid, gamma);
        for (int var = 0; var != 3; var++){
            flux[var] = 0.5*(LF[var] + RI[var]);
        }

    }
    return flux;
}



mat split(mat &u, mat &flux, double dt, double &dx, int &nxCells, double &gamma, string method){
    mat uL, uR;
    if (method == "SLIC"){
        tie(uL, uR) = dataReconstruct(u);
        halfTimeUpdate(uL, uR, gamma, dx, dt);
        //getXFLux
        for (int i = 0; i != nxCells+1; i++){
                flux[i] = getXFlux(uR[i], uL[i+1], gamma, dx, dt, "FORCE");
        }
        for (int i = 1; i != nxCells+1; i++){
            for (int k = 0; k != 3; k++){
                u[i][k] = u[i][k] - dt/dx*(flux[i][k] - flux[i-1][k]);
            }
        }
    }
    else {
        //getXFLux
        for (int i = 0; i != nxCells+1; i++){
                flux[i] = getXFlux(u[i], u[i+1], gamma, dx, dt, method);
        }
        for (int i = 1; i != nxCells+1; i++){
            for (int k = 0; k != 3; k++){
                u[i][k] = u[i][k] - dt/dx*(flux[i][k] - flux[i-1][k]);
            }

        }
    }

    //Need ensure boundary condition here
    transmissiveBC(u, nxCells);
    return u;
}



void solver(mat &u0, double T, double x0, double x1, int nxCells = 100, double gamma = 1.4, string method = "FORCE"){
    double dx = (x1-x0)/nxCells;
    double t = 0;
    double dt = 0;
    mat flux;
    mat u = u0;
    flux.resize(nxCells+1);
    u = primitiveToConserved(u, gamma);

    do{
        dt = calc_dt(u, gamma, dx);
        t += dt;
        //Transmissive BC
        transmissiveBC(u, nxCells);
        //Dimensional Split 1
        mat u_temp = u;
        // split(u_temp, flux, dt, dx, dy, nxCells, nyCells, gamma, "x", method);
        // split(u_temp, flux, dt, dx, dy, nxCells, nyCells, gamma, "y", method);
        split(u, flux, dt, dx, nxCells, gamma, method);
        // for (int i = 0; i != nxCells + 2; i++){
        //     for (int j = 0; j != nyCells + 2; j++){
        //         for (int var = 0; var != 4; var++){
        //             u[i][j][var] = 0.5*(u[i][j][var] + u_temp[i][j][var]);
        //         }
        //     }
        // }

        //Dimensionial Split 2
        // split(u, flux, dt/2, dx, dy, nxCells, nyCells, gamma, "x", method);
        // split(u, flux, dt, dx, dy, nxCells, nyCells, gamma, "y", method);
        // split(u, flux, dt/2, dx, dy, nxCells, nyCells, gamma, "x", method);

        //cout << t << endl;
    }while(t<T);

    u = conservedToPrimitive(u, gamma);
    ofstream output1("rho.dat");
    ofstream output2("v.dat");
    ofstream output3("p.dat");


    for (int i = 1; i != nxCells + 1; i++){
        double x = x0 + (i-0.5)*dx;
        output1 << x << " " << u[i][0] << endl;
        output2 << x << " " << u[i][1] << endl;
        output3 << x << " " << u[i][2] << endl;
        }

    output1.close();
    output2.close();
    output3.close();

}



int main()
{
    double x0 = 0, x1 = 1, T = 0.25, gamma = 1.4, nxCells = 100, nyCells = 100;
    double dx = (x1 - x0)/nxCells;
    mat u0;
    u0.resize(nxCells+2);
    for (int i = 0; i != nxCells+2; i++){
        double x = x0 + (i-0.5)*dx;
        u0[i][0] = ((x < 0.5) ? 1 : 0.125);
        u0[i][1] = ((x < 0.5) ? 0 : 0);
        u0[i][2] = ((x < 0.5) ? 1 : 0.1);
    }

    solver(u0, T, x0, x1,nxCells, gamma, "SLIC");
    return 0;

}