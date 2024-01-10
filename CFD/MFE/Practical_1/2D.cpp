#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;
using mat = vector<vector<array<double, 4>>>;   // u[i][j][var] = var[i][j]
using myTuple = tuple<vector<double>, vector<double>, vector<double>, vector<double>>;


array<double, 4> indi_primitiveToConserved(array<double, 4> &u_pri, double &gamma);
array<double, 4> indi_conservedToPrimitive(array<double, 4> &u_con, double &gamma);
mat conservedToPrimitive(mat &u, double &gamma);
mat primitiveToConserved(mat &u, double &gamma);
double calc_dt(mat &u_con, double &gamma, double &dx, double &dy);
array<double, 4> f(array<double, 4> &u, double &gamma);
array<double, 4> g(array<double, 4> &u, double &gamma);
array<double, 4> getXFlux(array<double, 4> &uL, array<double, 4> &uR); //u_left and u_right
array<double, 4> getYFlux(array<double, 4> &uD, array<double, 4> &uU); //u_down and u_up




array<double, 4> indi_primitiveToConserved(array<double, 4> &u_pri, double &gamma){
    array<double, 4> u_con = u_pri;
    double rho, vx, vy, p;
    tie(rho, vx, vy, p) = tie(u_pri[0], u_pri[1], u_pri[2], u_pri[3]);
    double rhoVx = rho*vx;
    double rhoVy = rho*vy;
    double epsilon = p/(gamma - 1)/rho;
    double E = rho*epsilon + rho*(vx*vx + vy*vy)/2;
    tie(u_con[0], u_con[1], u_con[2], u_con[3]) = tie(rho, rhoVx, rhoVy, E);
    return u_con;
}

array<double, 4> indi_conservedToPrimitive(array<double, 4> &u_con, double &gamma){
    array<double, 4> u_pri = u_con;
    double rho, rhoVx, rhoVy, E;
    tie(rho, rhoVx, rhoVy, E) = tie(u_con[0], u_con[1], u_con[2], u_con[3]);
    double vx = rhoVx/rho;
    double vy = rhoVy/rho;
    double epsilon = E/rho - 0.5*(vx*vx + vy*vy);
    double p = (gamma-1)*rho*epsilon;
    tie(u_pri[0], u_pri[1], u_pri[2], u_pri[3]) = tie(rho, vx, vy, p);
    return u_pri;
}

mat conservedToPrimitive(mat &u, double &gamma){
    //u would be changed
    for (int i = 0; i != u.size(); i++){
        for (int j = 0; j != u[0].size(); j++){
            u[i][j] = indi_conservedToPrimitive(u[i][j], gamma);
        }
    }
    return u;
}

mat primitiveToConserved(mat &u, double &gamma){
    //u would be changed
    for (int i = 0; i != u.size(); i++){
        for (int j = 0; j != u[0].size(); j++){
            u[i][j] = indi_primitiveToConserved(u[i][j], gamma);
        }
    }
    return u;
}

double calc_dt(mat &u_con, double &gamma, double &dx, double &dy){
    double delta = ((dx < dy) ? dx : dy);
    double C = 0.8;
    double a_max = 0;
    int m = u_con.size(), n = u_con[0].size();
    for (int i = 0; i != m; i++){
        for (int j = 0; j != n; j++){
            double rho, vx, vy, p;
            array<double, 4> u_pri_indi = indi_conservedToPrimitive(u_con[i][j], gamma);
            tie(rho, vx, vy, p) = tie(u_pri_indi[0], u_pri_indi[1], u_pri_indi[2], u_pri_indi[3]);
            double a = sqrt(vx*vx + vy*vy) + sqrt(gamma*p/rho);
            a_max = ((a_max < a) ? a : a_max);
        }
    }
    return a_max;
}

array<double, 4> f(array<double, 4> &u, double &gamma){
    double rho, vx, vy, p;
    double rhoVx, rhoVy, E;
    tie(rho, rhoVx, rhoVy, E) = tie(u[0], u[1], u[2], u[3]);
    vx = rhoVx/rho;
    vy = rhoVy/rho;
    p = (gamma - 1)*(E - 0.5*rhoVx*vx);
    array<double, 4> f_res;
    f_res[0] = rhoVx;
    f_res[1] = rhoVx*vx + p;
    f_res[2] = rhoVx*vy;
    f_res[3] = (E + p)*vx;
    return f_res;
}

array<double, 4> g(array<double, 4> &u, double &gamma){
    double rho, vx, vy, p;
    double rhoVx, rhoVy, E;
    tie(rho, rhoVx, rhoVy, E) = tie(u[0], u[1], u[2], u[3]);
    vx = rhoVx/rho;
    vy = rhoVy/rho;
    p = (gamma - 1)*(E - 0.5*rhoVx*vx);
    array<double, 4> g_res;
    g_res[0] = rhoVy;
    g_res[1] = rhoVx*vy;
    g_res[2] = rhoVy*vy + p;
    g_res[3] = (E + p)*vy;
    return g_res;
}

array<double, 4> getXFlux(array<double, 4> &uL, array<double, 4> &uR, double &gamma, double &dx, double &dt, string method){ //u_left and u_right
    array<double, 4> flux, u_mid;
    if (method == "FORCE"){
        array<double, 4> fL, fR, LF, RI;
        fL = f(uL, gamma);
        fR = f(uR, gamma);
        for (int var = 0; var != 4; var++){
            LF[var] = 0.5*dx/dt*(uL[var] - uR[var]) + 0.5*(fL[var] + fR[var]);
            u_mid[var] = 0.5*(uL[var] + uR[var]) - 0.5*dt/dx*(fR[var] - fL[var]);
        }
        RI = f(u_mid, gamma);
        for (int var = 0; var != 4; var++){
            flux[var] = 0.5*(LF[var] + RI[var]);
        }

    }
    return flux;
}

array<double, 4> getYFlux(array<double, 4> &uL, array<double, 4> &uR, double &gamma, double &dy, double &dt, string method){ //u_down and u_up
    array<double, 4> flux, u_mid;
    if (method == "FORCE"){
        array<double, 4> fL, fR, LF, RI;
        fL = g(uL, gamma);
        fR = g(uR, gamma);
        for (int var = 0; var != 4; var++){
            LF[var] = 0.5*dy/dt*(uL[var] - uR[var]) + 0.5*(fL[var] + fR[var]);
            u_mid[var] = 0.5*(uL[var] + uR[var]) - 0.5*dt/dy*(fR[var] - fL[var]);
        }
        RI = f(u_mid, gamma);
        for (int var = 0; var != 4; var++){
            flux[var] = 0.5*(LF[var] + RI[var]);
        }
    }
    return flux;
}


void solver(mat &u0, double T, double x0, double x1, double y0, double y1, int nxCells = 100, int nyCells = 100, double gamma = 1.4, string method = "FORCE"){
    double dx = (x1-x0)/nxCells;
    double dy = (y1-y0)/nyCells;
    double t = 0;
    double dt = 0;
    mat flux;
    mat u = u0;
    flux.resize(nxCells+1, vector<array<double, 4>> (nyCells+1));
    u = primitiveToConserved(u, gamma);

    do{
        dt = calc_dt(u, gamma, dx, dy);
        t += dt;
        //Transmissive BC
        for (int i = 0; i != nxCells+2; i++){
            for (int var = 0; var != 4; var++){
                u[i][0][var] = u[i][1][var];
                u[i][nyCells+1][var] = u[i][nyCells][var];
            }
        }
        for (int j = 0; j != nyCells+2; j++){
            for (int var = 0; var != 4; var++){
                u[0][j][var] = u[1][j][var];
                u[nxCells+1][j][var] = u[nxCells][j][var];
            }
        }
        //Dimensional Split
        mat uBar = u;
        //getXFLux
        for (int i = 0; i != nxCells+1; i++){
            for (int j = 0; j != nyCells+1; j++){
                flux[i][j] = getXFlux(u[i][j], u[i+1][j], gamma, dx, dt, method);
            }
        }

        for (int i = 1; i != nxCells+1; i++){
            for (int j = 1; j != nyCells+1; j++){
                for (int k = 0; k != 4; k++){
                    uBar[i][j][k] = u[i][j][k] - dt/dx*(flux[i][j][k] - flux[i-1][j][k]);
                }
            }
        }
        
        //getYFlux
        for (int i = 0; i != nxCells+1; i++){
            for (int j = 0; j != nyCells+1; j++){
                flux[i][j] = getYFlux(uBar[i][j], uBar[i][j+1], gamma, dy, dt, method);
            }
        }

        for (int i = 1; i != nxCells+1; i++){
            for (int j = 1; j != nyCells+1; j++){
                for (int k = 0; k != 4; k++){
                    u[i][j][k] = uBar[i][j][k] - dt/dx*(flux[i][j][k] - flux[i][j-1][k]);
                }
            }
        }

        cout << t << endl;
    }while(t<T);

    u = conservedToPrimitive(u, gamma);
    ofstream output1("rho.dat");
    ofstream output2("vx.dat");
    ofstream output3("vy.dat");
    ofstream output4("p.dat");
    for (int i = 1; i != nxCells + 1; i++){
        double x = x0 + (i-0.5)*dx;
        
            double y = y0 + (i-0.5)*dy;
            output1 << x << " " << u[i][1][1] << endl;
            output2 << x << " " << u[i][1][2] << endl;
            output3 << x << " " << u[i][1][3] << endl;
            output4 << x << " " << u[i][1][4] << endl;


    }

    output1.close();
    output2.close();
    output3.close();
    output4.close();
}



int main()
{
    double x0 = 0, x1 = 1, y0 = 0, y1 = 1, T = 0.25, gamma = 1.4, nxCells = 100, nyCells = 100;
    double dx = (x1 - x0)/nxCells;
    double dy = (y1 - y0)/nyCells;
    mat u0;
    u0.resize(nxCells+2, vector<array<double, 4>> (nyCells+2));
    for (int i = 0; i != nxCells+2; i++){
        double x = x0 + (i-0.5)*dx;
        for (int j = 0; j != nyCells+2; j++){
            double y = y0 + (j-0.5)*dy;
            u0[i][j][0] = ((x < 0.5) ? 1 : 0.125);
            u0[i][j][1] = ((x < 0.5) ? 0 : 0);
            u0[i][j][2] = ((x < 0.5) ? 0 : 0);
            u0[i][j][3] = ((x < 0.5) ? 1 : 0.1);
        }
    }
    solver(u0, T, x0, x1, y0, y1, nxCells, nyCells, gamma, "FORCE");
    return 0;

}