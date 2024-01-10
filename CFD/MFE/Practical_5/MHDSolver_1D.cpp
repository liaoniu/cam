#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;
using mat = vector<array<double, 7>>;   // u[i][var] = var[i]


array<double, 7> indi_primitiveToConserved(const array<double, 7> &u_pri, double &gamma, double Bx);
array<double, 7> indi_conservedToPrimitive(const array<double, 7> &u_con, double &gamma, double Bx);
mat conservedToPrimitive(mat &u, double &gamma, double Bx);
mat primitiveToConserved(mat &u, double &gamma, double Bx);
double calc_dt(mat &u_con, double &gamma, double &dx, double &dy, double Bx);
array<double, 7> f(array<double, 7> &u, double &gamma, double Bx);
array<double, 7> s(array<double, 7> &u, double r, double gamma, int alpha);
double zeta(double delta_minus, double delta_plus, string limiter);
tuple<mat, mat> dataReconstruct(mat &u, string direction);
array<double, 7> getXFlux(array<double, 7> &uL, array<double, 7> &uR); //u_left and u_right
mat split(mat &u, mat &flux, double dt, double &dx, int &nxCells, double &gamma, string method);
mat operatorSplit(mat &u, double &dt, double &r0, double dr, int &nxCells, double &gamma, string method);
bool inCircle(double x, double y, double R);
void transmissiveBC(mat& u, int nxCells);
void halfReflectiveBC(mat&u, int nxCells);

bool inCircle(double x, double y, double R){
    return ((x-1)*(x-1) + (y-1)*(y-1) < R) ? true : false;
}

void transmissiveBC(mat &u, int nxCells){
    u[0] = u[1];
    u[nxCells+1] = u[nxCells];
}

void halfReflectiveBC(mat &u, int nxCells){
    // u[0] = u[3];
    // u[1] = u[2];
    // u[nxCells+2] = u[nxCells+1];
    // u[0][1] = -u[0][1];
    // u[1][1] = -u[1][1];
}

array<double, 7> indi_primitiveToConserved(const array<double, 7> &u_pri, double &gamma, double Bx=0){
    array<double, 7> u_con = u_pri;
    double rho, vx, vy, vz, p, By, Bz;
    tie(rho, vx, vy, vz, p, By, Bz) = tie(u_pri[0], u_pri[1], u_pri[2], u_pri[3], u_pri[4], u_pri[5], u_pri[6]);
    double rhoVx = rho*vx, rhoVy = rho*vy, rhoVz = rho*vz;
    double epsilon = p/(gamma - 1)/rho;
    double U = rho*epsilon + rho*(vx*vx + vy*vy + vz*vz)/2 + (Bx*Bx + By*By + Bz*Bz)/2;
    tie(u_con[0], u_con[1], u_con[2], u_con[3], u_con[4], u_con[5], u_con[6]) = tie(rho, rhoVx, rhoVy, rhoVz, U, By, Bz);
    return u_con;
}

array<double, 7> indi_conservedToPrimitive(const array<double, 7> &u_con, double &gamma, double Bx=0){
    array<double, 7> u_pri = u_con;
    double rho, rhoVx, rhoVy, rhoVz, U, By, Bz;
    tie(rho, rhoVx, rhoVy, rhoVz, U, By, Bz) = tie(u_con[0], u_con[1], u_con[2], u_con[3], u_con[4], u_con[5], u_con[6]);
    double vx = rhoVx/rho, vy = rhoVy/rho, vz = rhoVz/rho;
    double epsilon = (U - (Bx*Bx + By*By + Bz*Bz)/2 - (vx*vx + vy*vy + vz*vz)*rho/2)/rho;
    double p = (gamma-1)*rho*epsilon;
    tie(u_pri[0], u_pri[1], u_pri[2], u_pri[3], u_pri[4], u_pri[5], u_pri[6]) = tie(rho, vx, vy, vz, p, By, Bz);
    return u_pri;
}

mat conservedToPrimitive(mat &u, double &gamma, double Bx=0){
    //u would be changed
    for (int i = 0; i != u.size(); i++){
        u[i] = indi_conservedToPrimitive(u[i], gamma, Bx);
    }
    return u;
}

mat primitiveToConserved(mat &u, double &gamma, double Bx=0){
    //u would be changed

    for (int i = 0; i != u.size(); i++){
        u[i] = indi_primitiveToConserved(u[i], gamma, Bx);
    }
    
    return u;
}

double calc_dt(mat &u_con, double &gamma, double &dx, double Bx=0){
    double delta = dx;
    double C = 0.8;
    double a_max = 0;
    int m = u_con.size();
    for (int i = 0; i != m; i++){
            double rho, vx, vy, vz, p, By, Bz;
            array<double, 7> u_pri_indi = indi_conservedToPrimitive(u_con[i], gamma, Bx);
            tie(rho, vx, vy, vz, p, By, Bz) = tie(u_pri_indi[0], u_pri_indi[1], u_pri_indi[2], u_pri_indi[3], u_pri_indi[4], u_pri_indi[5], u_pri_indi[6]);
            double ca = fabs(Bx)/sqrt(rho);
            double cs = sqrt(gamma*p/rho);
            double cf = sqrt(0.5*(cs*cs + ca*ca + sqrt(pow(cs*cs + ca*ca, 2) - 4*cs*cs*Bx*Bx/rho)));
            double a = sqrt(vx*vx + vy*vy + vz*vz) + cf;
            a_max = ((a_max < a) ? a : a_max);
    }
    double dt = C*delta/a_max;
    return dt;
}

array<double, 7> f(array<double, 7> &u, double &gamma, double Bx=0){
    double rho, vx, vy, vz, p, By, Bz;
    double rhoVx, rhoVy, rhoVz, U;
    tie(rho, rhoVx, rhoVy, rhoVz, U, By, Bz) = tie(u[0], u[1], u[2], u[3], u[4], u[5], u[6]);
    double B = sqrt(Bx*Bx + By*By + Bz*Bz);
    array<double, 7> u_con = indi_conservedToPrimitive(u, gamma, Bx);
    tie(ignore, vx, vy, vz, p, ignore, ignore) = tie(u_con[0], u_con[1], u_con[2], u_con[3], u_con[4], u_con[5], u_con[6]);
    array<double, 7> f_res;
    f_res[0] = rhoVx;
    f_res[1] = rhoVx*vx + p + B*B/2 - Bx*Bx;
    f_res[2] = rhoVx*vy - Bx*By;
    f_res[3] = rhoVx*vz - Bx*Bz;
    f_res[4] = (U + p + B*B/2)*vx - (vx*Bx + vy*By + vz*Bz)*Bx;
    f_res[5] = By*vx - Bx*vy;
    f_res[6] = Bz*vx - Bx*vz;
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
        for (int var = 0; var != 7; var++){
            double delta_minus = u[i][var] - u[i-1][var];
            double delta_plus = u[i+1][var] - u[i][var];
            double delta = 0.5*(1+w)*delta_minus + 0.5*(1-w)*delta_plus;
            uL[i][var] = u[i][var] - 0.5*zeta(delta_minus, delta_plus)*delta;
            uR[i][var] = u[i][var] + 0.5*zeta(delta_minus, delta_plus)*delta;
        }
    }
    return tie(uL, uR);
}

tuple<mat, mat> halfTimeUpdate(mat &uL, mat &uR, double gamma, double dx, double dt, double Bx=0){
    int nxCells = uL.size()-2;
    for (int i = 0; i != nxCells+2; i++){
        array<double, 7> fL = f(uL[i], gamma, Bx);
        array<double, 7> fR = f(uR[i], gamma, Bx);
        for (int var = 0; var != 7; var++){
            uL[i][var] = uL[i][var] - 0.5*dt/dx*(fR[var] - fL[var]);
            uR[i][var] = uR[i][var] - 0.5*dt/dx*(fR[var] - fL[var]);
        }
        
    }
    return tie(uL, uR);
}


array<double, 7> getXFlux(array<double, 7> &uL, array<double, 7> &uR, double &gamma, double &dx, double &dt, string method, double Bx=0){ //u_left and u_right
    array<double, 7> flux, u_mid;
    if (method == "FORCE"){
        array<double, 7> fL, fR, LF, RI;
        fL = f(uL, gamma, Bx);
        fR = f(uR, gamma, Bx);
        for (int var = 0; var != 7; var++){
            LF[var] = 0.5*dx/dt*(uL[var] - uR[var]) + 0.5*(fL[var] + fR[var]);
            u_mid[var] = 0.5*(uL[var] + uR[var]) - 0.5*dt/dx*(fR[var] - fL[var]);
        }
        RI = f(u_mid, gamma);
        for (int var = 0; var != 7; var++){
            flux[var] = 0.5*(LF[var] + RI[var]);
        }

    }
    return flux;
}



mat split(mat &u, mat &flux, double dt, double &dx, int &nxCells, double &gamma, string method, double Bx=0){
    mat uL, uR;
    if (method == "SLIC"){
        tie(uL, uR) = dataReconstruct(u);
        halfTimeUpdate(uL, uR, gamma, dx, dt);
        //getXFLux
        for (int i = 0; i != nxCells+1; i++){
                flux[i] = getXFlux(uR[i], uL[i+1], gamma, dx, dt, "FORCE", Bx);
        }
        for (int i = 1; i != nxCells+1; i++){
            for (int k = 0; k != 7; k++){
                u[i][k] = u[i][k] - dt/dx*(flux[i][k] - flux[i-1][k]);
            }
        }
    }
    else {
        //getXFLux
        for (int i = 0; i != nxCells+1; i++){
                flux[i] = getXFlux(u[i], u[i+1], gamma, dx, dt, method, Bx);
        }
        for (int i = 1; i != nxCells+1; i++){
            for (int k = 0; k != 7; k++){
                u[i][k] = u[i][k] - dt/dx*(flux[i][k] - flux[i-1][k]);
            }

        }
    }

    //Need ensure boundary condition here
    halfReflectiveBC(u, nxCells);
    return u;
}


void solver(mat &u0, double T, double x0, double x1, int nxCells = 100, double gamma = 1.4, string method = "FORCE", double Bx=0){
    double dx = (x1-x0)/nxCells;
    double t = 0;
    double dt = 0;
    mat flux;
    mat u = u0;
    flux.resize(nxCells+1);
    u = primitiveToConserved(u, gamma, Bx);

    do{
        dt = calc_dt(u, gamma, dx, Bx);
        t += dt;
        //Transmissive BC
        transmissiveBC(u, nxCells);
        split(u, flux, dt, dx, nxCells, gamma, method, Bx);
        
    }while(t<T);

    u = conservedToPrimitive(u, gamma, Bx);
    ofstream output1("rho.dat");
    ofstream output2("vx.dat");
    ofstream output3("vy.dat");
    ofstream output4("vz.dat");
    ofstream output5("p.dat");
    ofstream output6("By.dat");
    ofstream output7("Bz.dat");

    for (int i = 1; i != nxCells + 1; i++){
        double x = x0 + (i-0.5)*dx;
        output1 << x << " " << u[i][0] << endl;
        output2 << x << " " << u[i][1] << endl;
        output3 << x << " " << u[i][2] << endl;
        output4 << x << " " << u[i][3] << endl;
        output5 << x << " " << u[i][4] << endl;
        output6 << x << " " << u[i][5] << endl;
        output7 << x << " " << u[i][6] << endl;

        }

    output1.close();
    output2.close();
    output3.close();

}



int main()
{
    double x0 = 0, x1 = 800, T = 80, gamma = 2, nxCells = 1000, nyCells = 1000, Bx = 0.75;
    double dx = (x1 - x0)/nxCells;
    mat u0;
    u0.resize(nxCells+2);
    for (int i = 0; i != nxCells+2; i++){
        double x = x0 + (i-0.5)*dx;
        u0[i][0] = ((x < 400) ? 1 : 0.125);
        u0[i][1] = ((x < 400) ? 0 : 0);
        u0[i][2] = ((x < 400) ? 1 : 0.2);
        u0[i][3] = ((x < 400) ? 0 : 0);
        u0[i][4] = ((x < 400) ? 1 : 0.1);
        u0[i][5] = ((x < 400) ? 1 : -1);
        u0[i][6] = ((x < 400) ? 0 : 0);

    }

    solver(u0, T, x0, x1, nxCells, gamma, "SLIC", Bx);
    return 0;

}