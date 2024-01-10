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
double zeta(double delta_minus, double delta_plus, string limiter);
tuple<mat, mat> dataReconstruct(mat &u, string direction);
array<double, 4> getXFlux(array<double, 4> &uL, array<double, 4> &uR); //u_left and u_right
array<double, 4> getYFlux(array<double, 4> &uD, array<double, 4> &uU); //u_down and u_up
mat split(mat &u, mat &flux, double dt, double &dx, double &dy, int &nxCells, int &nyCells, double &gamma, string direction, string method);
bool inCircle(double x, double y, double R);
void transmissiveBC(mat& u, int nxCells, int nyCells);


bool inCircle(double x, double y, double R){
    return ((x-1)*(x-1) + (y-1)*(y-1) < R*R) ? true : false;
}

void transmissiveBC(mat &u, int nxCells, int nyCells){
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
}

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
    double dt = C*delta/a_max;
    return dt;
}

array<double, 4> f(array<double, 4> &u, double &gamma){
    double rho, vx, vy, p;
    double rhoVx, rhoVy, E;
    tie(rho, rhoVx, rhoVy, E) = tie(u[0], u[1], u[2], u[3]);
    vx = rhoVx/rho;
    vy = rhoVy/rho;
    p = (gamma - 1)*(E - 0.5*(rhoVx*vx + rhoVy*vy));
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
    p = (gamma - 1)*(E - 0.5*(rhoVx*vx + rhoVy*vy));
    array<double, 4> g_res;
    g_res[0] = rhoVy;
    g_res[1] = rhoVx*vy;
    g_res[2] = rhoVy*vy + p;
    g_res[3] = (E + p)*vy;
    return g_res;
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

tuple<mat, mat> dataReconstruct(mat &u, string direction){
    mat uL = u, uR = u;
    int nxCells = u.size()-2;
    int nyCells = u[0].size()-2;
    double w = 0;
    if (direction == "y"){
        for (int j = 1; j != nyCells+1; j++){
            for (int i = 0; i != nxCells+2; i++){
                for (int var = 0; var != 4; var++){
                    double delta_minus = u[i][j][var] - u[i][j-1][var];
                    double delta_plus = u[i][j+1][var] - u[i][j][var];
                    double delta = 0.5*(1+w)*delta_minus + 0.5*(1-w)*delta_plus;
                    uL[i][j][var] = u[i][j][var] - 0.5*zeta(delta_minus, delta_plus)*delta;
                    uR[i][j][var] = u[i][j][var] + 0.5*zeta(delta_minus,delta_plus)*delta;
                }
            }
        }
    }
    else if (direction == "x"){
        for (int i = 1; i != nxCells+1; i++){
            for (int j = 0; j != nyCells+2; j++){
                for (int var = 0; var != 4; var++){
                    double delta_minus = u[i][j][var] - u[i-1][j][var];
                    double delta_plus = u[i+1][j][var] - u[i][j][var];
                    double delta = 0.5*(1+w)*delta_minus + 0.5*(1-w)*delta_plus;
                    uL[i][j][var] = u[i][j][var] - 0.5*zeta(delta_minus, delta_plus)*delta;
                    uR[i][j][var] = u[i][j][var] + 0.5*zeta(delta_minus,delta_plus)*delta;
                }
            }
        }
    }
    return tie(uL, uR);
}

tuple<mat, mat> halfTimeUpdate(mat &uL, mat &uR, double gamma, double dx, double dy, double dt, string direction){
    int nxCells = uL.size()-2;
    int nyCells = uL[0].size()-2;
    if (direction == "x"){
        for (int i = 0; i != nxCells+2; i++){
            for (int j = 0; j != nyCells+2; j++){
                array<double, 4> fL = f(uL[i][j], gamma);
                array<double, 4> fR = f(uR[i][j], gamma);
                for (int var = 0; var != 4; var++){
                    uL[i][j][var] = uL[i][j][var] - 0.5*dt/dx*(fR[var] - fL[var]);
                    uR[i][j][var] = uR[i][j][var] - 0.5*dt/dx*(fR[var] - fL[var]);
                }
            }
        }
    }
    if (direction == "y"){
        for (int i = 0; i != nxCells+2; i++){
            for (int j = 0; j != nyCells+2; j++){
                array<double, 4> gL = g(uL[i][j], gamma);
                array<double, 4> gR = g(uR[i][j], gamma);
                for (int var = 0; var != 4; var++){
                    uL[i][j][var] = uL[i][j][var] - 0.5*dt/dy*(gR[var] - gL[var]);
                    uR[i][j][var] = uR[i][j][var] - 0.5*dt/dy*(gR[var] - gL[var]);
                }
            }
        }
    }
    
    return tie(uL, uR);
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
        RI = g(u_mid, gamma);
        for (int var = 0; var != 4; var++){
            flux[var] = 0.5*(LF[var] + RI[var]);
        }
    }
    
    return flux;
}

mat split(mat &u, mat &flux, double dt, double &dx, double &dy, int &nxCells, int &nyCells, double &gamma, string direction, string method){
    mat uL, uR;
    if (method == "SLIC"){
        tie(uL, uR) = dataReconstruct(u, direction);
        halfTimeUpdate(uL, uR, gamma, dx, dy, dt, direction);
        if (direction == "x"){
            //getXFLux
            for (int i = 0; i != nxCells+1; i++){
                for (int j = 0; j != nyCells+1; j++){
                    flux[i][j] = getXFlux(uR[i][j], uL[i+1][j], gamma, dx, dt, "FORCE");
                }
            }
            for (int i = 1; i != nxCells+1; i++){
                for (int j = 1; j != nyCells+1; j++){
                    for (int k = 0; k != 4; k++){
                        u[i][j][k] = u[i][j][k] - dt/dx*(flux[i][j][k] - flux[i-1][j][k]);
                    }
                }
            }
        }
        else if (direction == "y"){
            //getYFlux
            for (int i = 0; i != nxCells+1; i++){
                for (int j = 0; j != nyCells+1; j++){
                    flux[i][j] = getYFlux(uR[i][j], uL[i][j+1], gamma, dy, dt, "FORCE");
                }
            }

            for (int i = 1; i != nxCells+1; i++){
                for (int j = 1; j != nyCells+1; j++){
                    for (int k = 0; k != 4; k++){
                        u[i][j][k] = u[i][j][k] - dt/dy*(flux[i][j][k] - flux[i][j-1][k]);
                    }
                }
            }
        }
    }
    else {
        if (direction == "x"){
            //getXFLux
            for (int i = 0; i != nxCells+1; i++){
                for (int j = 0; j != nyCells+1; j++){
                    flux[i][j] = getXFlux(u[i][j], u[i+1][j], gamma, dx, dt, method);
                }
            }
            for (int i = 1; i != nxCells+1; i++){
                for (int j = 1; j != nyCells+1; j++){
                    for (int k = 0; k != 4; k++){
                        u[i][j][k] = u[i][j][k] - dt/dx*(flux[i][j][k] - flux[i-1][j][k]);
                    }
                }
            }
        }
        if (direction == "y"){
            //getYFlux
            for (int i = 0; i != nxCells+1; i++){
                for (int j = 0; j != nyCells+1; j++){
                    flux[i][j] = getYFlux(u[i][j], u[i][j+1], gamma, dy, dt, method);
                }
            }

            for (int i = 1; i != nxCells+1; i++){
                for (int j = 1; j != nyCells+1; j++){
                    for (int k = 0; k != 4; k++){
                        u[i][j][k] = u[i][j][k] - dt/dy*(flux[i][j][k] - flux[i][j-1][k]);
                    }
                }
            }
        }
    }
    //Need ensure boundary condition here
    transmissiveBC(u, nxCells, nyCells);
    return u;
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
        transmissiveBC(u, nxCells, nyCells);
        //Dimensional Split 1
        mat u_temp = u;
        split(u_temp, flux, dt, dx, dy, nxCells, nyCells, gamma, "x", method);
        split(u_temp, flux, dt, dx, dy, nxCells, nyCells, gamma, "y", method);
        split(u, flux, dt, dx, dy, nxCells, nyCells, gamma, "y", method);
        split(u, flux, dt, dx, dy, nxCells, nyCells, gamma, "x", method);
        for (int i = 0; i != nxCells + 2; i++){
            for (int j = 0; j != nyCells + 2; j++){
                for (int var = 0; var != 4; var++){
                    u[i][j][var] = 0.5*(u[i][j][var] + u_temp[i][j][var]);
                }
            }
        }

        //Dimensionial Split 2
        // split(u, flux, dt/2, dx, dy, nxCells, nyCells, gamma, "x", method);
        // split(u, flux, dt, dx, dy, nxCells, nyCells, gamma, "y", method);
        // split(u, flux, dt/2, dx, dy, nxCells, nyCells, gamma, "x", method);

        //cout << t << endl;
    }while(t<T);

    u = conservedToPrimitive(u, gamma);
    ofstream output1("rho.dat");
    ofstream output2("vx.dat");
    ofstream output3("vy.dat");
    ofstream output4("p.dat");

    for (int i = 1; i != nxCells + 1; i++){
        double x = x0 + (i-0.5)*dx;
        
        for (int j = 1; j != nyCells + 1; j++){
            double y = y0 + (j-0.5)*dy;
            output1 << x << " " << y << " " << u[i][j][0] << endl;
            output2 << x << " " << y << " " << u[i][j][1] << endl;
            output3 << x << " " << y << " " << u[i][j][2] << endl;
            output4 << x << " " << y << " " << u[i][j][3] << endl;
        }


    }

    output1.close();
    output2.close();
    output3.close();
    output4.close();


    // ofstream output1("rho.dat");
    // ofstream output2("vx.dat");
    // ofstream output3("vy.dat");
    // ofstream output4("p.dat");

    // for (int i = 150; i != nxCells + 1; i++){
    //     double x = x0 + (i-0.5)*dx;
    //     output1 << x << " " << u[i][151][0] << endl;
    //     output2 << x << " " << u[i][151][1] << endl;
    //     output3 << x << " " << u[i][151][2] << endl;
    //     output4 << x << " " << u[i][151][3] << endl;
    // }
    // output1.close();
    // output2.close();
    // output3.close();
    // output4.close();
}



int main()
{
    // double x0 = 0, x1 = 1, y0 = 0, y1 = 1, T = 0.25, gamma = 1.4, nxCells = 200, nyCells = 200;
    // double dx = (x1 - x0)/nxCells;
    // double dy = (y1 - y0)/nyCells;
    // mat u0;
    // u0.resize(nxCells+2, vector<array<double, 4>> (nyCells+2));
    // for (int i = 0; i != nxCells+2; i++){
    //     double x = x0 + (i-0.5)*dx;
    //     for (int j = 0; j != nyCells+2; j++){
    //         double y = y0 + (j-0.5)*dy;
    //         u0[i][j][0] = ((x < 0.5) ? 1 : 0.125);
    //         u0[i][j][1] = ((x < 0.5) ? 0 : 0);
    //         u0[i][j][2] = ((x < 0.5) ? 0 : 0);
    //         u0[i][j][3] = ((x < 0.5) ? 1 : 0.1);
    //     }
    // }


    //Diagonal Initial Condition
    // for (int i = 0; i != nxCells+2; i++){
    //     double x = x0 + (i-0.5)*dx;
    //     for (int j = 0; j != nyCells+2; j++){
    //         double y = y0 + (j-0.5)*dy;
    //         u0[i][j][0] = ((x-y < 0) ? 1 : 0.125);
    //         u0[i][j][1] = ((x-y < 0) ? 0 : 0);
    //         u0[i][j][2] = ((x-y < 0) ? 0 : 0);
    //         u0[i][j][3] = ((x-y < 0) ? 1 : 0.1);
    //     }
    // }

    // Cylindrical Explosion
    double x0 = 0, x1 = 2, y0 = 0, y1 = 2, T = 0.25, gamma = 1.4, nxCells = 300, nyCells = 300;
    double dx = (x1 - x0)/nxCells;
    double dy = (y1 - y0)/nyCells;
    mat u0;
    u0.resize(nxCells+2, vector<array<double, 4>> (nyCells+2));
    for (int i = 0; i != nxCells+2; i++){
        double x = x0 + (i-0.5)*dx;
        for (int j = 0; j != nyCells+2; j++){
            double y = y0 + (j-0.5)*dy;
            u0[i][j][0] = (inCircle(x, y, 0.4) ? 1 : 0.125);
            u0[i][j][1] = (inCircle(x, y, 0.4) ? 0 : 0);
            u0[i][j][2] = (inCircle(x, y, 0.4) ? 0 : 0);
            u0[i][j][3] = (inCircle(x, y, 0.4) ? 1 : 0.1);
        }
    }


    //Liu and Lax
    // double x0 = 0, x1 = 1, y0 = 0, y1 = 1, T = 0.2, gamma = 1.4, nxCells = 400, nyCells = 400;
    // double dx = (x1 - x0)/nxCells;
    // double dy = (y1 - y0)/nyCells;
    // mat u0;
    // u0.resize(nxCells+2, vector<array<double, 4>> (nyCells+2));
    // // Test 1
    // for (int i = 0; i != nxCells+2; i++){
    //     double x = x0 + (i-0.5)*dx;
    //     for (int j = 0; j != nyCells+2; j++){
    //         double y = y0 + (j-0.5)*dy;
    //         if ((x > 0.5) && (y > 0.5))
    //             u0[i][j] = array<double, 4> {0.5313, 0, 0, 0.4};
            
    //         else if ((x <= 0.5) && (y > 0.5))
    //             u0[i][j] = array<double, 4> {1, 0.7276, 0, 1};
            
    //         else if ((x <= 0.5) && (y <= 0.5))
    //             u0[i][j] = array<double, 4> {0.8, 0, 0, 1};
            
    //         else
    //             u0[i][j] = array<double, 4> {1, 0, 0.7276, 1};
            
    //     }
    // }

    //Test 2
    // for (int i = 0; i != nxCells+2; i++){
    //     double x = x0 + (i-0.5)*dx;
    //     for (int j = 0; j != nyCells+2; j++){
    //         double y = y0 + (j-0.5)*dy;
    //         if ((x > 0.5) && (y > 0.5))
    //             u0[i][j] = array<double, 4> {1, 0.75, -0.5, 1};
            
    //         else if ((x <= 0.5) && (y > 0.5))
    //             u0[i][j] = array<double, 4> {2, 0.75, 0.5, 1};
            
    //         else if ((x <= 0.5) && (y <= 0.5))
    //             u0[i][j] = array<double, 4> {1, -0.75, 0.5, 1};
            
    //         else if ((x > 0.5) && (y <= 0.5))
    //             u0[i][j] = array<double, 4> {3, -0.75, -0.5, 1};
            
    //     }
    // }



    solver(u0, T, x0, x1, y0, y1, nxCells, nyCells, gamma, "SLIC");
    return 0;

}